# -*- coding: utf-8 -*-
import re
import warnings
from collections import namedtuple
from pathlib import Path

"""GOCAD Object file template:

https://gist.github.com/T4mmi/b0545c0dfd7b60f3e4f15f3b4e54f7e8

GOCAD <type> <version>
HEADER {
name: <name>
[<key>: <value>]
}
[PROPERTIES <name> ... <name>]
ATOM <ID> <X> <Y> <Z> [<PV> ...]
[<SUBSET_TYPE>]
[<CELL_TYPE> <ATOM> ... <ATOM>]
END

based on : http://paulbourke.net/dataformats/gocad/gocad.pdf
"""

# regex that matches objects "GOCAD ... END"
BLOCK = re.compile(r"((?<=^GOCAD ).+?(?=\s*^END\s*$))+", re.MULTILINE | re.DOTALL)


class AsciiReader:
    """A static class helper to read GOCAD ASCII objects files"""

    Object = namedtuple("GocadObject", ["info", "atoms", "cells", "properties"])
    Formats = ("python", "meshio", "vtk", "wkt")
    Objects = ("TSurf", "TSolid", "VSet")  # missing Well, PLine

    @staticmethod
    def open(filename: str, safe_mod: bool = True, throw_on_error: bool = True):
        """check/open file and split objects"""

        obj = Path(filename)
        raw = obj.read_text()
        try:
            if safe_mod:
                blocks = BLOCK.findall(raw)
            else:
                blocks = []
                _, _, raw = raw.partition("GOCAD ")
                while raw:
                    block, _, raw = raw.partition("GOCAD ")
                    blocks.append(block)
        except Exception as err:
            if throw_on_error:
                raise err
            else:
                return None
        else:
            return blocks

    @staticmethod
    def parse(block: str, throw_on_error: bool = True):
        """parse a single object"""

        try:
            obj = AsciiReader.decode(block)

        except Exception as err:
            if throw_on_error:
                raise err
            else:
                return None

        else:
            return obj

    @staticmethod
    def read(filename: str, safe_mod: bool = True, use_parallel: bool = True):
        """entry point for complete reading"""

        blocks = AsciiReader.open(filename, safe_mod)

        if use_parallel:
            from multiprocessing.pool import ThreadPool

            pool = ThreadPool(processes=len(blocks))
            blocks = pool.map(AsciiReader.parse, blocks)
        else:
            blocks = [AsciiReader.parse(block) for block in blocks]

        return blocks

    @staticmethod
    def decode(block: str, throw_on_error: bool = False):
        lines, l = block.splitlines(), 0

        line, l = lines[l].strip(), l + 1
        obj_type, *obj_version = line.split()

        if obj_type not in AsciiReader.Objects:
            if throw_on_error:
                raise NotImplementedError(f"Unsupported GOCAD object type: {obj_type}")
            else:
                warnings.warn(
                    f"Ignoring unsupported GOCAD object type: {obj_type}",
                    RuntimeWarning,
                )

        # header/metadata map
        info = {
            "type": obj_type,
            "version": 0,
        }
        # geometry parts
        atoms, cells = [], []
        # properties
        len_prop, tag_prop, nan_prop, properties = 0, (""), None, []

        while l < len(lines):

            line, l = lines[l].strip(), l + 1
            if not line or line.startswith("#"):
                continue
            elif line.startswith("HEADER {"):
                line, l = lines[l].strip(), l + 1
                while l < len(lines) and "}" not in line:
                    if not line or line.startswith("#"):
                        continue
                    if line.count(":") > 1:
                        lsplit = line.split(":")
                        key = lsplit[0]
                        value = ":".join(lsplit[1:])
                    else:
                        key, value = line.split(":")
                    info[key.strip("*")] = value.strip()
                    line, l = lines[l].strip(), l + 1

            el = line.split()

            if el[0] in ("END"):
                break  # used for unsafe reading

            elif el[0] in ("VRTX", "PVRTX"):
                atoms.append(tuple(float(f) for f in el[2:5]))
                if el[0].startswith("P"):
                    properties.append(tuple(f for f in el[5 : 5 + len_prop]))
                elif nan_prop:
                    properties.append(nan_prop)

            elif el[0] in ("ATOM", "PATOM"):
                atoms.append(atoms[int(el[2]) - 1])
                if el[0].startswith("P"):
                    properties.append(tuple(f for f in el[3 : 3 + len_prop]))
                elif nan_prop:
                    properties.append(nan_prop)

            elif el[0] in ("SEG", "TRGL", "TETRA"):
                cells.append(tuple(int(i) - 1 for i in el[1:]))

            elif el[0] in ("HDR"):  # single-line style header metadata
                key, value = "".join(el[1:]).split(":")
                info[key.strip("*")] = value.strip()

            elif el[0] in ("PROPERTIES"):
                assert not tag_prop  # there should be only one properties line
                assert not atoms  # properties should be declared prior to atoms
                tag_prop = tuple(el[1:])
                len_prop = len(tag_prop)
                tuple([None] * len_prop)

            elif el[0] in ("NO_DATA_VALUES"):
                nan_prop = tuple(el[1:])
                assert len(nan_prop) == len_prop

        if tag_prop:
            info["properties"] = {
                "names": tag_prop,
                "default": nan_prop,
            }

        return AsciiReader.Object(info, atoms, cells, list(zip(*properties)))
