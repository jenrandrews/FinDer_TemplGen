# stdlib imports
import os
import fnmatch
from math import floor, ceil
from numpy import sqrt
import geographiclib.geodesic as geo
import numpy as np
import pandas as pd
import create_site_mesh as csm

def rdRuptureParams(fname):
    '''
    Get centroid and magnitude from rupture file
    '''
    lats = []
    lons = []
    def getVal(instr):
        return float(instr.split('>')[1].split('<')[0])
    with open(fname, 'r') as fin:
        for l in fin:
            if l.find('magnitude') > -1:
                mag = getVal(l)
                continue
            if l.find('hypocenter') > -1:
                lon = float(l.split('lon="')[1].split('"')[0])
                lat = float(l.split('lat="')[1].split('"')[0])
                continue
            if l.find('<') == -1:
                fs = l.split()
                lons.append(float(fs[0]))
                lats.append(float(fs[1]))
    return mag, lat, lon, lats, lons

if __name__ == '__main__':
    rupdir = 'event_xmls'
    rlist = fnmatch.filter(os.listdir(rupdir), 'hik_*.xml')

    bounds = [-47.5, -34.0, 166.0, 179.0]
    absminlat = None
    absmaxlat = None
    absminlon = None
    absmaxlon = None

    maxdist = 200.*1000.

    rup_extents = {}
    for rup in rlist:
        mag, lat, lon, lats, lons = rdRuptureParams(os.path.join(rupdir, rup))

        maxlat = geo.Geodesic.WGS84.Direct(max(lats), max(lons), 0., maxdist)['lat2']
        maxlon = geo.Geodesic.WGS84.Direct(max(lats), max(lons), 90., maxdist)['lon2']
        minlat = geo.Geodesic.WGS84.Direct(min(lats), min(lons), 180., maxdist)['lat2']
        minlon = geo.Geodesic.WGS84.Direct(min(lats), min(lons), 270., maxdist)['lon2']
        if maxlon < 0:
            maxlon += 360.
        rup_extents[rup] = [minlat, maxlat, minlon, maxlon]

        if absminlat is None or minlat < absminlat:
            absminlat = minlat
        if absmaxlat is None or maxlat > absmaxlat:
            absmaxlat = maxlat
        if absminlon is None or minlon < absminlon:
            absminlon = minlon
        if absmaxlon is None or maxlon > absmaxlon:
            absmaxlon = maxlon

    templ_sets = {}
    rups_assigned = []
    absminlon = floor(absminlon)
    absminlat = floor(absminlat)
    absmaxlon = ceil(absmaxlon)
    absmaxlat = ceil(absmaxlat)
    max_size = max([absmaxlon - absminlon, absmaxlat - absminlat])
    for size in range(4, max_size+2, 2):
        for start_lat in range(absminlat, absminlat + max_size - size + 2, 2):
            for start_lon in range(absminlon, absminlon + max_size - size + 2, 2):
                tkey = f'{size}_{start_lat}_{start_lon}'
                if tkey not in templ_sets:
                    templ_sets[tkey] = []
                for rup in rlist:
                    if rup in rups_assigned:
                        continue
                    if rup_extents[rup][0] >= start_lat and \
                       rup_extents[rup][1] <= start_lat + size and \
                       rup_extents[rup][2] >= start_lon and \
                       rup_extents[rup][3] <= start_lon + size:
                           templ_sets[tkey].append(rup)
                           rups_assigned.append(rup)

    for tset in templ_sets:
        if len(templ_sets[tset]) == 0:
            continue
        size = float(tset.split('_')[0])
        minlat = float(tset.split('_')[1])
        minlon = float(tset.split('_')[2])
        maxlat = minlat + size
        maxlon = minlon + size

        # Apply bounds
        if minlat < bounds[0]:
            minlat = bounds[0]
        if maxlat < bounds[1]:
            maxlat = bounds[1]
        if minlon < bounds[2]:
            minlon = bounds[2]
        if maxlon < bounds[3]:
            maxlon = bounds[3]
       
        lats, lons = csm.mesh_from_bb(5., minlat, maxlat, minlon, maxlon)
        csm.vs30_mesh(tset, lats, lons)
        print(tset, len(templ_sets[tset]), templ_sets[tset])
