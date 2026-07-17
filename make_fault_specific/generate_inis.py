import glob
import os.path

from xml.etree import ElementTree as ET
ET.register_namespace("", "http://openquake.org/xmlns/nrml/0.5")

import configparser

def rdTemplateSetsFile(fname):
    '''
    Read the template set file
    '''
    tsets = {}
    with open(fname, 'r') as fin:
        for l in fin:
            fs = l.split()
            tfiles = [x.replace('[','').replace(']','').replace(',','').replace("'",'') for x in fs[2:]]
            tsets[fs[0]] = tfiles
    return tsets

def write_inversion_ini(idir: str, odir: str, name_i: str, file_dict: dict):
    '''
    This script generates the ini files for use in OQ, by reading and adapting the template ini files.
    '''
    inst = file_dict['inst']
    description = (f"Hazard maps for {inst} component of {name_i}")
    config = configparser.ConfigParser()
    config.read(f"{os.path.dirname(os.path.abspath(__file__))}/logic_tree_template.ini")
    config["general"]["description"] = description
    config["general"]["calculation_mode"] = "scenario"
    if inst == 'crustal':
        config["general"]["random_seed"] = "1234"
    else:
        config["general"]["random_seed"] = "12345"
    config["logic_tree"]["number_of_logic_tree_samples"] = "0"
    config["rupture"]["rupture_model_file"] = f"../{idir}/{name_i}.xml"
    config["rupture"]["rupture_mesh_spacing"] = "2.0"
    config["site_params"]["site_model_file"] = file_dict['site_model']
    config["hazard_calculation"]["gsim_logic_tree_file"] = file_dict['gsim_logic_tree_file']
    config["hazard_calculation"]["intensity_measure_types"] = "PGA"
    config["hazard_calculation"]["maximum_distance"] = "2000.0"
    config["hazard_calculation"]["truncation_level"] = "0"
    config["hazard_calculation"]["number_of_ground_motion_fields"] = "1"
    config["output"]["export_dir"] = "outputs"
    config.write(open(f"{odir}/{name_i}.ini", "w"))
    return

if __name__ == "__main__":
    odir = 'event_inis'
    idir = 'event_xmls'

    if not os.path.isdir(odir):
        os.makedirs(odir)

    # Fault-specific templates
    file_dict = {
            'inst': 'subduction',
            'gsim_logic_tree_file': '../../NZ_NSHM_GMM_av_subduction.xml',
    }


    #### Create the OpenQuake ini files
    # Read template sets information to set the site meshes
    tsets = rdTemplateSetsFile('tsets.dat')
    for xml in glob.glob(f'{idir}/*.xml'):
        name = os.path.basename(xml).split(".xml")[0]
        tset = None
        for t in tsets:
            if os.path.basename(xml) in tsets[t]:
                tset = t
                break
        if tset is None:
            print(f'Could not find the template set for {xml}')
            continue
        file_dict['site_model'] = f'../site_meshes/{tset}_sites.csv'
        write_inversion_ini(idir, odir, name, file_dict)

