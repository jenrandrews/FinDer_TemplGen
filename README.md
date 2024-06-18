# FinDer_TemplGen
FinDer Template Generation scripts

# Description

These are simple scripts to create template input files for FinDer. Information on FinDer can be found here: 
 * http://www.seismo.ethz.ch/en/research-and-teaching/products-software/EEW/finite-fault-rupture-detector-finder/

# Requirements

The template generation script is built on top of ShakeMap, which in turn uses OpenQuake. Documentation for ShakeMap can be found here:
 * https://usgs.github.io/shakemap/index.html 
 * https://github.com/usgs/shakemap/wiki

For ShakeMap versions 4.2 and greater, the easiest approach is to create a conda environment and
install esi-shakelib:
 * https://pypi.org/project/esi-shakelib/

2024/06/17: create a python 3.9 environment and pip install esi-shakelib.

# Usage

The easiest way to run the scripts is to leverage your shakelib (if you did a pip install into
a conda environment) or shakemap environment:
```
conda run -n <your shake env> python <FinDer_TemplGen script choice>.py -e event.conf -g gmpe.conf -c calc.conf
```
Example configuration files are provided and should be edited before running the script. 

# Available Scripts for Computing Ground Motion at Position/Distance
There are four wrapper scripts available designed to create FinDer templates or FinDer scenario input data (data_ files). The different wrappers allow for different input specifications (fault geometry or parameters) and/or different outputs.
 * **makeFinDerTemplates**: designed to create generic FinDer templates
 * **makeFaultSpecificFinDerTemplates**: designed to create fault-specific FinDer template sets
 * **makeEventData**: designed to create templates or scenario data for a single fault rupture 
 * **makeScenarioStnData**: designed to create ground motion data for generic fault scenarios

## makeFinDerTemplates
When using the FinDer template wrapper (makeFinDerTemplates.py):
 * In calc.conf 
   * Set the configuration parameter ["grid"]["compute"] to True.
   * Supply the ["grid"]["griddkm"] value, i.e. the grid spacing in km.
   * Supply the fixed ["grid"]["vs30"] to be used. 
   * If asymmetric templates are required *as well as* symmetric templates, set ["grid"]["asym"] to True. 
   * The additional info files such as template_info, as used by FinDer v3, can be created by setting the "rupinfo" configuration option to True.
 * In event.conf
   * Set event parameters with strike as 0.

The script works internally with (lat, lon) values, as required by ShakeMap's structures, so the fault is placed approximately at the equator to make geographic to distance conversions simpler. The script will always create symmetric templates, and optionally *also* asymmetric templates if the "asym" configuration option is set to True. 

## makeFaultSpecificFinDerTemplates
When using the fault-specific FinDer template wrapper (makeFaultSpecificFinDerTemplates.py):
 * In calc.conf 
   * Set the configuration parameter ["grid"]["compute"] to True.
   * Supply the ["grid"]["griddkm"] value, i.e. the grid spacing in km.
   * Supply the fixed ["grid"]["vs30"] to be used. 
   * The ["fault-specific"]["centroid_polygon_distance"] parameter is used to create the centroid polygon for the template set file, and is a radius in km. 
   * The scripts switch between geographic and cartesian coordinates to perform the walk along fault to create sub-faults in the required magnitude range. That switch requires the user to set the epsg code for the cartesian transformation using ["fault-specific"]["epsg"]. 
   * The additional info files such as latitude, longitude, template_info and template set, as used by FinDer v3, can be created by setting the "rupinfo" configuration option to True.
   * Template sets are used to group the templates for efficiency when running FinDer. Grouping is based on configuration parameters ["fault-specific"]["maxperset"] and ["fault-specific"]["minperset"]. The former sets the target maximum number of smallest magnitude templates in a template set, and the latter sets the target minimum number of largest magnitude templates in a set (unless the total number for that magnitude is smaller). Setting a very large "maxperset" will therefore result in a single set for all templates. 
   * If you wish to apply a station mask to the template, supply a ["points"]["points_file"] option, and supply a ["points"]["stnmaskdist"] distance (in km). This will create a mask allowing only grid points within 'stnmaskdist' of a supplied point in the 'points_file' to have a non-zero value. In this case the 'points_file' should have space delimited fields: lan lon dummy dummy at minimum since the vs30 and name values used elsehwere are not needed.
 * In event.conf
   * The rupture is specified via a geojson file to define the surface trace using ["evmech"]["geometry"]. N.B. The file *must* have a .geojson file extension. The fault plane is constructed using the ["evmech"]["dip"] and the ["evmech"]["strike"] is used to determine the dip direction. When creating sub-faults along the surface trace, the overlap between centroids is currently fixed at min([10% fault length, 10km]). 

## makeEventData
When creating files for a specific rupture (makeEventData.py):
 * In event.conf
   * Set ["evmech"]["geometry"] as a rupture.json file path or directory path with a number of rupture_XX.json files. The rupture.json files should follow the ShakeMap format and *must* have a .json file extension. Note that the input is currently only expecting a single polygon description with the points 1:N/2 specifying the top edge vertices and points N/2:N-1 specifying the bottom edge vertices, which are used as edges in the constructor of a ComplexFaultSurface. Magnitude and hypocenter location are also read from the rupture.json, so mag and centroid ranges in calc.conf are unused. This could be extended in future. 
 * In calc.conf
   * Specify either ["grid"] or ["points"] to True to create either template or scenario files.

An additional calculation option is available for points at given distances from a rupture by using ["pt"]["compute"] set to True in calc.conf. The array of Joyner-Boore distances (Rjb in km) are specified using ["pt"]["rjblogmin"], ["pt"]["rjblogmax"] and ["pt"]["rjblogstep"] in calc.conf, and other distance measures are assumed: Rx = Rjb and Ry0 = 0. A fixed vs30 is specified by ["pt"]["vs30"] in calc.conf.

## makeScenarioStnData
When using the scenario wrapper (makeScenarioStnData.py) follow the instructions for
makeFinDerTemplate in setting up generic rupture parameters. Additionally set the following:
 * In calc.conf
   * Set the parameter ["points"]["compute"] to True.
   * Set the path for a points file with parameter ["points"]["points_file"]. The file should contain lines with space delimited fields: lat lon vs30 name. The first three columns will be used as input to the ground motion computation. The name will be used when generating a FinDer data_ file.
 * Modify makeScenarioStnData.py directly to set centroid values for your scenarios. This script is
   not yet generalised to use only configuration file options and will likely be further developed!

An additional calculation option is available for points at given distances from a rupture by using ["pt"]["compute"] set to True in calc.conf. The array of Joyner-Boore distances (Rjb in km) are specified using ["pt"]["rjblogmin"], ["pt"]["rjblogmax"] and ["pt"]["rjblogstep"] in calc.conf, and other distance measures are assumed: Rx = Rjb and Ry0 = 0. A fixed vs30 is specified by ["pt"]["vs30"] in calc.conf.



# Scaling Relations
Currently available scaling relations:
 * Leonard2014_Interplate
 * WC1994 (i.e. Wells & Coppersmith (1994))
 * Strasser2010_Interface
 * Strasser2010_Intraslab
 * Blaser2010
 * Skarlatoudis2016


Other scaling relations used in ShakeMap/OQ have not been extended to have magnitude-length and magnitude-width relations. Note that rake information may be used to select the faulting style relation.

# GMPEs
The setup of GMPEs uses the same concept of weighted sets as used in ShakeMap, but in gmpe.conf the user explicitly creates the set and weights (i.e. tectonic weighting is not altered with location/depth as the template set is not designed to be geographically positioned).

# Earthquake Depth
Configure ["evloc"]["hypo_depth"] in event.conf as the value to be used as centroid in the rupture plane. Once fault width exceeds this depth and/or the seismogenic depth, the centroid depth may be altered within the script, with a note written to log.

# Rupinfo.tbl
The rupinfo.tbl file that is optionally created is not needed by FinDer, but provides values for the input hypocentral depth and dip used, and for each magnitude it provides fault plane length, fault plane width, fault plane top depth and fault plane bottom depth.
