# FinDer_TemplGen
FinDer Template Generation scripts

# Description

These are simple scripts to create template input files for FinDer. Information on FinDer can be found here: 
 * http://www.seismo.ethz.ch/en/research-and-teaching/products-software/EEW/finite-fault-rupture-detector-finder/

# Requirements

The template generation script is built on top of ShakeMap, which in turn uses OpenQuake. Documentation for ShakeMap can be found here:
 * https://usgs.github.io/shakemap/index.html 
 * https://github.com/usgs/shakemap/wiki
 
# Usage

The easiest way to run the script is to leverage the shakemap environment:
```
conda run -n shakemap python makeFinDerTemplates.py -e event.conf -g gmpe.conf -c calc.conf
```
```
conda run -n shakemap python makeScenarioStnData.py -e event.conf -g gmpe.conf -c calc.conf
```
Example configuration files are provided and should be edited before running the script. 

## Scaling Relations
Currently available scaling relations:
 * Leonard2014_Interplate
 * WC1994 (i.e. Wells & Coppersmith (1994))
 * Strasser2010_Interface
 * Strasser2010_Intraslab
 * Blaser2010
 * Skarlatoudis2016

Other scaling relations used in ShakeMap/OQ have not been extended to have magnitude-length and magnitude-width relations. Note that rake information may be used to select the faulting style relation.

## GMPEs
The setup of GMPEs uses the same concept of weighted sets as used in ShakeMap, but in gmpe.conf the user explicitly creates the set and weights (i.e. tectonic weighting is not altered with location/depth as the template set is not designed to be geographically positioned).

## Earthquake Depth
Configure the 'hypo_depth' as the value to be used as centroid in the rupture plane. Once fault width exceeds this depth and/or the seismogenic depth the centroid depth may be altered within the script. In this case a note written to is written to log.

## Computing Ground Motion at Position/Distance
There are two wrapper scripts available, one designed to create FinDer templates, the other designed to create ground motion scenario data.

When using the FinDer template wrapper (makeFinDerTemplates.py), set the configuration parameter ["grid"]["compute"] to True and supply the griddkm value, i.e. the grid spacing in km. The script works internally with (lat, lon) values, as required by ShakeMap's structures, so the fault is placed approximately at the equator to make geographic to distance conversions simpler. A fixed vs30 is used and is configured in the configuration file also. The script will always create symmetric templates, and optionally *also* asymmetric templates if the "asym" configuration option is set to True. The rupture info files, as used by FinDer v3, can be created by setting the "rupinfo" configuration option to True.

When using the scenario wrapper (makeScenarioStnData.py), set the configuration parameter ["points"]["compute"] to True and set the path for a points file. That file should contain lines with space delimited fields: lat lon vs30 name. The first three columns will be used as input to the ground motion computation. The name will be used when generating a FinDer data_ file.

A third calculation option is for a single point (["pt"]["compute"] set to True) which simply computes the PGA for a single point at a specified Rjb (km) alongside the fault (i.e. Rx = Rjb and Ry0 = 0) and a specified vs30.
