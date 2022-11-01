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

## Depth
Configure the 'hypo_depth' as the value to be used as centroid in the rupture plane. Once fault width exceeds this depth and/or the seismogenic depth, the centroid depth may be altered within the script, with a note written to log.

## Position/Distance
The user only needs to supply the griddkm value, i.e. the grid spacing in km. The script works internally with (lat, lon) values, as required by ShakeMap's structures, so the fault is placed approximately at the equator to make geographic to distance conversions simpler.

## Asymmetric/Symmetric Templates
The script will always create symmetric templates, and optionally *also* asymmetric templates if the "asym" configuration option is set to True.
