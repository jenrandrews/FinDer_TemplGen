# Creating Fault-Specific Templates for New Zealand

This is an initial check-in of ungeneralised scripts to create bespoke fault-specific templates for
New Zealand. 

## Create rupture files: in = surface; out = rupture xml files
### extract_fault_3d.py
1. Read NZ Community Fault Model surface file
2. Read 3D geometry as a shapely MultiPoint, convert to cartesian (NZTM) from lat/lon
3. Clip the area to remove some weird edge effects by creating a convex hull and trimming by a small factor
4. Use tricontour to contour the surface into depth contours at 0.5 km intervals between 6km and 80 km depth. Save each contour as a shapely LineString
5. For magnitudes in magnitude range:
  a. Use scaling relation to get area and aspect ratio to get length and width. NOTE only one aspect ratio is used.
  b. Find the top contour for each set of downdip ruptures:
    i. Using the shallowest contour as the starting point, compute the down dip distance to all contours. If magnitude requires a wider fault than available, adjust trial length and set width to max (i.e. aspect ratio not preserved).
    ii. Set the downdip step using the larger of 10 km or 10% rupture width.
    iii. Find the set of top contours that closest match the desired set of downdip steps
  c. For each downdip step:
    i. Find the number of lateral steps using the step of larger of 10 km or 10% rupture length and computing the number that fit along the top contour length.
    ii. For lateral step:
      1. Find the start and end points (x, y, z) on the top contour using desired rupture length along the contour. Use all the coords on the contour between desired start and end to form a LineString.
      2. Find the bottom contour that most closely matches the desired rupture width from the top contour and set the start point as the closest point on the bottom contour from the start point on the top contour. Set the end point at desired rupture length from the start point on the bottom contour. Use all the coords on the bottom contour between desired start and end to form a LineString.
      3. Use the top and bottom coords to form a Polygon, do some tweaks to remove narrow features and spikes.
      4. Convert rupture patch edge coords back to lat/lon and use shapely checks for validity (i.e. single polygon). 
      5. Compute a centroid (x, y, z).
      6. Create an openquake complexFaultRupture xml file and specify: magnitude, rake, hypocentre (centroid), and all contours on the rupture from our contour set (at 0.5 km downdip spacing). These are specified as faultTopEdge, a set of intermediateEdge and a faultBottomEdge.

## Create computation point files (i.e. template mesh files): in = rupture xml files; out = template sets files, mesh files
### create_tset_site_meshes.py and create_site_meshes.py
1. Find the min/max lat/lon extents of each of the rupture files, applying a 200 km buffer around the rupture (I suspect I use 200 km because of the scale of NZ).
2. Template set sizes range from 4 * 4 degrees to max size from (1) + 2 degrees, with 2 degree steps. Walk through the sizes starting smallest.
  a. Walk through lat range and the lon range, step size of 2 degrees.
  b. If rupture extent from (1) sits within this set extent, assign it.
  c. Walk through the template sets and, if there are ruptures assigned, crop the region for NZ bounds (no point having templates beyond where we have observation data).
  d. Create the mesh using the min/max lat/lon and desired resolution, e.g. 5 km.

## Create OpenQuake .ini files: in = template set and mesh; out = sites and .ini files for openquake
### generate_inis.py
1. For each template set, create a sites file from the mesh, which is a csv file with the lat/lon points and vs30
2. Create the .ini file for a scenario calculation in OpenQuake. The rupture mesh spacing is specified here: 2 km, and xml rupture file given. The logic tree is specified and can be set as using an average and mean ground motion. Set intensity measure type.

## Run OpenQuake: in = input openquake files; out = openquake output computation files
### run_oq_segmented.sh
1. Run the oq engine with the ini files created.

## Create FinDer input files: in = openquake files; out = FinDer files
### out2templ.py
1. Convert the OQ output to FinDer configuration files.
2. Create additional files, e.g. create the centroid polygon

