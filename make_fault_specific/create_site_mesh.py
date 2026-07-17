# stdlib imports
import os
from math import floor, ceil
from numpy import sqrt
import geographiclib.geodesic as geo
import numpy as np
import pandas as pd
import shapely as shp

def getNZborder():
    '''
    Extract NZ border
    '''
    from cartopy.io import shapereader
    import shapefile

    resolution = '10m'
    category = 'cultural'
    name = 'admin_0_countries'

    shpfname = shapereader.natural_earth(resolution, category, name)
    recs = shapereader.Reader(shpfname).records()
    for rec in recs:
        if rec.attributes['ADMIN'] == 'New Zealand':
            poly = rec.geometry
    return poly

def vs30_mesh(prefix, lats, lons):
    '''
    For the input lats and lons, find the vs30 and write the site file
    Args:
    - prefix: prefix to file name
    - lats: latitudes
    - lons: longitudes
    '''
    poly = getNZborder()
    # Grid vs30 data
    gdf = pd.read_csv(f'{os.path.dirname(os.path.abspath(__file__))}/sites/National_grid_1km.csv')
    # Create the output file
    with open(f'{prefix}_sites.csv', 'w') as fout:
        fout.write('lon,lat,vs30,z1pt0,z2pt5,vs30measured,backarc\n')
        for lat in lats:
            for lon in lons:
                if poly.contains(shp.Point(lon, lat)):
                    gdf['distance'] = sqrt((gdf['latitude'] - lat) ** 2 + (gdf['longitude'] - lon) ** 2)
                    mindist = min(gdf['distance'].tolist())
                    site = gdf.loc[gdf['distance'] == mindist]
                    vs30 = site.iloc[0]['Vs30']
                    z1pt0 = site.iloc[0]['Z1pt0']
                    z2pt5 = site.iloc[0]['Z2pt5']
                else:
                    vs30 = 760.
                    z1pt0 = -1.
                    z2pt5 = -1.
                ostr = f"{lon},{lat},{vs30},{z1pt0},{z2pt5},False,False\n"
                fout.write(ostr)
    return

def mesh_from_bb(dkm, minlat, maxlat, minlon, maxlon):
    '''
    Make a mesh for use as distance context based on extent
    Args:
    - dkm
    - minlat: southern latitude
    - maxlat: northern latitude
    - minlon: western longitude
    - maxlon: eastern longitude
    Return:
    - mesh
    '''
    # Longitude
    xlim = geo.Geodesic.WGS84.Inverse(minlat, minlon, minlat, maxlon)['s12']/1000.
    nx = ceil(xlim/dkm)
    xdist = nx * dkm
    xmax = geo.Geodesic.WGS84.Direct(minlat, minlon, 90., xdist*1000.)['lon2']
    if xmax < 0.:
        xmax += 360.
    # Latitude
    ylim = geo.Geodesic.WGS84.Inverse(minlat, minlon, maxlat, minlon)['s12']/1000.
    ny = ceil(ylim/dkm)
    ydist = ny * dkm
    ymax = geo.Geodesic.WGS84.Direct(minlat, minlon, 0., ydist*1000.)['lat2']
    # Grid
    dlat = (ymax - minlat)/ny
    lats = [minlat + n*dlat for n in range(ny)]
    dlon = (xmax - minlon)/nx
    #lons = [minlon + n*dlon for n in range(nx) if minlon + n*dlon < 180.]
    lons = [minlon + n*dlon for n in range(nx)]
    # Sites 
    #mlons, mlats = np.meshgrid(lons, lats)
    return lats, lons


