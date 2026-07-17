import os
import fnmatch
import gzip
from numpy import savetxt, zeros
from math import log10
import shapely as shp
import shapely.ops as ops
import pyproj

def rdRupFile(fname):
    '''
    Read the rupture file
    Args:
        fname: file name
    Return:
        rup_lines: dictionary of rupture lines, key is depth, value is contour lat, lon
    '''
    rup_lines = {}
    if not os.path.isfile(fname):
        return rup_lines
    with open(fname, 'r') as fin:
        bLine = False
        for l in fin:
            if l.startswith('<'):
                continue
            else:
                fs = l.split()
                z = fs[2]
                if z not in rup_lines:
                    rup_lines[z] = [[float(fs[0]), float(fs[1])]]
                else:
                    rup_lines[z].append([float(fs[0]), float(fs[1])])
    return rup_lines

def rdRupSummaryFile(fname):
    '''
    Read the log file from creation which has the 
    Args:
        fname: file name
    Return:
        rup_params: dictionary
    '''
    rup_params = {}
    with open(fname, 'r') as fin:
        headers = [v.rstrip().lstrip() for v in fin.readline().split(',')]
        for l in fin:
            fs = l.split(',')
            rup_params[fs[0].split(os.sep)[-1]] = {}
            for key, val in zip(headers, fs[1:]):
                rup_params[fs[0].split(os.sep)[-1]][key] = [float(v) for v in val.split(';')]
    return rup_params

def rdTemplSetFile(fname):
    '''
    Read the tset file
    Args:
        fname: file name
    Return:
        tsets: return dictionary of template set
    '''
    tsets = {}
    with open(fname, 'r') as fin:
        for l in fin:
            fs = l.split()
            tsets[fs[0]] = eval(' '.join(fs[2:]))
    return tsets

def rdGroundMotion(fname):
    '''
    Read the ground motion file
    Args:
        fname: file name
    Return:
        pga: list of pga values, units of g as output by OpenQuake
    '''
    pga = []
    if os.path.isfile(fname):
        with open(fname, 'r') as fin:
            fin.readline()
            fin.readline()
            for l in fin:
                fs = l.split(',')
                pga.append(float(fs[1]))
    elif os.path.isfile(f'{fname}.gz'):
        with gzip.open(f'{fname}.gz', 'rt') as fin:
            fin.readline()
            fin.readline()
            for l in fin:
                fs = l.split(',')
                pga.append(float(fs[1]))
    return pga

def rdSites(fname, b_setOnly = True):
    '''
    Read the sites file
    Args:
        fname: file name
    Return:
        lats: unique lats
        lons: unique lons
    '''
    lats = []
    lons = []
    with open(fname, 'r') as fin:
        fin.readline()
        for l in fin:
            fs = l.split(',')
            lon = float(fs[0])
            lat = float(fs[1])
            if not b_setOnly or (b_setOnly and lat not in lats):
                lats.append(lat)
            if not b_setOnly or (b_setOnly and lon not in lons):
                lons.append(lon)
    return lats,lons

def rdSiteMesh(fname):
    '''
    Read the site mesh file
    Args:
        fname: file name
    Return:
        sites: list of lat, lon 
    '''
    lats = []
    lons = []
    if not os.path.isfile(fname):
        return lats, lons
    with open(fname, 'r') as fin:
        fin.readline()
        for l in fin:
            fs = l.split(',')
            lats.append(float(fs[2]))
            lons.append(float(fs[1]))
    return lats, lons

def formatHeader(oname):
    '''
    Format the FinDer template header by replacing default output with FinDer specific format.
    Operation done on file.
    Args:
        oname: template file name
    '''
    fin = open(oname, 'r')
    text = fin.read()
    text = text.replace('# ','')
    text = text.replace('\n\n','\n')
    fin.close()
    fout = open(oname, 'w')
    fout.write(text)
    fout.close()
    return

def createCentroidPolygon(polygon_dist, templ_bbs):
    '''
    Create a polygon around an input set of templ_bbs
    Args:
    - polygon_dist: polygon distance from templ_bbs in km
    - templ_bbs: list of quadrilaterals making the bounding boxes of the templates
    Return:
    -  polygon
    '''
    if len(templ_bbs) < 2:
        # This should be a circle if 1 point
        return None

    flat_lats = list(sum([(x[0], x[1]) for x in templ_bbs], ()))
    flat_lons = list(sum([(x[2], x[3]) for x in templ_bbs], ()))
    minlat = min(flat_lats)
    maxlat = max(flat_lats)
    minlon = min(flat_lons)
    maxlon = max(flat_lons)

    # Set up the projections
    wgs84 = pyproj.Proj(init='epsg:4326')
    nz = pyproj.Proj(init=f'epsg:2193')
    project = pyproj.Transformer.from_proj(wgs84, nz)
    rev_project = pyproj.Transformer.from_proj(nz, wgs84)

    # Create line and buffer
    bb = shp.geometry.Polygon(((minlon, minlat), (minlon, maxlat), (maxlon, maxlat), (maxlon, minlat), (minlon, minlat)))
    bb_cart = ops.transform(project.transform, bb)
    poly_cart = bb_cart.buffer(polygon_dist * 1000, join_style=2)
    poly = ops.transform(rev_project.transform, poly_cart)
    npoly = poly.simplify(0.01, preserve_topology=False)

    if False:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter([p[0] for p in bb.exterior.coords], [p[1] for p in bb.exterior.coords], marker='x', c='r')
        plt.plot([p[0] for p in bb.exterior.coords], [p[1] for p in bb.exterior.coords], c='r')
        plt.plot([p[0] for p in npoly.exterior.coords], [p[1] for p in npoly.exterior.coords], c='b')
        plt.savefig('buffer.png')
    return npoly

def writeSetFile(name, templ_res, minflen, maxflen, minmag, maxmag, poly):
    '''
    Format the FinDer template set summary file
    Args:
       name: name for the template set
       templ_res: template resolution in km
       minflen: minimum fault length in the template set
       maxflen: maximum fault length in the template set
       minmag: minimum magnitude in the template set
       maxmag: maximum magnitude in the template set
       poly: centroid polygon
    '''
    ofile = f"{name}.txt"
    if os.path.isfile(ofile):
        return

    npstr = len(poly.exterior.coords)
    pstr = ' '.join([f'{p[1]:.6f},{p[0]:.6f}' for p in poly.exterior.coords])
    fout = open(ofile, 'w')
    fout.write(f"TEMPLATE_SET_NAME   {name}\n")
    fout.write("TEMPLATE_TYPE       SPECIFIC\n")
    fout.write(f"TEMPLATE_DIRECTORY  conf/Templates_PGA_{name}\n")
    fout.write(f"D_KM                {templ_res:.1f}\n")
    fout.write("LATITUDE_FILE       latitude.dat\n")
    fout.write("LONGITUDE_FILE      longitude.dat\n")
    fout.write(f"MIN_LENGTH          {minflen-1:.0f}\n")
    fout.write(f"MAX_LENGTH          {maxflen+1:.0f}\n")
    fout.write(f"MIN_MAG             {minmag-0.5:.1f}\n")
    fout.write(f"MAX_MAG             {maxmag:.1f}\n")
    fout.write("FAST_MAG_RANGE      1.0\n")
    fout.write("FAST_LATLON_RANGE   1.0\n")
    fout.write("TEMPLATE_INFO_FILE  template_info.txt\n")
    fout.write(f"CENTROID_POLYGON    {npstr} {pstr}\n")
    fout.close()
    return

def writeLatLon(odir, lats, lons):
    '''
    Format the latitude and longitude files for the FinDer template set
    Args:
        lats: list of latitudes
        lons: list of longitudes
    '''
    ofile = os.path.join(odir, 'latitude.dat')
    if not os.path.isfile(ofile):
        with open(ofile, 'w') as fout:
            for lat in sorted(lats):
                fout.write(f'{lat:.6f}\n')
    ofile = os.path.join(odir, 'longitude.dat')
    if not os.path.isfile(ofile):
        with open(ofile, 'w') as fout:
            for lon in sorted(lons):
                fout.write(f'{lon:.6f}\n')
    return

def writeFiles(odir, rfile, rup_params, templ_res):
    '''
    Write the FinDer template files
    Args:
        odir: name of the template set
        rfile: name of the rupture file
        rup_params[rup]: rupture parameter summary
    '''

    # Write rupture file
    rup = rdRupFile(os.path.join('event_xmls', rfile))
    minz = min(rup.keys())
    maxz = max(rup.keys())
    coords = [x for x in rup[minz]]
    coords.extend(rup[maxz][::-1])
    coords.append(rup[minz][0])
    polygon = shp.Polygon(coords)
    centroid = polygon.centroid
    if False:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter([p[0] for p in polygon.exterior.coords], [p[1] for p in polygon.exterior.coords], marker='x', c='r')
        plt.plot([p[0] for p in polygon.exterior.coords], [p[1] for p in polygon.exterior.coords], c='r')
        plt.scatter(centroid.x, centroid.y, marker='x', c='k')
        plt.show()
    ofile_rup = f'rupture_L{max(rup_params["len"]):.6f}_W{max(rup_params["wid"]):.4f}_{centroid.y:.4f}_{centroid.x:.4f}.txt'
    ofile = os.path.join(odir, ofile_rup)
    if not os.path.isfile(ofile):
        with open(ofile, 'w') as foutr:
            for pt in rup[minz]:
                foutr.write(f'{pt[1]} {pt[0]} {minz}\n')
            for pt in rup[maxz][::-1]:
                foutr.write(f'{pt[1]} {pt[0]} {maxz}\n')
            foutr.write(f'{rup[minz][0][1]} {rup[minz][0][0]} {minz}\n')
    else:
        print(f'{ofile} already exists')
    if False:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter([p[0] for p in rup[minz]], [p[1] for p in rup[minz]], marker='x', c='r')
        plt.scatter([p[0] for p in rup[maxz][::-1]], [p[1] for p in rup[maxz][::-1]], marker='x', c='b')
        plt.scatter(rup[minz][0][0], rup[minz][0][1], marker='o', c='r')
        plt.show()

    # Write template file
    gms = rdGroundMotion(os.path.join('averaged_gmms', rfile.replace('.xml', '.csv')))

    ## Use sitemesh if available:
    sitemesh = os.path.join('site_meshes', f'{odir}_sites.csv')
    if os.path.isfile(sitemesh):
        lats, lons = rdSiteMesh(sitemesh)
    else:
        lats, lons = rdSites(os.path.join('site_meshes', f'{ts}_sites.csv'), b_setOnly = False)
    if len(lats) != len(gms):
        print(f'ERROR for {ofile}: {len(gms),} {len(lats)}')
        return

    templ_gm = {}
    for lat, lon, gm in zip(lats, lons, gms):
        key = f'{lon:.6f}/{lat:.6f}'
        templ_gm[key] = gm
    lats = sorted(set(lats))
    lons = sorted(set(lons))
    gm_array = zeros((len(lats), len(lons)))
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            key = f'{lon:.6f}/{lat:.6f}'
            gm_array[i][j] = log10(templ_gm[key] * 980.665) # convert from g to cm/s/s
    ofile_templ = f'template_L{max(rup_params["len"]):.6f}_W{max(rup_params["wid"]):.4f}_{centroid.y:.4f}_{centroid.x:.4f}.txt'
    ofile = os.path.join(odir, ofile_templ)
    if not os.path.isfile(ofile):
        hstr = f'{len(lons):d} {len(lats):d}\n{max(rup_params["len"]):f} 0 {templ_res:.1f}\n'
        savetxt(ofile, gm_array, fmt='%.6e', header=hstr)
        formatHeader(ofile)
    else:
        print(f'{ofile} already exists')

    # Write template_info.txt file
    ofile = os.path.join(odir, 'template_info.txt')
    with open(ofile, 'a') as fout:
        fout.write(f'{max(rup_params["len"]):.6f} {max(rup_params["wid"]):.6f} {centroid.y:.6f} {centroid.x:.6f} {rup_params["Mag"][0]:.1f} {ofile_templ} {ofile_rup}\n')

    return

if __name__ == '__main__':
    
    templ_res = 5.
    centroid_dist = 20.

    # Read tset file
    tsets = rdTemplSetFile('tsets.dat')

    # Read rup summary file
    rup_params = rdRupSummaryFile('templ_mag_area.log')

    # For each tset
    for ts in tsets:
        tset = f'Templates_PGA_{ts}'
        if not os.path.isdir(tset):
            os.mkdir(tset)
        # Write lat / lon files
        lats, lons = rdSites(os.path.join('site_meshes', f'{ts}_sites.csv'))
        writeLatLon(tset, lats, lons)
        # Write summary file
        flens = [rup_params[x]['len'] for t in tsets for x in tsets[t]]
        minflen = min([min(x) for x in flens])
        maxflen = max([max(x) for x in flens])
        mags = [float(x.split('_')[1]) for t in tsets for x in tsets[t]]
        minmag = min(mags)
        maxmag = max(mags)
        coords = [rup_params[x]['lat/lon'] for t in tsets for x in tsets[t]]
        poly = createCentroidPolygon(centroid_dist, coords)
        writeSetFile(ts, templ_res, minflen, maxflen, minmag, maxmag, poly)

        # For each templ
        for rup in tsets[ts]:
            writeFiles(tset, rup, rup_params[rup], templ_res)
