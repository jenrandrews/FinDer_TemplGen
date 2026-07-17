import os
from numpy import arange, sqrt, zeros
#import fiona
#from fiona.model import to_dict
import matplotlib.pyplot as plt
import json
import shapely as shp
import shapely.ops as ops
import pyproj
from math import floor, ceil, degrees, radians, log10

# Local imports
import gocadascii as gocad

wgs84 = pyproj.Proj(init='epsg:4326')
nz = pyproj.Proj(init=f'epsg:2193')
project = pyproj.Transformer.from_proj(wgs84, nz)
rev_project = pyproj.Transformer.from_proj(nz, wgs84)

def CY_AspectRatio(mag, FNM = 0, FRV = 0):
    log_AR = (0.01752 - (0.00472 * FNM) - (0.01099 * FRV)) * pow(mag-4., 3.097)
    return pow(10., log_AR), pow(10., log_AR + 0.16), pow(10., log_AR - 0.16)

def scaling_relations(Mw):
    '''
    Opting for Murotani as pretty central in linear relation lines
    Murotani gives middle areas at low mag and smaller areas at high mag
    '''
    ### Murotani et al., 2013
    ### log10(M0) = log10(A) - log10(1.34e-10) * 3 / 2
    ### Mw = (log10(M0) + 7) - 16.05) / 1.5
    return pow(10., Mw - 3.83956187)

def scaling_relations_rev(area):
    return log10(area) + 3.83956187

################################################
# Work on NZ CFM Jen 2025
################################################

indir = 'Tsurf_individual'
name = 'hikurangi'
prefix = 'hik'
bPlot = False

if not os.path.isdir('event_xmls'):
    os.mkdir('event_xmls')

################################################
# Read in the 3D fault geometry from CFM
################################################

gcrd = gocad.AsciiReader()

all_blocks = []
for f in os.listdir(indir):
    if f.lower().find(name) == -1:
        continue
    if f.lower().find('plateau') > -1:
        continue
    blocks = gcrd.read(os.path.join(indir, f))
    all_blocks.append([b.atoms for b in blocks])

fault_surf = shp.geometry.MultiPoint([[p[0], p[1], p[2]] for block in all_blocks for b in block for p in b])
fault_cart = ops.transform(rev_project.transform, fault_surf)
bb = fault_cart.convex_hull
new_bb = bb.buffer(-0.01)
tight_bb = bb.buffer(-0.02)

x = [p.x for p in fault_cart.geoms if p.within(new_bb)]
y = [p.y for p in fault_cart.geoms if p.within(new_bb)]
z = [p.z/-1000. for p in fault_cart.geoms if p.within(new_bb)]

# Geographic
plt.gca()
plt.figure(figsize = (12,12))
trics = plt.tricontour(x, y, z, levels = arange(6,80,0.5), linewidths=0.5, colors='k')
#plt.plot([p[0] for p in bb.exterior.coords], [p[1] for p in bb.exterior.coords], c='g', lw=0.5)
plt.plot([p[0] for p in new_bb.exterior.coords], [p[1] for p in new_bb.exterior.coords], c='cyan', lw=0.5)
plt.plot([p[0] for p in tight_bb.exterior.coords], [p[1] for p in tight_bb.exterior.coords], c='magenta', lw=0.5)
cb = plt.scatter(x, y, c=z, marker=',', s=1)
plt.colorbar(cb)
#plt.xlim(xmin=177)
#plt.ylim(ymin=-38)
plt.savefig(f'{name}_geo.png')
#plt.show()

# Cartesian contours
plt.gca()
plt.figure(figsize = (12,12))
contour_dict = {}
contour_dict_cart = {}
for level, paths in zip(trics.levels, trics.get_paths()):
    if len(paths) == 0:
        continue
    cline = shp.geometry.LineString([[p[0], p[1], level] for p in paths.vertices if shp.geometry.Point(p[0], p[1]).within(tight_bb)])
    contour_dict[level] = cline
    contour_dict_cart[level] = ops.transform(project.transform, contour_dict[level])

    plt.plot([p[0] for p in cline.coords], [p[1] for p in cline.coords], lw=2, c='b')
    #plt.plot([p[0] for p in cline_simplify.coords], [p[1] for p in cline_simplify.coords], lw=1, c='r')
    #plt.scatter([c[0] for c in cline.coords], [c[1] for c in cline.coords], marker='x')
plt.savefig(f'{name}_contcompare.png')
#plt.show()

zlist = list(contour_dict_cart.keys())

# Need magnitude range: M7 - M9
jj = 0
templates = {}
print('Mag, area, len, wid, lat/lon')
for mw in arange(7., 9.05, 0.1):
    if bPlot:
        plt.figure(figsize = (8,8))
        for level in contour_dict:
            cline = contour_dict[level]
            plt.plot([p[0] for p in cline.coords], [p[1] for p in cline.coords], lw=0.5, c='b')

    templates[mw] = []
    
    # Get the scaling information for magnitude
    area = scaling_relations(mw)
    ar, ar_psd, ar_msd = CY_AspectRatio(mw, 0, 1)
    fwid = round(sqrt(area / ar))
    flen = area / fwid

    # Find the depth
    # Need to respect the coupling depth
    top = contour_dict_cart[zlist[0]]
    top_start = shp.geometry.Point([top.coords[0][0], top.coords[0][1], top.coords[0][2]])
    downdip_dists = [contour_dict_cart[lvl].distance(top_start)/1000. for lvl in contour_dict_cart]    
    # For largest mag, adjust AR for max width
    if fwid > max(downdip_dists):
        fwid = max(downdip_dists)
        flen = area / fwid
    downdip_step = max(10., fwid/10.)
    n_step = floor((downdip_dists[-1] - fwid)/downdip_step)
    downdip_step = (downdip_dists[-1] - fwid) / n_step
    downdip_contours = []
    for i in arange(n_step+1):
        dds = downdip_step * i
        dds_diffs = [abs(x-dds) for x in downdip_dists]
        downdip_contours.append(zlist[dds_diffs.index(min(dds_diffs))])

    # Down-dip steps
    for downdip_top in downdip_contours:
        top_contour = contour_dict_cart[downdip_top]
        tc_len = top_contour.length/1000.
        if tc_len < flen:
            print(f'Surface too small for event magnitude {mw}')
            continue
        lateral_step = max(10., flen/10.)
        n_step = floor((tc_len - flen)/lateral_step)
        if n_step > 0:
            lateral_step = (tc_len - flen) / n_step 
        else:
            lateral_step = 0.

        # Lateral steps
        for i in arange(n_step+1):

            jj += 1

            tc_start = top_contour.interpolate(i * lateral_step * 1000.)
            if tc_start is None:
                continue
            tc_start = shp.geometry.Point([tc_start.x, tc_start.y, tc_start.z])
            tc_end = top_contour.interpolate(((i * lateral_step) + flen) * 1000.)
            cdists = [abs(shp.geometry.Point(c).distance(tc_start)) for c in top_contour.coords]
            s = cdists.index(min(cdists))
            cdists = [abs(shp.geometry.Point(c).distance(tc_end)) for c in top_contour.coords]
            e = cdists.index(min(cdists))
            coords = [tc_start]
            coords.extend(top_contour.coords[s+1:e])
            coords.append(tc_end)
            tc_line = shp.geometry.LineString(coords)

            # Find the bottom contour
            zind = zlist.index(downdip_top)
            cdists = [abs((contour_dict_cart[lvl].distance(tc_start)/1000.)-fwid) for lvl in zlist[zind:]]
            bottom_zind = zind + cdists.index(min(cdists))
            bottom_contour = contour_dict_cart[zlist[bottom_zind]]
            bc_len = bottom_contour.length/1000.
            if bc_len < flen:
                print(f'Surface too small for event magnitude {mw}')
                continue
            bc_start = bottom_contour.interpolate(bottom_contour.project(tc_start))
            # Use fault length 
            bc_end = bottom_contour.interpolate((flen * 1000.) + bottom_contour.project(tc_start))
            # Use down-dip projection
            #bc_end = bottom_contour.interpolate(bottom_contour.project(tc_end))
            cdists = [shp.geometry.Point(c).distance(bc_start) for c in bottom_contour.coords]
            s = cdists.index(min(cdists))
            cdists = [shp.geometry.Point(c).distance(bc_end) for c in bottom_contour.coords]
            e = cdists.index(min(cdists))
            coords = [bc_start]
            coords.extend(bottom_contour.coords[s+1:e])
            coords.append(bc_end)
            bc_line = shp.geometry.LineString(coords)

            templ_coords = [c for c in tc_line.coords]
            templ_coords.extend([c for c in bc_line.coords[::-1]])
            templ = shp.geometry.Polygon(templ_coords)
            # try to remove narrow features / spikes
            small_templ = templ.buffer(-0.002)
            big_templ = small_templ.buffer(0.004)
            big_templ_cart = ops.transform(rev_project.transform, big_templ)
            tmp_templ_coords = [p for p in tc_line.coords if shp.geometry.Point(p).within(big_templ)]
            tmp_templ_coords.extend([p for p in bc_line.coords[::-1] if shp.geometry.Point(p).within(big_templ)])
            tmp_templ = shp.geometry.Polygon(tmp_templ_coords)
            tmp_templ_cart = ops.transform(rev_project.transform, tmp_templ)
            templ = tmp_templ
            ###

            flen_1 = tc_line.length
            flen_2 = bc_line.length
            fwid_1 = tc_line.distance(shp.Point(bc_line.coords[0]))
            fwid_2 = tc_line.distance(shp.Point(bc_line.coords[-1]))
            tc_cart = ops.transform(rev_project.transform, tc_line)
            bc_cart = ops.transform(rev_project.transform, bc_line)
            minlat = min([min([x[1] for x in tc_cart.coords]), min([x[1] for x in bc_cart.coords])])
            maxlat = max([max([x[1] for x in tc_cart.coords]), max([x[1] for x in bc_cart.coords])])
            minlon = min([min([x[0] for x in tc_cart.coords]), min([x[0] for x in bc_cart.coords])])
            maxlon = max([max([x[0] for x in tc_cart.coords]), max([x[0] for x in bc_cart.coords])])

            if not templ.is_valid or isinstance(templ, shp.geometry.multipolygon.MultiPolygon):
                if isinstance(templ, shp.geometry.polygon.Polygon):
                    plt.cla()
                    plt.plot([p[0] for p in templ.exterior.coords], 
                             [p[1] for p in templ.exterior.coords], lw=1, c='g')
                    templ = templ.buffer(0)
                if isinstance(templ, shp.geometry.multipolygon.MultiPolygon):
                    for geom in templ.geoms:
                        plt.plot([p[0] for p in geom.exterior.coords], 
                                 [p[1] for p in geom.exterior.coords], lw=1, c='b')
                    templ = shp.geometry.Polygon(templ.geoms[0])
                if bPlot:
                    plt.plot([p[0] for p in templ.exterior.coords], 
                             [p[1] for p in templ.exterior.coords], lw=1, c='r')
                    #plt.savefig(f'templ_{jj}.png')
                    plt.show()
                templ = shp.geometry.Polygon(templ.exterior.coords)
                templ_cart = ops.transform(rev_project.transform, templ)
            else:
                templ_cart = ops.transform(rev_project.transform, templ)

            if bPlot:
                plt.cla()
                plt.plot([p[0] for p in templ_cart.exterior.coords], 
                         [p[1] for p in templ_cart.exterior.coords], lw=1, c='r')
                plt.plot([p[0] for p in tmp_templ_cart.exterior.coords], 
                         [p[1] for p in tmp_templ_cart.exterior.coords], lw=1, c='b')
                plt.plot([p[0] for p in big_templ_cart.exterior.coords], 
                         [p[1] for p in big_templ_cart.exterior.coords], lw=1, c='g', ls=':')
                #plt.savefig('templ.png')
                plt.show()

#########################
            # Write the xml rupture file as a complex fault rupture (Open Quake)
#########################
            rake = 90.
            # Centroid
            cntrd = templ_cart.centroid
            ctrdists = {}
            for z in range(zind+1, bottom_zind):
                contour = contour_dict[zlist[z]]
                ctrdists[cntrd.distance(contour)] = zlist[z]
            sorted_ctrdists = sorted(ctrdists)
            pa = sorted_ctrdists[0] / (sorted_ctrdists[0] + sorted_ctrdists[1])
            cntrdz = ctrdists[sorted_ctrdists[0]] + (pa * (ctrdists[sorted_ctrdists[1]] - ctrdists[sorted_ctrdists[0]]))

            ofname = os.path.join('event_xmls', f'{prefix}_{mw:.1f}_{cntrd.x:.6f}_{cntrd.y:.6f}.xml')
            print(f'{ofname}, {mw:.1f};{scaling_relations_rev(templ.area/1e6):.1f}, {templ.area/1e6:.1f};{area:.1f}, {flen:.1f};{flen_1/1000.:.1f};{flen_2/1000.:.1f}, {fwid};{fwid_1/1000.:.1f};{fwid_2/1000.:.1f}, {minlat:.6f};{maxlat:.6f};{minlon:.6f};{maxlon:.6f}')

            ### Created the averaged gmms
            if os.path.isfile(os.path.join('hik_template_set/averaged_gmms', f'{prefix}_{mw:.1f}_{cntrd.x:.6f}_{cntrd.y:.6f}.csv.gz')):
                continue

            # Get all fault mesh points within the template bounds
            if True:
                plt.cla()
                templ_pts = [pt for pt in fault_cart.geoms if templ_cart.contains(pt)]
                templ_ext_pts = templ_cart.segmentize(max_segment_length = 0.01)
                zmin = min([min([p.z/-1000. for p in templ_pts]), min([p[2] for p in templ_ext_pts.exterior.coords])])
                zmax = max([max([p.z/-1000. for p in templ_pts]), max([p[2] for p in templ_ext_pts.exterior.coords])])
                plt.scatter([p.x for p in templ_pts], [p.y for p in templ_pts],
                            c=[p.z/-1000. for p in templ_pts], cmap='viridis',
                            marker=',', vmin=zmin, vmax=zmax)
                plt.scatter([p[0] for p in templ_ext_pts.exterior.coords],
                            [p[1] for p in templ_ext_pts.exterior.coords],
                            c=[p[2] for p in templ_ext_pts.exterior.coords],
                            marker=',', cmap='viridis', vmin=zmin, vmax=zmax)
                plt.scatter(templ_cart.centroid.x, templ_cart.centroid.y, marker='x', c='k')
                fname = f'rupture_{mw:.1f}_{cntrd.x:.6f}_{cntrd.y:.6f}.png'
                plt.savefig(fname)
                #plt.show()


            if os.path.isfile(ofname):
                #print(f'{ofname} already exists! DO NOT OVERWRITE')
                continue
            with open(ofname, 'w') as fout:
                contours = []
                fout.write("<?xml version='1.0' encoding='utf-8'?>\n")
                fout.write('<nrml xmlns:gml="http://www.opengis.net/gml" xmlns="http://openquake.org/xmlns/nrml/0.5">\n')
                fout.write('<complexFaultRupture>\n')
                fout.write(f'<magnitude>{mw:.1f}</magnitude>\n')
                fout.write(f'<rake>{rake}</rake>\n')
                fout.write(f'<hypocenter lon="{cntrd.x:.3f}" lat="{cntrd.y:.3f}" depth="{cntrdz:.3f}"/>\n')
                fout.write('<complexFaultGeometry>\n')
                tc_line_cart = ops.transform(rev_project.transform, tc_line)
                contours.append(tc_line_cart)
                fout.write('<faultTopEdge>\n')
                fout.write('<gml:LineString>\n')
                fout.write('<gml:posList>\n')
                for p in tc_line_cart.coords[::-1]:
                    fout.write(f'{p[0]:.6f} {p[1]:.6f} {downdip_top:.3f}\n')
                fout.write('</gml:posList>\n')
                fout.write('</gml:LineString>\n')
                fout.write('</faultTopEdge>\n')
                for z in range(zind+1, bottom_zind):
                    contour = contour_dict[zlist[z]]
                    templ_cont_cart = templ_cart.intersection(contour)
                    if isinstance(templ_cont_cart, shp.geometry.multilinestring.MultiLineString):
                        templ_cont_cart = shp.geometry.LineString(templ_cont_cart.geoms[0])
#                    if bPlot:
#                        plt.plot([p[0] for p in templ_cont_cart.coords], 
#                                 [p[1] for p in templ_cont_cart.coords], lw=0.5, c='r')
                    if templ_cont_cart.length < tc_line_cart.length * 0.7:
                        continue
                    fout.write('<intermediateEdge>\n')
                    fout.write('<gml:LineString>\n')
                    fout.write('<gml:posList>\n')
                    for p in templ_cont_cart.coords[::-1]:
                        fout.write(f'{p[0]:.6f} {p[1]:.6f} {zlist[z]:.3f}\n')
                    fout.write('</gml:posList>\n')
                    fout.write('</gml:LineString>\n')
                    fout.write('</intermediateEdge>\n')
                bc_line_cart = ops.transform(rev_project.transform, bc_line)
                fout.write('<faultBottomEdge>\n')
                fout.write('<gml:LineString>\n')
                fout.write('<gml:posList>\n')
                for p in bc_line_cart.coords[::-1]:
                    fout.write(f'{p[0]:.6f} {p[1]:.6f} {zlist[bottom_zind]:.3f}\n')
                fout.write('</gml:posList>\n')
                fout.write('</gml:LineString>\n')
                fout.write('</faultBottomEdge>\n')
                fout.write('</complexFaultGeometry>\n')
                fout.write('</complexFaultRupture>\n')
                fout.write('</nrml>\n')

        #plt.show()
#    if bPlot:
#        plt.savefig(f'{name}_{mw:.1f}_templates.png')

