#!/usr/bin/env python

#############
# N.B. This script is based heavily on shakemap_src/tests/shakelib/multigmpe_test.py
#############

# stdlib imports
import copy
import os
import sys
import logging
from math import log10, sqrt, exp, floor, sin, cos, tan, radians, ceil, atan, degrees
import numpy as np

# openquake imports
from openquake.hazardlib import imt, const
from openquake.hazardlib.gsim.base import RuptureContext
from openquake.hazardlib.gsim.base import DistancesContext
from openquake.hazardlib.gsim.base import SitesContext
from openquake.hazardlib.geo import geodetic, Point, Line, Mesh, RectangularMesh
from openquake.hazardlib.geo import utils
from openquake.hazardlib.geo.surface import PlanarSurface, SimpleFaultSurface, ComplexFaultSurface

# shakemap imports
from shakelib.conversions.imc.boore_kishida_2017 import BooreKishida2017
from shakelib.multigmpe import MultiGMPE, set_sites_depth_parameters
import shakelib.sites as sites
from shakelib.sites import Sites

# local imports
import Leonard2014_Interplate_Ext
import WC1994_Ext
import Blaser_2010
import Skarlatoudis_2016
import Strasser2010_Interface_Ext
import Strasser2010_Intraslab_Ext

logger = logging.getLogger(__name__)

def lng2cm(inx):
    '''
    Convert PGA in ln(g) to log10(cm/s/s)
    Args:
    - Numpy array of PGA values in ln(g)
    Return:
    - Numpy array of PGA values in log10(cm/s/s)
    '''
    inx = np.log10(np.exp(inx)) + log10(980.665)
    return inx


def getScalingRelation(conf):
    '''
    Returns the ScalingRelation class using configuration file input.
    Allowing for extension to other scaling relations. However, all will require an Ext class to 
    provide length and width rather than area.
    Args:
    - conf: configuration dictionary
    Return:
    - ScalingRelation class object
    '''
    if conf['scaling_relation']['scalrel'] == 'Leonard2014_Interplate':
        return Leonard2014_Interplate_Ext.Leonard2014_Interplate_Ext()
    elif conf ['scaling_relation']['scalrel'] == 'WC1994':
        return WC1994_Ext.WC1994_Ext()
    elif conf['scaling_relation']['scalrel'] == 'Strasser2010_Interface':
        return Strasser2010_Interface_Ext.Strasser2010_Interface_Ext()
    elif conf['scaling_relation']['scalrel'] == 'Strasser2010_Intraslab':
        return Strasser2010_Intraslab_Ext.Strasser2010_Intraslab_Ext()
    elif conf['scaling_relation']['scalrel'] == 'Blaser2010':
        return Blaser_2010.Blaser2010_Interface()
    elif conf['scaling_relation']['scalrel'] == 'Skarlatoudis2016':
        return Skarlatoudis_2016.Skarlatoudis_2016()
    else:
        return Leonard2014_Interplate_Ext.Leonard2014_Interplate_Ext()


def importConfig(fname):
    '''
    Imports a configuration file into a python dictionary. 
    Configuration file format expected (json) so that 'eval' works.
    Lines beginning with # are ignored.
    Args:
    - fname: filename
    Return:
    - Dictionary with contents of configuration file
    '''
    fin = open(fname, 'r')
    text = fin.read()
    ntext = ''
    fin.close()
    for l in text.split('\n'):
        l = l.rstrip().lstrip()
        if l.startswith('#'):
            continue
        ntext += l
    ddict = eval(ntext)
    return ddict


def importRuptureContext(evconf):
    '''
    Creates RuptureContext based on configuration file input for a fault mesh
    Args:
    - evconf: event configuration
    Return:
    - rctx: RuptureContext for rupture
    - faultplane: PlanarSurface for rupture
    '''
    from json import JSONDecoder
    with open(evconf['evmech']['geometry'], 'r') as fin:
        rup = JSONDecoder().decode(fin.read())
    evconf['mag'] = rup['metadata']['mag']
    if 'evloc' not in evconf:
        evconf['evloc'] = {}
    evconf['evloc']['centroid_lat'] = rup['metadata']['lat']
    evconf['evloc']['centroid_lon'] = rup['metadata']['lon']
    evconf['evloc']['hypo_depth'] = rup['metadata']['depth']
    ptlist = []
    allpts = len(rup['features'][0]['geometry']['coordinates'][0][0])
    for pt in rup['features'][0]['geometry']['coordinates'][0][0]:
        ptlist.append(Point(depth=pt[2], latitude=pt[1], longitude=pt[0]))
#    faultplane = SimpleFaultSurface.from_fault_data(fault_trace=Line(ptlist[:allpts//2]), 
#            upper_seismogenic_depth=ptlist[0].depth,
#            lower_seismogenic_depth=ptlist[-2].depth, 
#            dip=evconf['evmech']['dip'], 
#            mesh_spacing=1.0)
    topedge = ptlist[:allpts//2]
    bottomedge = ptlist[allpts//2:-1][::-1]
    faultplane = ComplexFaultSurface.from_fault_data(
            edges=[Line(topedge), Line(bottomedge)],
            mesh_spacing=1.0)
    rctx = RuptureContext()
    rctx.strike = faultplane.get_strike()
    evconf['evmech']['strike'] = rctx.strike
    rctx.rake = evconf['evmech']['rake']
    rctx.dip = evconf['evmech']['dip']
    rctx.mag = evconf['mag']
    rctx.hypo_depth = evconf['evloc']['hypo_depth']
    rctx.width = faultplane.get_width()
#    rctx.width = geodetic.distance(ptlist[0].longitude, ptlist[0].latitude, ptlist[0].depth, 
#            ptlist[-2].longitude, ptlist[-2].latitude, ptlist[-2].depth)
    rctx.ztor = 0.1
    if False:
        import matplotlib.pyplot as plt
        plt.plot([p.longitude for p in ptlist], [p.latitude for p in ptlist])
        plt.scatter(faultplane.mesh.lons, faultplane.mesh.lats)
        plt.savefig(evconf['evmech']['geometry'].replace('.json', '.png'))
    return evconf, [rctx], [faultplane]


def getKaklamanosCentroidDepth(rctx):
    '''
    Centroid using Kaklamanos et al. (2010)
    Args:
    - rctx: rupture context
    Return:
    - centroid depth
    '''
    if rctx.rake is None:
        return 7.08 + (0.61 * rctx.mag) # Kaklamanos eq.4 general
    if (-45 <= rctx.rake <= 45) or (rctx.rake >= 135) or (rctx.rake <= -135):
        return 5.63 + (0.68 * rctx.mag) # Kaklamanos eq.4 SS
    else:
        return 11.24 - (0.2 * rctx.mag) # Kaklamanos eq.4 non-SS


def adjustRCTXparams(rctx, seismogenic_depth, centroid_depth):
    '''
    Adjust the RCTX parameters so that depths are sensible
    Use initial hypo depth as centroid depth unless width requires us to move it
    Hypo depth empirically found to be at about 60% depth, but simpler to treat as 50% here
    Args:
    - rctx: the current rupture context
    - seismogenic_depth: seismogenic_depth from configuration
    - centroid_depth: centroid depth
    Return:
    - rctx: rupture context with modifications
    - maxz: maximum depth
    - offset: horizontal offset due to dip
    - centroid_depth: updated centroid depth
    '''
    diprad = radians(rctx.dip)
    if rctx.width > seismogenic_depth * sin(diprad):
        rctx.width = seismogenic_depth * sin(diprad)
        rctx.ztor = 0.1
        maxz = seismogenic_depth
        centroid_depth = (maxz-rctx.ztor)/2.
        logger.info('N.B. new centroid depth (too wide): %.2f' % (maxz/2.)) 
    elif centroid_depth < sin(diprad) * rctx.width * 0.5:
        rctx.ztor = 0.1
        maxz = sin(diprad) * rctx.width
        if maxz > seismogenic_depth:
            rctx.width = seismogenic_depth * sin(diprad)
            maxz = seismogenic_depth
        centroid_depth = (maxz-rctx.ztor)/2.
        logger.info('N.B. new centroid depth (too shallow): %.2f' % (maxz/2.)) 
    elif centroid_depth + (sin(diprad) * rctx.width * 0.5) > seismogenic_depth:
        maxz = seismogenic_depth
        rctx.ztor = maxz - (sin(diprad) * rctx.width)
        if rctx.ztor < 0.1:
            rctx.width = seismogenic_depth * sin(diprad)
            rctx.ztor = 0.1
        centroid_depth = (maxz-rctx.ztor)/2.
        logger.info('N.B. new centroid depth (too deep): %.2f' % (maxz/2.)) 
    else:
        offset = (sin(diprad) * rctx.width * 0.5)
        rctx.ztor = centroid_depth - offset
        maxz = centroid_depth + offset
    return rctx, maxz, centroid_depth


def createCentroidPolygon(calcconf, centroids):
    '''
    Create a polygon around an input set of centroids
    Args:
    - calcconf: calculation configuration
    - centroids: list of [lon, lat] points
    Return:
    -  polygon
    '''
    if len(centroids) < 2:
        # This should be a circle if 1 point
        return None
    if 'centroid_polygon_dist' not in calcconf['fault-specific']:
        return None

    import shapely as shp
    import shapely.ops as ops
    import pyproj
    # Set up the projections
    wgs84 = pyproj.Proj(init='epsg:4326')
    nz = pyproj.Proj(init=f'epsg:{calcconf["fault-specific"]["epsg"]}')
    project = pyproj.Transformer.from_proj(wgs84, nz)
    rev_project = pyproj.Transformer.from_proj(nz, wgs84)
    # Create line and buffer
    line = shp.geometry.LineString(centroids)
    line_cart = ops.transform(project.transform, line)
    xcorr = calcconf['fault-specific']['centroid_polygon_dist']
    poly_cart = line_cart.buffer(xcorr * 1000, cap_style=1)
    poly = ops.transform(rev_project.transform, poly_cart)
    npoly = poly.simplify(0.01, preserve_topology=False)
    return npoly


def createSubFaultRuptureContexts(evconf, calcconf):
    '''
    Create fault-specific ruptures for magnitudes in range
    Args:
    - evconf: event config
    - calcconf: calculation config
    Return:
    - lists, RuptureContexts and faultplanes
    '''
    import shapely as shp
    import shapely.ops as ops
    import pyproj
    from json import JSONDecoder

    # Set up the projections
    wgs84 = pyproj.Proj(init='epsg:4326')
    nz = pyproj.Proj(init=f'epsg:{calcconf["fault-specific"]["epsg"]}')
    project = pyproj.Transformer.from_proj(wgs84, nz)
    rev_project = pyproj.Transformer.from_proj(nz, wgs84)

    # Read in the fault geojson
    with open(evconf['evmech']['geometry'], 'r') as fin:
        infault = JSONDecoder().decode(fin.read())
    fault = shp.geometry.LineString([(p[0], p[1]) for p in infault['coordinates']])
    fault_cart = ops.transform(project.transform, fault)

    # Set parameters that won't change
    seismogenic_depth = evconf['seisstruc']['seismogenic_depth']
    scalrel = getScalingRelation(calcconf)
    flen = scalrel.get_median_length(evconf['mag'], evconf['evmech']['rake'])
    fwid = scalrel.get_median_width(evconf['mag'], evconf['evmech']['rake'])
    if flen > fault_cart.length / 1000.:
        logger.error(f'Fault is too small {fault_cart.length:.4f} for a M{evconf["mag"]:.1f} event {flen:.4f}')
        return None, None
    # Note that mag range extends width relation beyond validity, so fix assumes
    # len = wid when wid > len, i.e. aspect ratio 1.
    if fwid > flen: 
        fwid = flen
    farea = scalrel.get_median_area(evconf['mag'], evconf['evmech']['rake'])
    logger.info(f'Fault dimensions: {farea}, {flen}, {fwid}, {flen*fwid}, {flen/fwid}')
    diprad = radians(evconf['evmech']['dip'])
    xcorr = cos(diprad) * fwid

    # Set the iteration for multiple overlapping fault patches
    # Overlap is 10% or 20 km
    dist_step = 0
    overlap = min([20.*1000., flen*1000.*0.1])
    n_subfaults = round(((fault_cart.length)-(flen*1000.))/overlap)
    if n_subfaults > 0:
        dist_step = ((fault_cart.length)-(flen*1000.)) / n_subfaults
    l_rctx = []
    l_faultplane = []
    for sind in range(n_subfaults + 1):
        # Create rupture context
        rctx = RuptureContext()
        rctx.rake = evconf['evmech']['rake']
        rctx.dip = evconf['evmech']['dip']
        rctx.mag = evconf['mag']
        rctx.hypo_depth = evconf['evloc']['hypo_depth']
        rctx.width = fwid
        centroid_depth = evconf['evloc']['hypo_depth']
        if centroid_depth is None:
            centroid_depth = getKaklamanosCentroidDepth(rctx)
        rctx, maxz, centroid_depth = adjustRCTXparams(rctx, seismogenic_depth, centroid_depth)

        # Step through and create sub fault
        start_dist = sind * dist_step 
        end_dist = start_dist + (flen * 1000.)
        sbf = ops.substring(fault_cart, start_dist, end_dist)
        dx = sbf.coords[-1][0] - sbf.coords[0][0]
        dy = sbf.coords[-1][1] - sbf.coords[0][1]
        theta = degrees(atan(dx/dy))
        if sbf.coords[0][1] > sbf.coords[-1][1]:
            approx_strike = 180. + theta
        elif sbf.coords[0][0] > sbf.coords[-1][0]:
            approx_strike = 360. + theta
        else:
            approx_strike = theta
        # Start each with positive
        xcorr = abs(xcorr)
        side = 'left'
        delta_strike = abs(approx_strike - evconf['evmech']['strike'])
        delta_strike = delta_strike if delta_strike < 180. else 360. - delta_strike 
        if delta_strike < 90. and xcorr > 0.:
            xcorr *= -1.
            side = 'right'
        if int(shp.__version__.split('.')[0]) >= 2:
            b_sbf = sbf.offset_curve(xcorr * 1000., join_style=2)
        else:
            b_sbf = sbf.parallel_offset(abs(xcorr * 1000), side, join_style=2)
        sb_geo = ops.transform(rev_project.transform, sbf)
        b_sb_geo = ops.transform(rev_project.transform, b_sbf)
        topedge = [Point(depth=rctx.ztor, latitude=pt[1], longitude=pt[0]) for pt in sb_geo.coords]
        bottomedge = [Point(depth=maxz, latitude=pt[1], longitude=pt[0]) for pt in b_sb_geo.coords]
        faultplane = ComplexFaultSurface.from_fault_data(
            edges=[Line(topedge), Line(bottomedge)],
            mesh_spacing=1.0)
        l_faultplane.append(faultplane)
        l_rctx.append(rctx)
    return l_rctx, l_faultplane


def createRuptureContext(evconf, calcconf):
    '''
    Creates RuptureContext based on configuration file input for a simple, single planar surface
    defined by centroid latitude and longitude, strike, dip and rake. Magnitude and scaling 
    relation (from configuration file) are used to determine fault length and width. A seismogenic
    depth (from configuration file) is used to limit width and can modify centroid depth for
    large faults. To compute rupture plane corners from centroid, openquake geodetic functions
    are used, but note that errors accumulate for large distances.
    Args:
    - evconf: event configuration
    - calcconf: calculation configuration
    Return:
    - rctx: RuptureContext for rupture
    - faultplane: PlanarSurface for rupture
    '''
    seismogenic_depth = evconf['seisstruc']['seismogenic_depth']
    scalrel = getScalingRelation(calcconf)
    rctx = RuptureContext()
    rctx.rake = evconf['evmech']['rake']
    rctx.dip = evconf['evmech']['dip']
    rctx.mag = evconf['mag']
    rctx.hypo_depth = evconf['evloc']['hypo_depth']
    centroid_depth = evconf['evloc']['hypo_depth']
    diprad = radians(rctx.dip)
    if centroid_depth is None:
        centroid_depth = getKaklamanosCentroidDepth(rctx)
    flen = scalrel.get_median_length(rctx.mag, rctx.rake)
    fwid = scalrel.get_median_width(rctx.mag, rctx.rake)
    # Note that mag range extends width relation beyond validity, so fix assumes
    # len = wid when wid > len, i.e. aspect ratio 1.
    if fwid > flen: 
        rctx.width = flen
    else:
        rctx.width = fwid
    farea = scalrel.get_median_area(rctx.mag, rctx.rake)
    rctx, maxz, centroid_depth = adjustRCTXparams(rctx, seismogenic_depth, centroid_depth)
    logger.info(f'Fault dimensions: {farea}, {flen}, {rctx.width}, {flen*rctx.width}, {flen/rctx.width}')

    # Fault plane defined as 
    # lat: -0.5 fault_len to +0.5 fault_len; 
    # lon: -0.5 projected width to +0.5 projected width
    # given centroid lat, lon
    rctx.strike = evconf['evmech']['strike']
    e1_lon, e1_lat = geodetic.point_at(evconf['evloc']['centroid_lon'], 
            evconf['evloc']['centroid_lat'], evconf['evmech']['strike'], flen/2.) 
    e2_lon, e2_lat = geodetic.point_at(evconf['evloc']['centroid_lon'], 
            evconf['evloc']['centroid_lat'], evconf['evmech']['strike']+180., flen/2.) 
    xcorr = cos(diprad) * 0.5 * rctx.width
    tl_lon, tl_lat = geodetic.point_at(e2_lon, e2_lat, evconf['evmech']['strike']+270., xcorr) 
    tr_lon, tr_lat = geodetic.point_at(e1_lon, e1_lat, evconf['evmech']['strike']+270., xcorr) 
    bl_lon, bl_lat = geodetic.point_at(e2_lon, e2_lat, evconf['evmech']['strike']+90., xcorr) 
    br_lon, br_lat = geodetic.point_at(e1_lon, e1_lat, evconf['evmech']['strike']+90., xcorr) 
    logger.info('Fault corners: (%.4f, %.4f) (%.4f, %.4f) (%.4f, %.4f) (%.4f, %.4f)' % 
            (tl_lon, tl_lat, tr_lon, tr_lat, bl_lon, bl_lat, br_lon, br_lat))
    faultplane = PlanarSurface(strike=evconf['evmech']['strike'], dip=rctx.dip, 
            top_left=Point(depth=rctx.ztor, latitude=tl_lat, longitude=tl_lon),
            top_right=Point(depth=rctx.ztor, latitude=tr_lat, longitude=tr_lon),
            bottom_left=Point(depth=maxz, latitude=bl_lat, longitude=bl_lon),
            bottom_right=Point(depth=maxz, latitude=br_lat, longitude=br_lon))
    #fptest = PlanarSurface.from_hypocenter(Point(0., 0., rctx.hypo_depth), scalrel, rctx.mag, flen/rctx.width, 
    #    evconf['evmech']['strike'], rctx.dip, rctx.rake, ztor=None)
    logger.info(f'faultplane: {faultplane.get_surface_boundaries()}') 
    #logger.info(f'fptest: {fptest.get_surface_boundaries()}') 
    return [rctx], [faultplane]


def make_pga_lop(evconf, calcconf, rctx, faultplane, bPlots = False):
    '''
    Create Distance and Sites Contexts for a list of points
    Args:
    - evconf: event configuration
    - calcconf: calculation configuration
    - rctx: RuptureContext for fault
    - faultplane: PlanarSurface for fault
    - bPlots: boolean to control plot generation
    Return:
    - dctx: DistanceContext for the list of points
    - sctx: SitesContext for the list of points
    '''
    if not os.path.isfile(calcconf['points']['points_file']):
        return None, None
    # Expects a file with space separated fields: lat lon vs30 name
    fin = open(calcconf['points']['points_file'], 'r')
    lats = []
    lons = []
    lvs30 = []
    for l in fin:
        fs = l.split()
        lons.append(float(fs[1]))
        lats.append(float(fs[0]))
        lvs30.append(float(fs[2]))
    fin.close()
    mesh = Mesh(lons=np.asarray(lons), lats=np.asarray(lats), depths=np.zeros_like(np.asarray(lons)))
    dctx = DistancesContext()
    dctx.rjb = faultplane.get_joyner_boore_distance(mesh)
    if dctx.rjb.shape[0] != mesh.shape[0]:
        np.reshape(dctx.rjb, mesh.shape)
    ### JADEBUG 
#   # Noticed significant errors between computed rjb and distance between points for large fault
#   # sizes (several 100s km) and certain azimuths when siting in NZ. Not sure how general these error
#   # accumulations are, but avoiding these scenarios in testing for now!
#    tmp = zip(dctx.rjb, mesh.lons, mesh.lats)
#    for rjb,x,y in sorted(tmp):
#        print(x,y,rjb, geodetic.distance(x, y, 0., faultplane.top_right.longitude, 
#                faultplane.top_right.latitude, rctx.ztor))
    ### JADEBUG
    dctx.rjb_var = None
    dctx.rrup_var = None
    dctx.rvolc = np.zeros_like(dctx.rjb) # no correction for travel path in volcanic region
    rhyp = []
    for x, y in zip(lons, lats):
        rhyp.append(geodetic.distance(x, y, 0., evconf['evloc']['centroid_lon'], 
                evconf['evloc']['centroid_lat'], rctx.hypo_depth))
    dctx.rhypo = np.asarray(rhyp).reshape(len(lats))
    dctx.rx = faultplane.get_rx_distance(mesh)
    dctx.ry0 = faultplane.get_ry0_distance(mesh)
    dctx.rrup = faultplane.get_min_distance(mesh).reshape(mesh.shape)
    if bPlots:
        import matplotlib.pyplot as plt
        flat = []
        flon = []
        if isinstance(faultplane, PlanarSurface):
            for x in [faultplane.top_left, faultplane.top_right, faultplane.bottom_right, faultplane.bottom_left, faultplane.top_left]:
                flat.append(x.latitude)
                flon.append(x.longitude)
        else:
            top = faultplane.surface_nodes[0].nodes[0].nodes[0].nodes[0]
            bottom = faultplane.surface_nodes[0].nodes[-1].nodes[0].nodes[0]
            flon.extend([float(x) for x in top.to_str().split('[')[1].split(']')[0].split(',')[::3]])
            flat.extend([float(x) for x in top.to_str().split('[')[1].split(']')[0].split(',')[1::3]])
            flon.extend([float(x) for x in bottom.to_str().split('[')[1].split(']')[0].split(',')[::3]])
            flat.extend([float(x) for x in bottom.to_str().split('[')[1].split(']')[0].split(',')[1::3]])
        for v, lbl in zip([dctx.rrup, dctx.rjb, dctx.rx, dctx.ry0, dctx.rhypo], ['rrup', 'rjb', 'rx', 'ry0', 'rhypo']):
            cb = plt.scatter(lons, lats, c=[log10(x) if x>0 else 0 for x in v])
            plt.scatter(evconf['evloc']['centroid_lon'], evconf['evloc']['centroid_lat'], marker='*', s=80)
            plt.plot(flon, flat)
            plt.colorbar(cb)
            plt.savefig('%s_M%.1f_Azi%d.png' % (lbl, evconf['mag'], rctx.strike))
            plt.close()
    # --------------------------------------------------------------------------
    # Site context
    # --------------------------------------------------------------------------
    sctx = SitesContext()
    sctx.sids = np.arange(len(dctx.rjb))
    sctx.vs30 = np.asarray(lvs30)
    sctx.vs30measured = np.full_like(dctx.rjb, False, dtype="bool")
    sctx.z1pt0_ask14_cal = sites.Sites._z1pt0_from_vs30_ask14_cal(sctx.vs30)
    sctx.z1pt0_cy14_cal = sites.Sites._z1pt0_from_vs30_cy14_cal(sctx.vs30)
    sctx.z2pt5_cb14_cal = sites.Sites._z2pt5_from_vs30_cb14_cal(sctx.vs30) / 1000.0
    sctx.backarc = np.zeros_like(dctx.rjb, dtype=np.int16) # forearc/backarc is unknown
    return dctx, sctx


def make_pga_pt(evconf, calcconf, rctx, faultplane, bPlots = False):
    '''
    Create Distance and Sites Contexts for a point at fixed Rjb alongside fault and with
    fixed vs30. Other distance measures are computed.
    Args:
    - evconf: event configuration
    - calcconf: calculation configuration
    - rctx: RuptureContext for fault
    - faultplane: PlanarSurface for fault
    - bPlots: boolean to control plot generation
    Return:
    - dctx: DistanceContext for the list of points
    - sctx: SitesContext for the list of points
    '''
    # --------------------------------------------------------------------------
    # Distance context
    # --------------------------------------------------------------------------
    rjb = calcconf['pt']['rjb']
    dctx = DistancesContext()
    dctx.rjb = np.full((1, 1), rjb)
    dctx.rjb_var = None
    dctx.rrup_var = None
    dctx.rvolc = np.zeros_like(dctx.rjb) # no correction for travel path in volcanic region
    dctx.rhypo = np.full((1, 1), np.sqrt(pow(rjb, 2) + pow(rctx.hypo_depth, 2)))
    dctx.rx = dctx.rjb
    dctx.ry0 = np.zeros_like(dctx.rjb)
    dctx.rrup = faultplane.get_min_distance(mesh).reshape(mesh.shape)
    # --------------------------------------------------------------------------
    # Site context
    # --------------------------------------------------------------------------
    sctx = SitesContext()
    sctx.sids = np.arange(1).reshape(1,1)
    sctx.vs30 = np.ones_like(dctx.rjb) * calcconf['pt']['vs30']
    sctx.vs30measured = np.full_like(dctx.rjb, False, dtype="bool")
    sctx.z1pt0_ask14_cal = sites.Sites._z1pt0_from_vs30_ask14_cal(sctx.vs30)
    sctx.z1pt0_cy14_cal = sites.Sites._z1pt0_from_vs30_cy14_cal(sctx.vs30)
    sctx.z2pt5_cb14_cal = sites.Sites._z2pt5_from_vs30_cb14_cal(sctx.vs30) / 1000.0
    sctx.backarc = np.zeros_like(dctx.rjb, dtype=np.int16) # forearc/backarc is unknown
    return dctx, sctx


def mesh_from_bb(calcconf, minlat, maxlat, minlon, maxlon):
    '''
    Make a mesh for use as distance context based on extent
    Args:
    - calcconf
    - minlat: southern latitude
    - maxlat: northern latitude
    - minlon: western longitude
    - maxlon: eastern longitude
    Return:
    - mesh
    '''
    xlim = geodetic.distance(minlon, maxlat, 0., maxlon, maxlat, 0.)
    ylim = geodetic.distance(minlon, minlat, 0., minlon, maxlat, 0.)
    dkm = calcconf['grid']['griddkm']
    nx = ceil(xlim/dkm)
    xdist = nx*dkm
    inlons, lats, depths = geodetic.npoints_towards(minlon, maxlat, 0., 90., xdist, 0., nx+1)
    ny = ceil(ylim/dkm)
    ydist = ny*dkm
    lons, inlats, depths = geodetic.npoints_towards(minlon, minlat, 0., 0., ydist, 0., ny+1)
    lons, lats = np.meshgrid(inlons, inlats)
    mesh = Mesh(lons=lons.ravel(), lats=lats.ravel(), depths=None)
    return mesh


def make_mesh(evconf, calcconf, faultplane):
    '''
    Make a mesh for use as distance context
    Args:
    - evconf: event configuration
    - calcconf: calc configuration
    - faultplane: complex fault plane
    Return:
    - mesh
    '''
    # ----
    # Empirical estimate for maximum distance needed in template grid
    # ----
    maxdist = max([evconf['mag']*115. - 400., evconf['mag']*40 - 95])
    flen = faultplane.get_area()/faultplane.get_width()
    fwid = faultplane.get_width()
    bb = faultplane.get_bounding_box()
    clon = bb.west + ((bb.east - bb.west) / 2.)
    clat = bb.south + ((bb.north - bb.south) / 2.)
    xlim = geodetic.distance(clon, clat, 0., bb.east, clat, 0.) + maxdist
    ylim = geodetic.distance(clon, clat, 0., clon, bb.north, 0.) + maxdist
    logger.info(r'Mag %.1f  Max dist %.4f Rup len %.4f Rup wid %.4f xlim %.4f ylim %.4f ztor %.2f' % \
        (evconf['mag'], maxdist, flen, fwid, xlim, ylim, faultplane.get_top_edge_depth()))
    dkm = calcconf['grid']['griddkm']
    nx = ceil(xlim/dkm)
    xdist = nx*dkm
    lons, lats, depths = geodetic.npoints_towards(clon, clat, 0., 270., xdist, 0., nx+1)
    lonsr = np.flip(lons)
    lons, lats, depths = geodetic.npoints_towards(clon, clat, 0., 90., xdist, 0., nx+1)
    inlons = np.concatenate([lonsr[:-1], lons])
    ny = ceil(ylim/dkm)
    ydist = ny*dkm
    lons, lats, depths = geodetic.npoints_towards(clon, clat, 0., 180., ydist, 0., ny+1)
    latsr = np.flip(lats)
    lons, lats, depths = geodetic.npoints_towards(clon, clat, 0., 0., ydist, 0., ny+1)
    inlats = np.concatenate([latsr[:-1], lats])
    lons, lats = np.meshgrid(inlons, inlats)
    mesh = Mesh(lons=lons.ravel(), lats=lats.ravel(), depths=None)
    return mesh


def make_pga_grid(evconf, calcconf, rctx, faultplane, bPlots = False):
    '''
    Create Distance and Sites Contexts for a grid of points, fixed vs30
    Args:
    - evconf: event configuration
    - calcconf: calculation configuration
    - rctx: RuptureContext for fault
    - faultplane: PlanarSurface for fault
    - bPlots: boolean to control plot generation
    Return:
    - dctx: DistanceContext for the list of points
    - sctx: SitesContext for the list of points
    '''
    if 'mesh' not in calcconf:
        mesh = make_mesh(evconf, calcconf, faultplane)
    else:
        mesh = calcconf['mesh']
    # --------------------------------------------------------------------------
    # Distance context
    # --------------------------------------------------------------------------
    dctx = DistancesContext()
    dctx.rjb = faultplane.get_joyner_boore_distance(mesh)
    dctx.rjb_var = None
    dctx.rrup_var = None
    dctx.rvolc = np.zeros_like(dctx.rjb) # no correction for travel path in volcanic region
    rhyp = []
    bb = faultplane.get_bounding_box()
    clon = bb.west + ((bb.east - bb.west) / 2.)
    clat = bb.south + ((bb.north - bb.south) / 2.)
    for y, x in zip(mesh.lats, mesh.lons):
        rhyp.append(geodetic.distance(x, y, 0., clon, clat, rctx.hypo_depth))
    dctx.rhypo = np.asarray(rhyp)
    dctx.rx = faultplane.get_rx_distance(mesh)
    dctx.ry0 = faultplane.get_ry0_distance(mesh)
    dctx.rrup = faultplane.get_min_distance(mesh)

    nr = len(set(mesh.lats))
    nc = len(set(mesh.lons))
    dctx.rjb = np.reshape(dctx.rjb, (nr, nc))
    dctx.rvolc = np.reshape(dctx.rvolc, (nr, nc))
    dctx.rhypo = np.reshape(dctx.rhypo, (nr, nc))
    dctx.rx = np.reshape(dctx.rx, (nr, nc))
    dctx.ry0 = np.reshape(dctx.ry0, (nr, nc))
    dctx.rrup = np.reshape(dctx.rrup, (nr, nc))
    if bPlots:
        import matplotlib.pyplot as plt
        for v, lbl in zip([dctx.rrup, dctx.rjb, dctx.rx, dctx.ry0, dctx.rhypo], ['rrup', 'rjb', 'rx', 'ry0', 'rhypo']):
            cb = plt.imshow(v, origin='lower')
            plt.colorbar(cb)
            plt.savefig('%s_M%.1f.png' % (lbl, evconf['mag']))
            plt.close()
    # --------------------------------------------------------------------------
    # Site context
    # --------------------------------------------------------------------------
    sctx = SitesContext()
    sctx.sids = np.arange(nc*nr).reshape(dctx.rjb.shape)
    #sctx.sids = np.arange(len(dctx.rjb))
    sctx.vs30 = np.ones_like(dctx.rjb) * calcconf['grid']['vs30']
    sctx.vs30measured = np.full_like(dctx.rjb, False, dtype="bool")
    sctx.z1pt0_ask14_cal = sites.Sites._z1pt0_from_vs30_ask14_cal(sctx.vs30)
    sctx.z1pt0_cy14_cal = sites.Sites._z1pt0_from_vs30_cy14_cal(sctx.vs30)
    sctx.z2pt5_cb14_cal = sites.Sites._z2pt5_from_vs30_cb14_cal(sctx.vs30) / 1000.0
    sctx.backarc = np.zeros_like(dctx.rjb, dtype=np.int16) # forearc/backarc is unknown
    return dctx, sctx


def createRuptures(evconf, calcconf):
    '''
    Create rupture based on geometry or basic fault parameters
    Args:
    - 
    Return:
    - 
    '''
    # --------------------------------------------------------------------------
    # Set rupture basics
    # --------------------------------------------------------------------------
    if 'geometry' in evconf['evmech'] and os.path.isfile(evconf['evmech']['geometry']) \
            and evconf['evmech']['geometry'].split('.')[-1] == 'json':
        logger.info('Computing for fault geometry, so ensure no magnitude loop')
        calcconf['magrange']['magmax'] = calcconf['magrange']['magmin']
    rups = {}
    for mag in np.arange(calcconf['magrange']['magmin'], calcconf['magrange']['magmax'] + \
            calcconf['magrange']['magstep']/2, calcconf['magrange']['magstep']):
        # --------------------------------------------------------------------------
        # Rupture context
        # --------------------------------------------------------------------------
        if 'geometry' in evconf['evmech'] and os.path.isfile(evconf['evmech']['geometry']):
            if evconf['evmech']['geometry'].split('.')[-1] == 'json':
                evconf, l_rctx, l_faultplane = importRuptureContext(evconf)
                mag = evconf['mag']
            elif evconf['evmech']['geometry'].split('.')[-1] == 'geojson': 
                evconf['mag'] = mag
                l_rctx, l_faultplane = createSubFaultRuptureContexts(evconf, calcconf)        
                if l_rctx is None:
                    return rups
        else:
            evconf['mag'] = mag
            l_rctx, l_faultplane = createRuptureContext(evconf, calcconf)
        rups[mag] = {}
        for rctx, faultplane in zip(l_rctx, l_faultplane): 
            centroid = faultplane.get_middle_point()
            rups[mag][(centroid.latitude, centroid.longitude)] = [rctx, faultplane]
    return rups


def computeGM(gmpeconf, evconf, calcconf):
    '''
    Compute ground motion using openquake functionality. Use configuration files to define
    the fault plane and earthquake parameters, points at which to calculate ground motion,
    GMPEs and calculation settings. Ground motion computed is PGA in log10(cm/s/s) for
    greater of two horizontals.
    Args:
    - gmpeconf: GMPE configuration
    - evconf: event configuration
    - calcconf: calculation configuration
    Return:
    - gm: dictionary with magnitude as key and values: median PGA in log10(cm/s/s), 
    - faultplane (PlanarSurface) object
    '''
    rups = createRuptures(evconf, calcconf)
    template_sets = None
    if 'fault-specific' in calcconf:
        template_sets = createTemplateSets([(m, clat, clon) for m in rups for (clat, clon) in rups[m]], calcconf)
        if 'grid' in calcconf and calcconf['grid']['compute']:
            for i in sorted(template_sets):
                l_mesh_lats = []
                l_mesh_lons = []
                for mag in sorted(template_sets[i]):
                    for [clat, clon] in template_sets[i][mag]:
                        evconf['mag'] = mag
                        bb = make_mesh(evconf, calcconf, rups[mag][(clat, clon)][1]).get_convex_hull().get_bbox()
                        l_mesh_lats.append(bb[1])
                        l_mesh_lats.append(bb[3])
                        l_mesh_lons.append(bb[0])
                        l_mesh_lons.append(bb[2])
                bb = utils.get_spherical_bounding_box(l_mesh_lons, l_mesh_lats)
                template_sets[i]['mesh'] = mesh_from_bb(calcconf, bb.south, bb.north, bb.west, bb.east)

    # --------------------------------------------------------------------------
    # Get multigmpe from config
    # --------------------------------------------------------------------------
    mgmpe = MultiGMPE.__from_config__(gmpeconf)
    IMT = imt.PGA() # will be in units of "g"
    gm = {}
    tset = None
    for mag in sorted(rups):
        for centroid_lat, centroid_lon in sorted(rups[mag]): 
            evconf['mag'] = mag
            rctx, faultplane = rups[mag][(centroid_lat, centroid_lon)]
            # --------------------------------------------------------------------------
            # Distance and Source contexts
            # --------------------------------------------------------------------------
            if 'grid' in calcconf and calcconf['grid']['compute']:
                if template_sets is not None:
                    for i in template_sets:
                        if mag in template_sets[i] and [centroid_lat, centroid_lon] in template_sets[i][mag]:
                            tset = i
                            calcconf['mesh'] = template_sets[i]['mesh']
                dctx, sctx = make_pga_grid(evconf, calcconf, rctx, faultplane, bPlots = calcconf['plots'])
            if 'points' in calcconf and calcconf['points']['compute']:
                dctx, sctx = make_pga_lop(evconf, calcconf, rctx, faultplane, bPlots = calcconf['plots'])
            if 'pt' in calcconf and calcconf['pt']['compute']:
                dctx, sctx = make_pga_pt(evconf, calcconf, rctx, faultplane, bPlots = calcconf['plots'])
            # --------------------------------------------------------------------------
            # Compute ground motion
            # --------------------------------------------------------------------------
            lmean_mgmpe, lmean_sd = mgmpe.get_mean_and_stddevs(sctx, rctx, dctx, IMT, [const.StdDev.TOTAL])
            # US ShakeAlert measures largest of two horizontals (FFD2)
            # Seiscomp measures peak of real-time root or sum or squares of two horizontals
            # Most GMPEs use geometric mean measures (GM, GMRotD50 etc.) which are smaller
            # than either of the two above measures. Note that no conversion to the latter 
            # exists, but to the former is available with either class (Beyer & Bommer or
            # Boore & Kishida). Since the former is smaller than the latter, converting to
            # the former gets us closer...
            if mgmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT != const.IMC.GREATER_OF_TWO_HORIZONTAL:
                logger.info('Converting component')
                bk17 = BooreKishida2017(
                    mgmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT, const.IMC.GREATER_OF_TWO_HORIZONTAL
                )
                lmean_mgmpe = bk17.convertAmpsOnce(IMT, lmean_mgmpe, dctx.rrup, rctx.mag)
            # Method get_mean_and_stddevs() of actual GMPE implementations is supposed to return the 
            # mean value as a natural logarithm of intensity.
            if mag not in gm:
                gm[mag] = {}
            gm[mag][(centroid_lat, centroid_lon)] = [lng2cm(lmean_mgmpe), faultplane]
    return gm, evconf, template_sets


def createTemplateSets(data, calcconf):
    '''
    Split a set of magnitudes and centroid locations in sets based on the allowable number using both
    max and min per magnitude
    Args:
    - data: list of [mag, centroid_lat, centroid_lat]
    - calcconf: configuration 
    Return:
    - 
    '''
    from math import ceil
    maxperset = calcconf['fault-specific']['maxperset']
    minperset = calcconf['fault-specific']['minperset']
    # Convert input list into dictionary with key as magnitude and centroid as value
    templates = {}
    for d in data:
        if d[0] not in templates:
            templates[d[0]] = [[float(d[1]), float(d[2])]]
        else:
            templates[d[0]].append([float(d[1]), float(d[2])])
    n_tot_sets = 0
    while len([t for m in templates for t in templates[m] if len(t) < 3]) > 0:
        logger.info(f'Templates without a set: {len([t for m in templates for t in templates[m] if len(t) < 3])}')
        # Minimum magnitude for this set
        minmag = None
        for mag in templates:
            t = [t for t in templates[mag] if len(t) == 2]
            if len(t) > 0:
                minmag = mag
                break
        if minmag is None:
            return templates, 0
        # Determine the number of minimum magnitude templates for each set
        sets = []
        ntot = len(templates[minmag])
        nsets = round(ntot/maxperset)
        if nsets == 0:
            nsets = 1
        basen = ntot//nsets
        left = ntot%basen
        extras = ceil((left)/nsets)
        logger.info(f'For minmag: {minmag}, # templates: {ntot}, # sets: {nsets}, min # in set {basen}, remainder: {left}, extras per set: {extras}')
        # Determine the start/stop indices for each set for the minimum magnitude templates
        rsum = 0
        for i in range(nsets):
            rsum += basen
            if left > 0:
                rsum += extras
                left -= extras
            sets.append(rsum)
        logger.info(f'Minmag template end indices {sets}')
        # Add templates from minimum magnitude
        setbb = []
        start = 0
        for i, end in enumerate(sets):
            tset = templates[minmag][start:end]
            for t in tset:
                t.append(i+n_tot_sets)
            sclat = tset[0][0]
            eclat = tset[-1][0]
            sclon = tset[0][1]
            eclon = tset[-1][1]
            setbb.append([sclat, eclat, sclon, eclon])
            start = end
        # Add templates from minimum magnitude
        maxmag = None
        tol = 0.01
        for i, bb in enumerate(setbb):
            if i == 0: 
                slat = bb[0]
                slon = bb[2]
            else:
                slat = (setbb[i-1][1] + bb[0])/2
                slon = (setbb[i-1][3] + bb[2])/2
            if i == len(setbb)-1:
                elat = bb[1]
                elon = bb[3]
            else:
                elat = (setbb[i+1][0] + bb[1])/2
                elon = (setbb[i+1][2] + bb[3])/2
            minclat = min([slat, elat]) - tol
            maxclat = max([slat, elat]) + tol
            minclon = min([slon, elon]) - tol
            maxclon = max([slon, elon]) + tol
            for mag in [m for m in templates if m > minmag]:
                if maxmag is not None and mag >= maxmag:
                    break
                # get all templates for this mag with centroids in min/max range
                n_tset = [t for t in templates[mag] if \
                        t[0] >= minclat and t[0] <= maxclat and \
                        t[1] >= minclon and t[1] <= maxclon]
                if maxmag is None and len(templates[mag]) > minperset and len(n_tset) < minperset:
                    logger.info(f'Set maxmag: {mag}')
                    maxmag = mag
                    break
                for t in n_tset:
                    t.append(i+n_tot_sets)
        n_tot_sets += len(sets)
    template_sets = {}
    for i in range(n_tot_sets):
        template_sets[i] = {}
        for mag in templates:
            templs = [t[0:2] for t in templates[mag] if t[-1] == i]
            if len(templs) == 0:
                continue
            template_sets[i][mag] = templs
    return template_sets


##################################################################################################
# ORIGINAL VERSION OF DCTX THAT EXPLOITED SYMMETRY FOR A VERTICAL FAULT
#        dctx = DistancesContext()
#        rjb = []
#        rhypo = []
#        ry0 = []
#        rx = []
#        rrup = []
#        for x in np.arange(0., maxdist+(dkm/2.), dkm):
#            topside = 0
#            for y in np.arange(0., maxdist+(dkm/2.), dkm):
#                topside += 1
#                rjb.append(sqrt(pow(x,2) + pow(y,2)))
#                rhypo.append(sqrt(pow(x,2.) + pow(y+(flen/2),2.) + pow(rctx.centroid_depth,2.)))
#                ry0.append(y)
#                rx.append(x + xcorr_fw)
#        for x in np.arange(0., maxdist+(dkm/2.), dkm):
#            bottomside = 0
#            for y in np.arange(dkm - (flen % dkm), maxdist+(dkm/2.), dkm):
#                bottomside += 1
#                rjb.append(sqrt(pow(x,2.) + pow(y,2.)))
#                rhypo.append(sqrt(pow(x,2.) + pow(y+(flen/2.),2.) + pow(rctx.centroid_depth,2.)))
#                ry0.append(y)
#                rx.append(x + xcorr_fw)
#        if rctx.ztor == 0.:
#            rrup = rjb
#        else:
#            rrup = [sqrt(pow(d,2.)+pow(rctx.ztor,2.)) for d in rjb]
#        dctx.rrup = np.asarray(rrup)
#        dctx.rjb = np.asarray(rjb)
#        dctx.rjb_var = None
#        dctx.rrup_var = None
#        dctx.rhypo = np.asarray(rhypo)
#        dctx.rx = np.asarray(rx) 
#        dctx.ry0 = np.asarray(ry0)

#       This exploits symmetry and is only applicable for vertical strike slip faults
#        print(pow(10., min(lmean_mgmpe[:topside])), pow(10., max(lmean_mgmpe[:topside])))
#        templ_top = np.ones((topside, (topside*2)-1), np.float64)
#        i = 0
#        for x in range(topside):
#            for y in range(x, topside):
#                templ_top[topside-1-y][topside-1+x] = lmean_mgmpe[i]
#                templ_top[topside-1-y][topside-1-x] = lmean_mgmpe[i]
#                i += 1
#        if flen > dkm: # duplicate the zeroth data along the length of the fault
#            templ_mid = np.ones((floor(flen/dkm), (topside*2)-1))
#            for x in range(topside):
#                for y in range(floor(flen/dkm)):
#                    templ_mid[y][topside-1-x] = lmean_mgmpe[x]
#                    templ_mid[y][topside-1+x] = lmean_mgmpe[x]
#        templ_bottom = np.ones((bottomside, (topside*2)-1), np.float64)
#        for x in range(topside):
#            for y in range(bottomside):
#                templ_bottom[y][topside-1+x] = lmean_mgmpe[i]
#                templ_bottom[y][topside-1-x] = lmean_mgmpe[i]
#                i += 1
#        if flen > dkm: # duplicate the zeroth data along the length of the fault
#            templ_all = np.vstack((templ_top, templ_mid, templ_bottom))
#        else:
#            templ_all = np.vstack((templ_top, templ_bottom))

###
### Code to compute rrup geometrically
# Originally for long faults, the get_min_distance function 
# appears to only use the corners to 
# compute horizontal distance, resulting in Rrup >> actual
###
#    if rctx.dip == 90.:
#        dctx.rrup = np.sqrt(np.square(dctx.rjb) + pow(rctx.ztor,2))
#    else:
#        diprad = radians(rctx.dip)
#        # Very frustrating that we have to do this!
#        # Eq 14-20 Kaklamanos et al. (2011)
#        Acond = rctx.ztor * tan(diprad)
#        A = pow(rctx.ztor,2)
#        Arrupp = np.where(dctx.rx < Acond, np.sqrt(np.square(dctx.rx) + A), 0)
#        Ccond = Acond + (rctx.width / cos(diprad))
#        A = rctx.width * cos(diprad)
#        B = pow(rctx.ztor + (rctx.width * sin(diprad)), 2)
#        Crrupp = np.where(dctx.rx > Ccond, np.sqrt(np.square(dctx.rx - A) + B), 0)
#        A = sin(diprad)
#        B = rctx.ztor * cos(diprad)
#        Brrupp = np.where(np.logical_and(dctx.rx >= Acond, dctx.rx <= Ccond), (dctx.rx * A) + B, 0)
#        rrupp = Arrupp + Brrupp + Crrupp
#        dctx.rrup = np.sqrt(np.square(rrupp) + np.square(dctx.ry0))
