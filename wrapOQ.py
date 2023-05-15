#!/usr/bin/env python

#############
# N.B. This script is based heavily on shakemap_src/tests/shakelib/multigmpe_test.py
#############

# stdlib imports
import copy
import os
import sys
import logging
from math import log10, sqrt, exp, floor, sin, cos, tan, radians, ceil
import numpy as np

# openquake imports
from openquake.hazardlib import imt, const
from openquake.hazardlib.gsim.base import RuptureContext
from openquake.hazardlib.gsim.base import DistancesContext
from openquake.hazardlib.gsim.base import SitesContext
from openquake.hazardlib.geo import geodetic, Point, PlanarSurface, Mesh, RectangularMesh

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
        Numpy array of PGA values in ln(g)
    Return:
        Numpy array of PGA values in log10(cm/s/s)
    '''
    inx = np.log10(np.exp(inx)) + log10(980.665)
    return inx


def getScalingRelation(conf):
    '''
    Returns the ScalingRelation class using configuration file input.
    Allowing for extension to other scaling relations. However, all will require an Ext class to 
    provide length and width rather than area.
    Args:
        conf: configuration dictionary
    Return:
        ScalingRelation class object
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
        return Blase2010.Blase2010()
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
        fname: filename
    Return:
        Dictionary with contents of configuration file
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


def createRuptureContext(evconf, calcconf):
    '''
    Creates RuptureContext based on configuration file input for a simple, single planar surface
    defined by centroid latitude and longitude, strike, dip and rake. Magnitude and scaling 
    relation (from configuration file) are used to determine fault length and width. A seismogenic
    depth (from configuration file) is used to limit width and can modify centroid depth for
    large faults. To compute rupture plane corners from centroid, openquake geodetic functions
    are used, but note that errors accumulate for large distances.
    Args:
        evconf: event configuration
        calcconf: calculation configuration
    Return:
        rctx: RuptureContext for rupture
        faultplane: PlanarSurface for rupture
        xcorr: fault half width distance projected at surface
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
    # ----
    # If no depth is given, use Kaklamanos el al. (2010) to set
    # ----
    if centroid_depth is None:
        if rctx.rake is None:
            centroid_depth = 7.08 + (0.61 * rctx.mag) # Kaklamanos eq.4 general
        if (-45 <= rctx.rake <= 45) or (rctx.rake >= 135) or (rctx.rake <= -135):
            centroid_depth = 5.63 + (0.68 * rctx.mag) # Kaklamanos eq.4 SS
        else:
            centroid_depth = 11.24 - (0.2 * rctx.mag) # Kaklamanos eq.4 non-SS
    flen = scalrel.get_median_length(rctx.mag, rctx.rake)
    fwid = scalrel.get_median_width(rctx.mag, rctx.rake)
    # Note that mag range extends width relation beyond validity, so fix assumes
    # len = wid when wid > len
    if fwid > flen: 
        rctx.width = flen
    else:
        rctx.width = fwid
    # For FinDer, centroid must be kept at center of template
    # Use initial hypo depth as centroid depth unless width requires us to move it
    # Hypo depth empirically found to be at about 60% depth, but simpler to treat as 50% here
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
    logging.info('Fault corners: (%.4f, %.4f) (%.4f, %.4f) (%.4f, %.4f) (%.4f, %.4f)' % 
            (tl_lon, tl_lat, tr_lon, tr_lat, bl_lon, bl_lat, br_lon, br_lat))
    faultplane = PlanarSurface(strike=evconf['evmech']['strike'], dip=rctx.dip, 
            top_left=Point(depth=rctx.ztor, latitude=tl_lat, longitude=tl_lon),
            top_right=Point(depth=rctx.ztor, latitude=tr_lat, longitude=tr_lon),
            bottom_left=Point(depth=maxz, latitude=bl_lat, longitude=bl_lon),
            bottom_right=Point(depth=maxz, latitude=br_lat, longitude=br_lon))
    return rctx, faultplane, xcorr


def make_pga_lop(evconf, calcconf, rctx, faultplane, xcorr, bPlots = False):
    '''
    Create Distance and Sites Contexts for a list of points
    Args:
        evconf: event configuration
        calcconf: calculation configuration
        rctx: RuptureContext for fault
        faultplane: PlanarSurface for fault
        xcorr: fault half width distance projected at surface
        bPlots: boolean to control plot generation
    Return:
        dctx: DistanceContext for the list of points
        sctx: SitesContext for the list of points
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
    if rctx.dip == 90.:
        dctx.rrup = np.sqrt(np.square(dctx.rjb) + pow(rctx.ztor,2))
    else:
        diprad = radians(rctx.dip)
        # Very frustrating that we have to do this!
        # Eq 14-20 Kaklamanos et al. (2011)
        Acond = rctx.ztor * tan(diprad)
        A = pow(rctx.ztor,2)
        Arrupp = np.where(dctx.rx < Acond, np.sqrt(np.square(dctx.rx) + A), 0)
        Ccond = Acond + (rctx.width / cos(diprad))
        A = rctx.width * cos(diprad)
        B = pow(rctx.ztor + (rctx.width * sin(diprad)), 2)
        Crrupp = np.where(dctx.rx > Ccond, np.sqrt(np.square(dctx.rx - A) + B), 0)
        A = sin(diprad)
        B = rctx.ztor * cos(diprad)
        Brrupp = np.where(np.logical_and(dctx.rx >= Acond, dctx.rx <= Ccond), (dctx.rx * A) + B, 0)
        rrupp = Arrupp + Brrupp + Crrupp
        dctx.rrup = np.sqrt(np.square(rrupp) + np.square(dctx.ry0))
    if bPlots:
        import matplotlib.pyplot as plt
        flat = []
        flon = []
        for x in [faultplane.top_left, faultplane.top_right, faultplane.bottom_right, faultplane.bottom_left, faultplane.top_left]:
            flat.append(x.latitude)
            flon.append(x.longitude)
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


def make_pga_pt(evconf, calcconf, rctx, faultplane, xcorr, bPlots = False):
    '''
    Create Distance and Sites Contexts for a point at fixed Rjb alongside fault and with
    fixed vs30. Other distance measures are computed.
    Args:
        evconf: event configuration
        calcconf: calculation configuration
        rctx: RuptureContext for fault
        faultplane: PlanarSurface for fault
        xcorr: fault half width projected at surface
        bPlots: boolean to control plot generation
    Return:
        dctx: DistanceContext for the list of points
        sctx: SitesContext for the list of points
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
    # For long faults, the get_min_distance function appears to only use the corners to 
    # compute horizontal distance, resulting in Rrup >> actual
    # If future libs fix this issue, reinstate following line:
    #dctx.rrup = faultplane.get_min_distance(mesh).reshape(mesh.shape)
    if rctx.dip == 90.:
        dctx.rrup = np.sqrt(np.square(dctx.rjb) + pow(rctx.ztor,2))
    else:
        diprad = radians(rctx.dip)
        # Very frustrating that we have to do this!
        # Eq 14-20 Kaklamanos et al. (2011)
        Acond = rctx.ztor * tan(diprad)
        A = pow(rctx.ztor,2)
        Arrupp = np.where(dctx.rx < Acond, np.sqrt(np.square(dctx.rx) + A), 0)
        Ccond = Acond + (rctx.width / cos(diprad))
        A = rctx.width * cos(diprad)
        B = pow(rctx.ztor + (rctx.width * sin(diprad)), 2)
        Crrupp = np.where(dctx.rx > Ccond, np.sqrt(np.square(dctx.rx - A) + B), 0)
        A = sin(diprad)
        B = rctx.ztor * cos(diprad)
        Brrupp = np.where(np.logical_and(dctx.rx >= Acond, dctx.rx <= Ccond), (dctx.rx * A) + B, 0)
        rrupp = Arrupp + Brrupp + Crrupp
        dctx.rrup = np.sqrt(np.square(rrupp) + np.square(dctx.ry0))
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


def make_pga_grid(evconf, calcconf, rctx, faultplane, xcorr, bPlots = False):
    '''
    Create Distance and Sites Contexts for a grid of points, fixed vs30
    Args:
        evconf: event configuration
        calcconf: calculation configuration
        rctx: RuptureContext for fault
        faultplane: PlanarSurface for fault
        xcorr: fault half width projected at surface
        bPlots: boolean to control plot generation
    Return:
        dctx: DistanceContext for the list of points
        sctx: SitesContext for the list of points
    '''
    # ----
    # Empirical estimate for maximum distance needed in template grid
    # ----
    flen = faultplane.get_area()/faultplane.get_width()
    fwid = faultplane.get_width()
    maxz = faultplane.bottom_left.depth
    maxdist = max([evconf['mag']*115. - 400., evconf['mag']*40 - 95])
    logger.info(r'Mag %.1f  Max dist %.4f Rup len %.4f Rup wid %.4f/%.4f ztor %.2f zmax %.2f' % \
        (evconf['mag'], maxdist, flen, fwid, rctx.width, rctx.ztor, maxz))
    xlim = xcorr + maxdist
    ylim = (flen/2.) + maxdist
    # --------------------------------------------------------------------------
    # Distance context
    # --------------------------------------------------------------------------
    dkm = calcconf['grid']['griddkm']
    nx = ceil(xlim/dkm)
    xdist = nx*dkm
    lons, lats, depths = geodetic.npoints_towards(evconf['evloc']['centroid_lon'], 
            evconf['evloc']['centroid_lat'], 0., 270., xdist, 0., nx+1)
    lonsr = np.flip(lons)
    lons, lats, depths = geodetic.npoints_towards(evconf['evloc']['centroid_lon'], 
            evconf['evloc']['centroid_lat'], 0., 90., xdist, 0., nx+1)
    inlons = np.concatenate([lonsr[:-1], lons])
    ny = ceil(ylim/dkm)
    ydist = ny*dkm
    lons, lats, depths = geodetic.npoints_towards(evconf['evloc']['centroid_lon'], 
            evconf['evloc']['centroid_lat'], 0., 180., ydist, 0., ny+1)
    latsr = np.flip(lats)
    lons, lats, depths = geodetic.npoints_towards(evconf['evloc']['centroid_lon'], 
            evconf['evloc']['centroid_lat'], 0., 0., ydist, 0., ny+1)
    inlats = np.concatenate([latsr[:-1], lats])
    lons, lats = np.meshgrid(inlons, inlats)
    mesh = RectangularMesh(lons=lons, lats=lats, depths=None)
    dctx = DistancesContext()
    dctx.rjb = faultplane.get_joyner_boore_distance(mesh)
    if dctx.rjb.shape[0] != mesh.shape[0]:
        dctx.rjb = np.reshape(dctx.rjb, mesh.shape)
    dctx.rjb_var = None
    dctx.rrup_var = None
    dctx.rvolc = np.zeros_like(dctx.rjb) # no correction for travel path in volcanic region
    rhyp = []
    for y in inlats:
        for x in inlons:
            rhyp.append(geodetic.distance(x, y, 0., evconf['evloc']['centroid_lon'], 
                evconf['evloc']['centroid_lat'], rctx.hypo_depth))
    dctx.rhypo = np.asarray(rhyp).reshape(len(inlats), len(inlons))
    dctx.rx = faultplane.get_rx_distance(mesh)
    if dctx.rx.shape[0] != mesh.shape[0]:
        dctx.rx = np.reshape(dctx.rx, mesh.shape)
    dctx.ry0 = faultplane.get_ry0_distance(mesh)
    if dctx.ry0.shape[0] != mesh.shape[0]:
        dctx.ry0 = np.reshape(dctx.ry0, mesh.shape)
    # For long faults, the get_min_distance function appears to only use the corners to 
    # compute horizontal distance, resulting in Rrup >> actual
    # If future libs fix this issue, reinstate following line:
    #dctx.rrup = faultplane.get_min_distance(mesh).reshape(mesh.shape)
    if rctx.dip == 90.:
        dctx.rrup = np.sqrt(np.square(dctx.rjb) + pow(rctx.ztor,2))
    else:
        diprad = radians(rctx.dip)
        # Very frustrating that we have to do this!
        # Eq 14-20 Kaklamanos et al. (2011)
        Acond = rctx.ztor * tan(diprad)
        A = pow(rctx.ztor,2)
        Arrupp = np.where(dctx.rx < Acond, np.sqrt(np.square(dctx.rx) + A), 0)
        Ccond = Acond + (rctx.width / cos(diprad))
        A = rctx.width * cos(diprad)
        B = pow(rctx.ztor + (rctx.width * sin(diprad)), 2)
        Crrupp = np.where(dctx.rx > Ccond, np.sqrt(np.square(dctx.rx - A) + B), 0)
        A = sin(diprad)
        B = rctx.ztor * cos(diprad)
        Brrupp = np.where(np.logical_and(dctx.rx >= Acond, dctx.rx <= Ccond), (dctx.rx * A) + B, 0)
        rrupp = Arrupp + Brrupp + Crrupp
        dctx.rrup = np.sqrt(np.square(rrupp) + np.square(dctx.ry0))
    if bPlots:
        import matplotlib.pyplot as plt
        for v, lbl in zip([dctx.rrup, dctx.rjb, dctx.rx, dctx.ry0, dctx.rhypo], ['rrup', 'rjb', 'rx', 'ry0', 'rhypo']):
            cb = plt.imshow(v)
            plt.colorbar(cb)
            plt.savefig('%s_M%.1f.png' % (lbl, evconf['mag']))
            plt.close()
    # --------------------------------------------------------------------------
    # Site context
    # --------------------------------------------------------------------------
    sctx = SitesContext()
    sctx.sids = np.arange(len(inlons)*len(inlats)).reshape(dctx.rjb.shape)
    #sctx.sids = np.arange(len(dctx.rjb))
    sctx.vs30 = np.ones_like(dctx.rjb) * calcconf['grid']['vs30']
    sctx.vs30measured = np.full_like(dctx.rjb, False, dtype="bool")
    sctx.z1pt0_ask14_cal = sites.Sites._z1pt0_from_vs30_ask14_cal(sctx.vs30)
    sctx.z1pt0_cy14_cal = sites.Sites._z1pt0_from_vs30_cy14_cal(sctx.vs30)
    sctx.z2pt5_cb14_cal = sites.Sites._z2pt5_from_vs30_cb14_cal(sctx.vs30) / 1000.0
    sctx.backarc = np.zeros_like(dctx.rjb, dtype=np.int16) # forearc/backarc is unknown
    return dctx, sctx


def computeGM(gmpeconf, evconf, calcconf):
    '''
    Compute ground motion using openquake functionality. Use configuration files to define
    the fault plane and earthquake parameters, points at which to calculate ground motion,
    GMPEs and calculation settings. Ground motion computed is PGA in log10(cm/s/s) for
    greater of two horizontals.
    Args:
        gmpeconf: GMPE configuration
        evconf: event configuration
        calcconf: calculation configuration
    Return:
        gm: dictionary with magnitude as key and values: median PGA in log10(cm/s/s), 
        faultplane (PlanarSurface) object, fault half width distance projected at 
        surface.
    '''
    # --------------------------------------------------------------------------
    # Get multigmpe from config
    # --------------------------------------------------------------------------
    mgmpe = MultiGMPE.__from_config__(gmpeconf)
    IMT = imt.PGA() # will be in units of "g"
    # --------------------------------------------------------------------------
    # Set rupture basics
    # --------------------------------------------------------------------------
    gm = {}
    for mag in np.arange(calcconf['magrange']['magmin'], calcconf['magrange']['magmax'] + \
            calcconf['magrange']['magstep']/2, calcconf['magrange']['magstep']):
        # --------------------------------------------------------------------------
        # Rupture context
        # --------------------------------------------------------------------------
        evconf['mag'] = mag
        rctx, faultplane, xcorr = createRuptureContext(evconf, calcconf)
        # --------------------------------------------------------------------------
        # Distance and Source contexts
        # --------------------------------------------------------------------------
        if 'grid' in calcconf and calcconf['grid']['compute']:
            dctx, sctx = make_pga_grid(evconf, calcconf, rctx, faultplane, xcorr, bPlots = calcconf['plots'])
        if 'points' in calcconf and calcconf['points']['compute']:
            dctx, sctx = make_pga_lop(evconf, calcconf, rctx, faultplane, xcorr, bPlots = calcconf['plots'])
        if 'pt' in calcconf and calcconf['pt']['compute']:
            dctx, sctx = make_pga_pt(evconf, calcconf, rctx, faultplane, xcorr, bPlots = calcconf['plots'])
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
        gm[mag] = [lng2cm(lmean_mgmpe), faultplane, xcorr]
    return gm


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
