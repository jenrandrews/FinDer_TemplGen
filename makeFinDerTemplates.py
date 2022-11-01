#!/usr/bin/env python

#############
# N.B. This script is based heavily on shakemap_src/tests/shakelib/multigmpe_test.py
#############

# stdlib imports
import copy
import os
import sys
import logging
import time as time
import argparse
from math import log10, sqrt, exp, floor, sin, cos, tan, radians, ceil
import numpy as np

# openquake imports
from openquake.hazardlib import imt, const
from openquake.hazardlib.gsim.base import RuptureContext
from openquake.hazardlib.gsim.base import DistancesContext
from openquake.hazardlib.gsim.base import SitesContext
from openquake.hazardlib.geo import Point, PlanarSurface, RectangularMesh

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


def lng2cm(inx):
    inx = np.log10(np.exp(inx)) + log10(980.665)
    return inx


def formatHeader(oname):
    fin = open(oname, 'r')
    text = fin.read()
    text = text.replace('# ','')
    text = text.replace('\n\n','\n')
    fin.close()
    fout = open(oname, 'w')
    fout.write(text)
    fout.close()
    return


def getScalingRelation(conf):
    '''
    Allowing for extension to other scaling relations. However, all will require an Ext class to 
    provide length and width rather than area.
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


def make_pga_grid(gmpeconf, evconf, calcconf):
    # --------------------------------------------------------------------------
    # Get multigmpe from config
    # --------------------------------------------------------------------------
    mgmpe = MultiGMPE.__from_config__(gmpeconf)
    IMT = imt.PGA() # will be in units of "g"

    # --------------------------------------------------------------------------
    # Set rupture basics
    # --------------------------------------------------------------------------
    scalrel = getScalingRelation(calcconf)
    strike = evconf['evmech']['strike']
    seismogenic_depth = evconf['seisstruc']['seismogenic_depth']
    dkm = calcconf['griddkm']
    chosen_vs30 = calcconf['vs30']
    for mag in np.arange(calcconf['magrange']['magmin'], calcconf['magrange']['magmax'] + \
            calcconf['magrange']['magstep']/2, calcconf['magrange']['magstep']):
        # ----
        # Empirical estimate for maximum distance needed in template
        # ----
        maxdist = max([mag*115. - 400., mag*40 - 95])
        # --------------------------------------------------------------------------
        # Rupture context
        # --------------------------------------------------------------------------
        rctx = RuptureContext()
        rctx.rake = evconf['evmech']['rake']
        rctx.dip = evconf['evmech']['dip']
        rctx.mag = mag
        rctx.hypo_depth = evconf['evloc']['hypo_depth']
        centroid_depth = evconf['evloc']['hypo_depth']
        diprad = radians(rctx.dip)
        # ----
        # If no depth is given, use Kaklamanos el al. (2010) to set
        # ----
        if centroid_depth is None:
            if rctx.rake is None:
                centroid_depth = 7.08 + (0.61 * mag) # Kaklamanos eq.4 general
            if (-45 <= rctx.rake <= 45) or (rctx.rake >= 135) or (rctx.rake <= -135):
                centroid_depth = 5.63 + (0.68 * mag) # Kaklamanos eq.4 SS
            else:
                centroid_depth = 11.24 - (0.2 * mag) # Kaklamanos eq.4 non-SS
        flen = scalrel.get_median_length(mag, rctx.rake)
        fwid = scalrel.get_median_width(mag, rctx.rake)
        # Note that mag range extends width relation beyond validity, so fix assumes
        # len = wid when wid > len
        if fwid > flen: 
            rctx.width = flen
        else:
            rctx.width = fwid
        # For FinDer, centroid must be kept at center of template
        # Use initial hypo depth as centroid depth unless width requires us to move it
        # Hypo depth empirically found to be at about 60% depth, but simpler to treat
        # as 50% here
        if rctx.width > seismogenic_depth:
            rctx.width = seismogenic_depth
            rctx.ztor = 0.1
            maxz = rctx.width
        elif centroid_depth < sin(diprad) * rctx.width * 0.5:
            # Note centroid_depth is now altered!
            rctx.ztor = 0.1
            maxz = sin(diprad) * rctx.width
            centroid_depth = maxz/2.
            logging.info('N.B. new centroid depth: %.2f' % (maxz/2.)) 
        else:
            offset = (sin(diprad) * rctx.width * 0.5)
            rctx.ztor = centroid_depth - offset
            maxz = centroid_depth + offset
        logging.info(r'Mag %.1f  Max dist %.4f Rup len %.4f Rup wid %.4f/%.4f ztor %.2f zmax %.2f' % \
                (mag, maxdist, flen, fwid, rctx.width, rctx.ztor, maxz))
        
        # --------------------------------------------------------------------------
        # Distance context
        # --------------------------------------------------------------------------
        dd2km = 111.194926645
        km2dd = 0.008993216
        xcorr = (cos(diprad) * 0.5 * rctx.width) * km2dd
        ylen = flen * km2dd
        # Fault plane defined as lat: 0 to fault_len; lon: -0.5 projected width to +0.5 projected width
        faultplane = PlanarSurface(
                0., rctx.dip, 
                top_left=Point(depth=rctx.ztor, latitude=0., longitude=-1.*xcorr),
                top_right=Point(depth=rctx.ztor, latitude=ylen, longitude=-1*xcorr),
                bottom_left=Point(depth=maxz, latitude=0., longitude=xcorr),
                bottom_right=Point(depth=maxz, latitude=ylen, longitude=xcorr))
        delta = dkm * km2dd
        min_lat = -1. * delta * ceil(maxdist / dkm)
        max_lat = delta * ceil((maxdist + flen)/ dkm)
        inlats = np.arange(min_lat, max_lat+(delta/2.), delta)
        min_lon = -1. * delta * ceil((maxdist + xcorr) / dkm)
        max_lon = delta * ceil((maxdist + xcorr) / dkm)
        inlons = np.arange(min_lon, max_lon+(delta/2.), delta)
        # Mask for asymmetric templates
        masklonind = min([i for i,x in enumerate(inlons) if x > -1.*xcorr])
        lons, lats = np.meshgrid(inlons, inlats)
        mesh = RectangularMesh(lons=lons, lats=lats, depths=None)

        dctx = DistancesContext()
        dctx.rjb = faultplane.get_joyner_boore_distance(mesh)
        dctx.rjb_var = None
        dctx.rrup_var = None
        dctx.rvolc = np.zeros_like(dctx.rjb) # no correction for travel path in volcanic region
        rhyp = []
        for y in inlats:
            for x in inlons:
                rhyp.append(sqrt(pow(x*dd2km,2.) + pow((y-(ylen/2.))*dd2km,2.) + pow(rctx.hypo_depth,2.)))
        dctx.rhypo = np.asarray(rhyp).reshape(len(inlats), len(inlons))
        dctx.rx = faultplane.get_rx_distance(mesh)
        dctx.ry0 = faultplane.get_ry0_distance(mesh)
        # For long faults, the get_min_distance function appears to only use the corners to 
        # compute horizontal distance, resulting in Rrup >> actual
        # If future libs fix this issue, reinstate following line:
        #dctx.rrup = faultplane.get_min_distance(mesh).reshape(mesh.shape)
        if rctx.dip == 90.:
            dctx.rrup = np.sqrt(np.square(dctx.rjb) + pow(rctx.ztor,2))
        else:
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

        if calcconf['plots']:
            import matplotlib.pyplot as plt
            plt.imshow(dctx.rrup)
            plt.colorbar()
            plt.savefig('rrup_M%.1f.png' % mag)
            plt.close()
            plt.imshow(dctx.rjb)
            plt.colorbar()
            plt.savefig('rjb_M%.1f.png' % mag)
            plt.close()
            plt.imshow(dctx.rx)
            plt.colorbar()
            plt.savefig('rx_M%.1f.png' % mag)
            plt.close()
            plt.imshow(dctx.ry0)
            plt.colorbar()
            plt.savefig('ry0_M%.1f.png' % mag)
            plt.close()
            plt.imshow(dctx.rhypo)
            plt.colorbar()
            plt.savefig('rhypo_M%.1f.png' % mag)
            plt.close()

        # --------------------------------------------------------------------------
        # Site context
        # --------------------------------------------------------------------------
        sctx = SitesContext()
        sctx.sids = np.arange(len(inlons)*len(inlats)).reshape(dctx.rjb.shape)
        #sctx.sids = np.arange(len(dctx.rjb))
        sctx.vs30 = np.ones_like(dctx.rjb) * chosen_vs30
        sctx.vs30measured = np.full_like(dctx.rjb, False, dtype="bool")
        sctx.z1pt0_ask14_cal = sites.Sites._z1pt0_from_vs30_ask14_cal(sctx.vs30)
        sctx.z1pt0_cy14_cal = sites.Sites._z1pt0_from_vs30_cy14_cal(sctx.vs30)
        sctx.z2pt5_cb14_cal = sites.Sites._z2pt5_from_vs30_cb14_cal(sctx.vs30) / 1000.0
        sctx.backarc = np.zeros_like(dctx.rjb, dtype=np.int16) # forearc/backarc is unknown

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
            logging.info('Converting component')
            bk17 = BooreKishida2017(
                mgmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT, const.IMC.GREATER_OF_TWO_HORIZONTAL
            )
            lmean_mgmpe = bk17.convertAmpsOnce(IMT, lmean_mgmpe, dctx.rrup, rctx.mag)

        # Method get_mean_and_stddevs() of actual GMPE implementations is supposed to return the 
        # mean value as a natural logarithm of intensity.
        lmean_mgmpe = lng2cm(lmean_mgmpe)
     
        oname = 'template_L%.6f_Azi0.txt' % flen
        hstr = '%d %d\n%f %d %.1f\n' % (lmean_mgmpe.shape[1], lmean_mgmpe.shape[0], flen, 0, dkm)
        np.savetxt(oname, lmean_mgmpe, fmt='%.6e', header=hstr)
        formatHeader(oname)

        if calcconf['plots']:
            import matplotlib.pyplot as plt
            plt.imshow(lmean_mgmpe)
            plt.colorbar()
            plt.savefig('templ_M%.1f.png' % mag)
            plt.close()

        if calcconf['asym']:
            oname = 'template_L%.6f_Azi0_asym.txt' % flen
            hstr = '%d %d\n%f %d %.1f\n' % (lmean_mgmpe.shape[1], lmean_mgmpe.shape[0], flen, 0, dkm)
            lmean_mgmpe[:,:masklonind] = -2.0
            np.savetxt(oname, lmean_mgmpe, fmt='%.6e', header=hstr)
            formatHeader(oname)

            if calcconf['plots']:
                import matplotlib.pyplot as plt
                plt.imshow(lmean_mgmpe)
                plt.colorbar()
                plt.savefig('templ_M%.1f_asym.png' % mag)
                plt.close()
    return


if __name__ == "__main__":
    # --------------------------------------------------------------------------
    # Minimal config dictionary
    # --------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Make FinDer templates using ShakeMap/openquake')
    parser.add_argument('--gmpe', '-g', action='store', required=True, 
        dest='gmpeconf', help='GMPE configuration file')
    parser.add_argument('--event', '-e', action='store', required=True, 
        dest='evconf', help='Event configuration file')
    parser.add_argument('--calc', '-c', action='store', required=True, 
        dest='calcconf', help='Calculation configuration file')
    args = parser.parse_args()

    gmpeconf = importConfig(args.gmpeconf)
    evconf = importConfig(args.evconf)
    calcconf = importConfig(args.calcconf)

    logging.basicConfig(filename='makeFinDerTemplates_%s.log'% \
      time.strftime('%y%m%dT%H%M%S', time.gmtime(time.time())), level=logging.INFO)

    make_pga_grid(gmpeconf, evconf, calcconf)


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
