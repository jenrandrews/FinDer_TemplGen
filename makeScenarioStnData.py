#!/usr/bin/env python

import sys
sys.path.append('/home/jena/RCET-FinDer/scripts')
import time as time
import argparse
import numpy as np
from math import log10

import calcdist
import wrapOQ as woq


if __name__ == "__main__":
    import logging.config
    from DEFLOG import DEFLOG
    DEFLOG['handlers']['fileHandler']['filename'] = \
            'makeFinDerTemplates_%s.log' % time.strftime('%y%m%dT%H%M%S', time.gmtime(time.time()))
    logging.config.dictConfig(DEFLOG)

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

    gmpeconf = woq.importConfig(args.gmpeconf)
    evconf = woq.importConfig(args.evconf)
    calcconf = woq.importConfig(args.calcconf)

    inStrikes = np.arange(0., 180., 15.).tolist()
    outStrikes = np.arange(0., 50., 15.).tolist()
    outStrikes.extend(np.arange(135., 230., 15.).tolist())
    outStrikes.extend(np.arange(315., 360., 15.).tolist())

    # Grid for event epicenters
    NZ_BB = [[-44.9, 165.4], [-38.0, 173.4], [-33.4, 171.0], [-33.4, 174.0], [-37.0, 179.6], 
            [-39.5, 179.6], [-48.1, 169.4], [-48.1, 165.4], [-44.9, 165.4]]
    NZbounds = [-48., -33.0, 165.5, 180.0]
    dDeg = 0.5
    glats = []
    glons = []
    polys = {}
    polys['lat'] = [x[0] for x in NZ_BB]
    polys['lon'] = [x[1] for x in NZ_BB]
    for lat in np.arange(NZbounds[0], NZbounds[1], dDeg):
        for lon in np.arange(NZbounds[2], NZbounds[3], dDeg):
            if calcdist.inRegion(polys, lon, lat):
                glats.append(lat)
                glons.append(lon)
    inNZ, outNZ, boolNZ = calcdist.inNZ(glats, glons)
   
    ### JADEBUG
    #import matplotlib.pyplot as plt
    #plt.scatter([x[0] for x in outNZ], [x[1] for x in outNZ])
    #plt.savefig('offshore.png')
    #exit()
    ### JADEBUG

    fin = open(calcconf['points']['points_file'], 'r')
    lats = []
    lons = []
    stnnames = []
    for l in fin:
        fs = l.split()
        lons.append(float(fs[1]))
        lats.append(float(fs[0]))
        stnnames.append(fs[-1])
    fin.close() 

    # Loop over epicenters
    ### symmetric: all strikes and inNZ locations
    ### asymmetric: restrict to relevant strikes and outNZ locations
    for locs, strikes in zip([inNz, outNZ], [inStrikes, outStrikes]):
        for loc in locs: 
            epilat = loc[1]
            epilon = loc[0]
            evconf['evloc']['centroid_lat'] = epilat
            evconf['evloc']['centroid_lon'] = epilon
            logging.info('Epicenter: %.6f %.6f' % (epilat, epilon))
            # Loop over strikes
            for strike in strikes:
                evconf['evmech']['strike'] = strike
                gm, dummy = woq.computeGM(gmpeconf, evconf, calcconf)
                logging.info('Strike: %d' % strike)
                for mag in gm:
                    for centroid_lat, centroid_lon in gm[mag]:
                        logging.info('Mag: %.2f' % mag)
                        lmean_mgmpe, faultplane, xcorr = gm[mag]
                        maxpga = max(lmean_mgmpe)
                        logging.info('Max PGA: %.4f' % maxpga)
                        if maxpga < log10(2.):
                            continue
                        oname = 'data_0_%.2f_%.2f_%.2f_%03d' % (epilat, epilon, mag, strike)
                        fout = open(oname, 'w')
                        fout.write('#  %.4f  %.4f  %.2f  %03d\n' % (epilat, epilon, mag, strike))
                        for pga, lat, lon, stnn in zip(lmean_mgmpe, lats, lons, stnnames):
        #                    fout.write('%.5f %.5f %s %.5f\n' % (lat, lon, stnn, pga))
                            fout.write('%.5f %.5f %.5f\n' % (lat, lon, pga))
                        fout.close()
                        # write data to file
                        if calcconf['plots']:
                            import matplotlib.pyplot as plt
                            flat = []
                            flon = []
                            for x in [faultplane.top_left, faultplane.top_right, faultplane.bottom_right, 
                                    faultplane.bottom_left, faultplane.top_left]:
                                flat.append(x.latitude)
                                flon.append(x.longitude)
                            plt.scatter(lons, lats, c=lmean_mgmpe)
                            plt.plot(flon, flat)
                            plt.scatter(epilon, epilat, marker='*', s=80)
                            plt.colorbar()
                            plt.savefig('%s.png' % oname)
                            plt.close()

