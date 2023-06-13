#!/usr/bin/env python

import os
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
    rdir = evconf['evmech']['geometry']
    rlist = sorted(os.listdir(rdir))
    for f in rlist:
        evconf['evmech']['geometry'] = os.path.join(rdir, f)
        nind = f.replace('rupture_','').replace('.json','')
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
        
        gm, evconf = woq.computeGM(gmpeconf, evconf, calcconf)
        mag = list(gm.keys())[0]
        lmean_mgmpe, faultplane, xcorr = gm[mag]
        maxpga = max(lmean_mgmpe)
        logging.info('Max PGA: %.4f' % maxpga)
        if maxpga < log10(2.):
            logging.info('All PGA below 2 cm/s/s threshold, no output files written')
            exit()
        oname = 'data_0_%s_%.4f_%.4f_%.2f_%03d' % (nind, evconf['evloc']['centroid_lat'], evconf['evloc']['centroid_lon'], mag, evconf['evmech']['strike'])
        fout = open(oname, 'w')
        fout.write('#  %.4f  %.4f  %.2f  %03d\n' % (evconf['evloc']['centroid_lat'], evconf['evloc']['centroid_lon'], mag, evconf['evmech']['strike']))
        for pga, lat, lon, stnn in zip(lmean_mgmpe, lats, lons, stnnames):
    #                    fout.write('%.5f %.5f %s %.5f\n' % (lat, lon, stnn, pga))
            fout.write('%.5f %.5f %.5f\n' % (lat, lon, pga))
        fout.close()
        # write data to file
        if calcconf['plots']:
            import matplotlib.pyplot as plt
            flat = []
            flon = []
            for x in evconf['flist']:
                flat.append(x.latitude)
                flon.append(x.longitude)
            plt.scatter(lons, lats, c=lmean_mgmpe)
            plt.plot(flon, flat)
            plt.scatter(evconf['evloc']['centroid_lon'], evconf['evloc']['centroid_lat'], marker='*', s=80)
            plt.colorbar()
            plt.savefig('%s.png' % oname)
            plt.close()

