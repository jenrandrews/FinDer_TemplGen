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

    if 'grid' in calcconf and calcconf['grid']['compute']:
        if 'rupinfo' in calcconf and calcconf['rupinfo']:
            fout = open('rupinfo.tbl', 'w')
            fout.write('Depth: {:.1f} Dip: {:.1f}\n'.format(evconf['evloc']['hypo_depth'], 
                evconf['evmech']['dip']))
            fout2 = open('template_info.txt', 'w')

    for f in rlist:
        if f.find('.json') == -1:
            continue
        evconf['evmech']['geometry'] = os.path.join(rdir, f)
        nind = f.replace('rupture_','').replace('.json','')
        gm, evconf = woq.computeGM(gmpeconf, evconf, calcconf)
        mag = list(gm.keys())[0]
        lmean_mgmpe, faultplane, xcorr = gm[mag]
        maxpga = np.amax(lmean_mgmpe)
        logging.info('Max PGA: %.4f' % maxpga)
        if maxpga < log10(2.):
            logging.info('All PGA below 2 cm/s/s threshold, no output files written')
            exit()

        # --------------------------------------------------------------------------
        # Write out data
        # --------------------------------------------------------------------------
        if 'points' in calcconf and calcconf['points']['compute']:
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
        
            oname = 'data_0_{}_{:.4f}_{:.4f}_{:.2f}_{:03d}'.format(nind, 
                    evconf['evloc']['centroid_lat'], 
                    evconf['evloc']['centroid_lon'], 
                    mag, 
                    evconf['evmech']['strike'])
            fout = open(oname, 'w')
            fout.write('#  {:.4f}  {:.4f}  {:.2f}  {:03d}\n'.format(evconf['evloc']['centroid_lat'], 
                    evconf['evloc']['centroid_lon'], 
                    mag, 
                    evconf['evmech']['strike']))
            for pga, lat, lon, stnn in zip(lmean_mgmpe, lats, lons, stnnames):
                #fout.write('%.5f %.5f %s %.5f\n' % (lat, lon, stnn, pga))
                fout.write('{:.5f} {:.5f} {:.5f}\n'.format(lat, lon, pga))
            fout.close()
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

        if 'grid' in calcconf and calcconf['grid']['compute']:
            dkm = calcconf['grid']['griddkm']
            flen = faultplane.get_area()/faultplane.get_width()
            oname = 'template_L%.6f_Azi0.txt' % flen
            if 'rupinfo' in calcconf and calcconf['rupinfo']:
                fout.write('{:.1f} {:.4f} {:.4f} {:.2f} {:.2f}\n'.format(
                    mag, 
                    flen, 
                    faultplane.get_width(), 
                    min(faultplane.top_left.depth, faultplane.top_right.depth), 
                    max(faultplane.bottom_left.depth, faultplane.bottom_right.depth)))
                fout2.write('{:.6f} {} {:.1f}\n'.format(flen, oname, mag))
            hstr = '%d %d\n%f %d %.1f\n' % (lmean_mgmpe.shape[1], lmean_mgmpe.shape[0], flen, 0, dkm)
            woq.np.savetxt(oname, lmean_mgmpe, fmt='%.6e', header=hstr)
            formatHeader(oname)

            if calcconf['plots']:
                import matplotlib.pyplot as plt
                plt.imshow(lmean_mgmpe, origin='lower')
                plt.colorbar()
                plt.savefig('templ_M%.1f.png' % mag)
                plt.close()

            if calcconf['grid']['asym']:
                oname = 'template_L%.6f_Azi0_asym.txt' % flen
                hstr = '%d %d\n%f %d %.1f\n' % (lmean_mgmpe.shape[1], lmean_mgmpe.shape[0], flen, 0, dkm)
                masklonind = woq.floor(lmean_mgmpe.shape[1]/2) - round(xcorr/dkm)
                lmean_mgmpe[:,:masklonind] = -2.0
                woq.np.savetxt(oname, lmean_mgmpe, fmt='%.6e', header=hstr)
                formatHeader(oname)

                if calcconf['plots']:
                    import matplotlib.pyplot as plt
                    plt.imshow(lmean_mgmpe, origin='lower')
                    plt.colorbar()
                    plt.savefig('templ_M%.1f_asym.png' % mag)
                    plt.close()

    if 'grid' in calcconf and calcconf['grid']['compute']:
        if 'rupinfo' in calcconf and calcconf['rupinfo']:
            fout.close()
            fout2.close()

