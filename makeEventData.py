#!/usr/bin/env python

import os
import sys
import time as time
import argparse
import numpy as np
from math import log10

import wrapOQ as woq
from openquake.hazardlib.geo.surface import ComplexFaultSurface


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
    if os.path.isdir(rdir):
        rlist = [os.path.join(rdir, r) for r in sorted(os.listdir(rdir))]
    elif os.path.isfile(rdir):
        rlist = [rdir]
    else:
        logging.error(f'Rupture geometry {rdir} is incorrectly specified')
        exit()

    if 'grid' in calcconf and calcconf['grid']['compute']:
        if 'rupinfo' in calcconf and calcconf['rupinfo']:
            fout = open('rupinfo.tbl', 'w')
            fout.write('Depth: {:.1f} Dip: {:.1f}\n'.format(evconf['evloc']['hypo_depth'], 
                evconf['evmech']['dip']))
            fout2 = open('template_info.txt', 'w')

    for f in sorted(rlist):
        if f.find('.json') == -1:
            continue
        evconf['evmech']['geometry'] = f
        nind = f.replace('rupture_','').replace('.json','')
        gm, evconf, dummy = woq.computeGM(gmpeconf, evconf, calcconf)
        for mag in gm:
            for (centroid_lat, centroid_lon) in gm[mag]:
#                lmean_mgmpe, faultplane = gm[mag][(centroid_lat, centroid_lon)]
                lmean_mgmpe, faultplane, rjb = gm[mag][(centroid_lat, centroid_lon)]
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
                    round(evconf['evmech']['strike']))
            fout3 = open(oname, 'w')
            fout3.write('#  {:.4f}  {:.4f}  {:.2f}  {:03d}\n'.format(evconf['evloc']['centroid_lat'], 
                    evconf['evloc']['centroid_lon'], 
                    mag, 
                    round(evconf['evmech']['strike'])))
            for pga, lat, lon, stnn, dist in zip(lmean_mgmpe, lats, lons, stnnames, rjb):
                #fout.write('%.5f %.5f %s %.5f\n' % (lat, lon, stnn, pga))
                #fout.write('{:.5f} {:.5f} {:.5f}\n'.format(lat, lon, pga))
                fout3.write('{:.5f} {:.5f} {:.5f} {:.2f}\n'.format(lat, lon, pga, dist))
            fout3.close()
            if calcconf['plots']:
                import matplotlib.pyplot as plt
                flat = []
                flon = []
                if isinstance(faultplane, ComplexFaultSurface):
                    top = faultplane.surface_nodes[0].nodes[0].nodes[0].nodes[0]
                    bottom = faultplane.surface_nodes[0].nodes[-1].nodes[0].nodes[0]
                    flon.extend([float(x) for x in top.to_str().split('[')[1].split(']')[0].split(',')[::3]])
                    flat.extend([float(x) for x in top.to_str().split('[')[1].split(']')[0].split(',')[1::3]])
                    flon.extend([float(x) for x in bottom.to_str().split('[')[1].split(']')[0].split(',')[::3]])
                    flat.extend([float(x) for x in bottom.to_str().split('[')[1].split(']')[0].split(',')[1::3]])
                else:
                    for x in [faultplane.top_left, faultplane.top_right, faultplane.bottom_right, faultplane.bottom_left, faultplane.top_left]:
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
                if isinstance(faultplane, ComplexFaultSurface):
                    top = faultplane.get_top_edge_depth()
                    bottom = float(faultplane.surface_nodes[0].nodes[-1].nodes[0].nodes[0].to_str().split('[')[1].split(']')[0].split(',')[-1])
                else:
                    top = min(faultplane.top_left.depth, faultplane.top_right.depth)
                    bottom = max(faultplane.bottom_left.depth, faultplane.bottom_right.depth)
                fout.write('{:.1f} {:.4f} {:.4f} {:.2f} {:.2f}\n'.format(
                    mag, 
                    flen, 
                    faultplane.get_width(), 
                    top,
                    bottom))
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

