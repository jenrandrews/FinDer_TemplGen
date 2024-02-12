#!/usr/bin/env python

import time as time
import argparse
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

    gm, dummy, dummy = woq.computeGM(gmpeconf, evconf, calcconf)

    dkm = calcconf['grid']['griddkm']
    if 'rupinfo' in calcconf and calcconf['rupinfo']:
        fout = open('rupinfo.tbl', 'w')
        fout.write('Depth: {:.1f} Dip: {:.1f}\n'.format(evconf['evloc']['hypo_depth'], 
            evconf['evmech']['dip']))
        fout2 = open('template_info.txt', 'w')
    for mag in gm:
        for (centroid_lat, centroid_lon) in gm[mag]:
            lmean_mgmpe, faultplane, rjb = gm[mag][(centroid_lat, centroid_lon)]
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
                plt.imshow(lmean_mgmpe)
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
                    plt.imshow(lmean_mgmpe)
                    plt.colorbar()
                    plt.savefig('templ_M%.1f_asym.png' % mag)
                    plt.close()
    if 'rupinfo' in calcconf and calcconf['rupinfo']:
        fout.close()
        fout2.close()

