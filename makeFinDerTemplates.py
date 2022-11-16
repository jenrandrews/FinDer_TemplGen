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

    gm = woq.computeGM(gmpeconf, evconf, calcconf)

    dkm = calcconf['grid']['griddkm']
    for mag in gm:
        lmean_mgmpe, faultplane, xcorr = gm[mag]
        flen = faultplane.get_area()/faultplane.get_width()
        oname = 'template_L%.6f_Azi0.txt' % flen
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

