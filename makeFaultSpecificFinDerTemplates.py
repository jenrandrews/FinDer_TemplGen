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

    gm, dummy = woq.computeGM(gmpeconf, evconf, calcconf)

#    for mag in gm:
#        for centroid_lat, centroid_lon in gm[mag]:
#            print(f'{mag:.1f}, {centroid_lat:.6f}, {centroid_lon:.6f}')
#    exit()

    dkm = calcconf['grid']['griddkm']
    if 'rupinfo' in calcconf and calcconf['rupinfo']:
        fout = open('rupinfo.tbl', 'w')
        fout.write('Depth: {:.1f} Dip: {:.1f}\n'.format(evconf['evloc']['hypo_depth'], 
            evconf['evmech']['dip']))
        fout2 = open('template_info.txt', 'w')
    for mag in gm:
        for centroid_lat, centroid_lon in gm[mag]:
            lmean_mgmpe, faultplane = gm[mag][(centroid_lat, centroid_lon)]
            fwid = faultplane.get_width()
            flen = faultplane.get_area()/fwid
            oname = f'template_L{flen:.6f}_W{fwid:.4f}_{centroid_lat:.4f}_{centroid_lon:.4f}.txt'
            hstr = '{:d} {:d}\n{:f} {:d} {:.1f}\n'.format(lmean_mgmpe.shape[1], lmean_mgmpe.shape[0], flen, 0, dkm)
            woq.np.savetxt(oname, lmean_mgmpe, fmt='%.6e', header=hstr)
            formatHeader(oname)

            if calcconf['plots']:
                import matplotlib.pyplot as plt
                plt.imshow(lmean_mgmpe, origin='lower')
                plt.colorbar()
                plt.savefig(f'templ_M{mag:.1f}_{centroid_lat:.4f}_{centroid_lon:.4f}.png')
                plt.close()

            ## Write out rupture_* files
            onamer = f'rupture_L{flen:.6f}_W{fwid:.4f}_{centroid_lat:.4f}_{centroid_lon:.4f}.txt'
            top = faultplane.surface_nodes[0].nodes[0].nodes[0].nodes[0]
            bottom = faultplane.surface_nodes[0].nodes[-1].nodes[0].nodes[0]
            with open(onamer, 'w') as foutr:
                for i in range(len(top.text)//3):
                    foutr.write(f'{top.text[i*3+1]} {top.text[i*3]} {top.text[i*3+2]}\n')
                for i in range(len(bottom.text)//3-1, -1, -1):
                    foutr.write(f'{bottom.text[i*3+1]} {bottom.text[i*3]} {bottom.text[i*3+2]}\n')
                foutr.write(f'{top.text[1]} {top.text[0]} {top.text[2]}\n')

            if 'rupinfo' in calcconf and calcconf['rupinfo']:
                fout.write('{:.1f} {:.4f} {:.4f} {:.2f} {:.2f}\n'.format(
                    mag, 
                    flen, 
                    faultplane.get_width(), 
                    min([top.text[i*3+2] for i in range(len(top.text)//3)]),
                    min([bottom.text[i*3+2] for i in range(len(bottom.text)//3)])))
                fout2.write('{:.6f} {:.6f} {:.6f} {:.6f} {:.1f} {} {}\n'.format(
                    flen, fwid, centroid_lat, centroid_lon, mag, oname, onamer
                    ))

    if 'rupinfo' in calcconf and calcconf['rupinfo']:
        fout.close()
        fout2.close()

