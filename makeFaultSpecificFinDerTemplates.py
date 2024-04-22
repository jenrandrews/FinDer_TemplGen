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


def writeSetFile(calcconf, name, minflen, maxflen, minmag, maxmag, coords):
    poly = woq.createCentroidPolygon(calcconf, coords)
    npstr = len(poly.exterior.coords)
    pstr = ' '.join([f'{p[1]:.6f},{p[0]:.6f}' for p in poly.exterior.coords])
    if False:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter([p[0] for p in coords], [p[1] for p in coords], marker='x', c='r')
        plt.plot([p[0] for p in poly.exterior.coords], [p[1] for p in poly.exterior.coords], c='b')
        plt.savefig('buffer.png')
    fout = open(f"{name}.txt", 'w')
    fout.write(f"TEMPLATE_SET_NAME   {name}\n")
    fout.write("TEMPLATE_TYPE       SPECIFIC\n")
    fout.write(f"TEMPLATE_DIRECTORY  conf/Templates_PGA_{name}\n")
    fout.write(f"D_KM                {calcconf['grid']['griddkm']:.1f}\n")
    fout.write("LATITUDE_FILE       latitude.dat\n")
    fout.write("LONGITUDE_FILE      longitude.dat\n")
    fout.write(f"MIN_LENGTH          {minflen-1:.0f}\n")
    fout.write(f"MAX_LENGTH          {maxflen+1:.0f}\n")
    fout.write(f"MIN_MAG             {minmag-0.5:.1f}\n")
    fout.write(f"MAX_MAG             {maxmag:.1f}\n")
    fout.write("FAST_MAG_RANGE      1.0\n")
    fout.write("FAST_LATLON_RANGE   1.0\n")
    fout.write("TEMPLATE_INFO_FILE  template_info.txt\n")
    fout.write(f"CENTROID_POLYGON    {npstr} {pstr}\n")
    fout.close()
    return


if __name__ == "__main__":
    import os
    import logging.config
    from DEFLOG import DEFLOG
    DEFLOG['handlers']['fileHandler']['filename'] = \
            'makeFinDerTemplates_%s.log' % time.strftime('%y%m%dT%H%M%S', time.gmtime(time.time()))
    logging.config.dictConfig(DEFLOG)
    logger = logging.getLogger(__name__)

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

    gm, dummy, templ_sets = woq.computeGM(gmpeconf, evconf, calcconf)
    for tset in sorted(templ_sets):
        odir = '_'.join(['Templates', 'PGA', calcconf['fault-specific']['name'], str(tset)])
        if not os.path.isdir(odir):
            os.mkdir(odir)

        dkm = calcconf['grid']['griddkm']
        if 'rupinfo' in calcconf and calcconf['rupinfo']:
            fout = open(os.path.join(odir, 'rupinfo.tbl'), 'w')
            fout.write('Depth: {:.1f} Dip: {:.1f}\n'.format(evconf['evloc']['hypo_depth'], 
                evconf['evmech']['dip']))
            fout2 = open(os.path.join(odir, 'template_info.txt'), 'w')
            with open(os.path.join(odir, 'latitude.dat'), 'w') as fout3:
                for lat in sorted(set(templ_sets[tset]['mesh'].lats)):
                    fout3.write(f'{lat:.6f}\n')
            with open(os.path.join(odir, 'longitude.dat'), 'w') as fout4:
                for lon in sorted(set(templ_sets[tset]['mesh'].lons)):
                    fout4.write(f'{lon:.6f}\n')
        minflen = maxflen = None
        for mag in sorted([m for m in templ_sets[tset] if m not in ['mesh', 'mask']]):
            for (centroid_lat, centroid_lon) in sorted(templ_sets[tset][mag]):
                lmean_mgmpe, faultplane, rjb = gm[mag][(centroid_lat, centroid_lon)]
                fwid = faultplane.get_width()
                flen = faultplane.get_area()/fwid
                if minflen is None or flen < minflen:
                    minflen = flen
                if maxflen is None or flen > maxflen:
                    maxflen = flen
                oname = f'template_L{flen:.6f}_W{fwid:.4f}_{centroid_lat:.4f}_{centroid_lon:.4f}.txt'
                hstr = '{:d} {:d}\n{:f} {:d} {:.1f}\n'.format(lmean_mgmpe.shape[1], lmean_mgmpe.shape[0], flen, 0, dkm)
                woq.np.savetxt(os.path.join(odir, oname), lmean_mgmpe, fmt='%.6e', header=hstr)
                formatHeader(os.path.join(odir, oname))

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
                with open(os.path.join(odir, onamer), 'w') as foutr:
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
            coords = [c for m in templ_sets[tset] for c in templ_sets[tset][m] if m not in ['mesh', 'mask']]
            mags = [m for m in templ_sets[tset] if m not in ['mesh', 'mask']]
            name = f'{calcconf["fault-specific"]["name"]}_{tset:d}'
            writeSetFile(calcconf, name, minflen, maxflen, min(mags), max(mags), [(c[1], c[0]) for c in coords])

