from math import log10
from openquake.hazardlib.scalerel.wc1994 import WC1994

class WC1994_Ext(WC1994):
    """
    Wells, D.L., and K.J. Coppersmith, 1994. New empirical relationships among magnitude, rupture
    length, rupture width, rupture area, and surface displacement. Bulletin of the Seismological 
    Society of America, 84(4), pp 974-1002.

    Extends the base class with magnitude-length and length-magnitude scaling relationships,
    plus magnitude-width and width-magnitude scaling relationships.
    """

    def get_median_SRlength(self, mag, rake=None):
        """
        Calculates median fault length from magnitude.
        Note this uses the SRL (surface) relation, not RLD (subsurface)
        """
        if rake is None:
            # Their 'ALL' case
            return pow(10., (0.69 * mag) - 3.22)
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            return pow(10., (0.74 * mag) - 3.55)
        elif rake > 0.:
            # thrust/reverse
            return pow(10., (0.63 * mag) - 2.86)
        else:
            # normal
            return pow(10., (0.50 * mag) - 2.01)

    def get_median_length(self, mag, rake=None):
        """
        Calculates median fault length from magnitude.
        Note this uses the RLD relation, not SRL
        """
        if rake is None:
            # Their 'ALL' case
            return pow(10., (0.59 * mag) - 2.44)
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            return pow(10., (0.62 * mag) - 2.57)
        elif rake > 0.:
            # thrust/reverse
            return pow(10., (0.58 * mag) - 2.42)
        else:
            # normal
            return pow(10., (0.50 * mag) - 1.88)


    def get_median_mag_from_length(self, flen, rake=None):
        """
        Calculates median magnitude from fault length.
        """
        if rake is None:
            # Their 'ALL' case
            return 4.38 + (1.49 * log10(flen))
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            return 4.33 + (1.49 * log10(flen))
        elif rake > 0.:
            # thrust/reverse
            return 4.49 + (1.49 * log10(flen))
        else:
            # normal
            return 4.34 + (1.54 * log10(flen))


    def get_median_width(self, mag, rake=None):
        """
        Calculates median fault width (RW) from magnitude.
        """
        if rake is None:
            # Their 'ALL' case
            return pow(10., (0.32 * mag) - 1.01)
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            return pow(10., (0.27 * mag) - 0.76)
        elif rake > 0.:
            # thrust/reverse
            return pow(10., (0.41 * mag) - 1.61)
        else:
            # normal
            return pow(10., (0.35 * mag) - 1.14)


    def get_median_mag_from_width(self, fwid, rake=None):
        """
        Calculates median magnitude from fault width (RW).
        """
        if rake is None:
            # Return average of strike-slip and dip-slip curves
            return 4.06 + (2.25 * log10(fwid))
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            return 3.80 + (2.59 * log10(fwid))
        elif rake > 0.:
            # thrust/reverse
            return 4.37 + (1.95 * log10(fwid))
        else:
            # normal
            return 4.04 + (2.11 * log10(fwid))


