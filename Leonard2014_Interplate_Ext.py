from math import log10
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_Interplate

class Leonard2014_Interplate_Ext(Leonard2014_Interplate):
    """
    Leonard, M., 2014. Self-consistent earthquake fault-scaling relations:
    Update and extension to stable continental strike-slip faults.
    Bulletin of the Seismological Society of America, 104(6), pp 2953-2965.

    Extends the base class with magnitude-length and length-magnitude scaling relationships,
    plus magnitude-width and width-magnitude scaling relationships.
    """

    def get_median_length(self, mag, rake=None):
        """
        Calculates median fault length from magnitude.
        """
        if rake is None:
            # Return average of strike-slip and dip-slip curves
            lenSS = self.get_median_length(mag, 0.)
            lenDS = self.get_median_length(mag, 60.)
            return (lenSS + lenDS) / 2.
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # Strike slip
            if mag > 6.9:
                return pow(10.0, mag - 5.27)
            else:
                return pow(10.0, (mag - 4.17) / 1.667)
        else:
            # Dip slip (thrust or normal)
            if mag > 5.5:
                return pow(10.0, (mag - 4.24) / 1.667)
            else:
                return pow(10.0, (mag - 4.) / 2.)


    def get_median_mag_from_length(self, flen, rake=None):
        """
        Calculates median magnitude from fault length.
        """
        if rake is None:
            # Return average of strike-slip and dip-slip curves
            magSS = self.get_median_mag_from_length(flen, 0.)
            magDS = self.get_median_mag_from_length(flen, 60.)
            return (magSS + magDS) / 2.
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # Strike slip
            if flen > 45.:
                return 5.27 + log10(flen)
            else:
                return 4.17 + (1.667 * log10(flen))
        else:
            # Dip slip (thrust or normal)
            if flen > 5.4:
                return 4.24 + (1.667 * log10(flen))
            else:
                return 4. + (2. * log10(flen))


    def get_median_width(self, mag, rake=None):
        """
        Calculates median fault width from magnitude.
        """
        if rake is None:
            # Return average of strike-slip and dip-slip curves
            widSS = self.get_median_width(mag, 0.)
            widDS = self.get_median_width(mag, 60.)
            return (widSS + widDS) / 2.
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # Strike slip
            return pow(10.0, (mag - 3.88) / 2.5)
        else:
            # Dip slip (thrust or normal)
            return pow(10.0, (mag - 3.63) / 2.5)


    def get_median_mag_from_width(self, fwid, rake=None):
        """
        Calculates median magnitude from fault width.
        """
        if rake is None:
            # Return average of strike-slip and dip-slip curves
            magSS = self.get_median_mag_from_width(fwid, 0.)
            magDS = self.get_median_mag_from_width(fwid, 60.)
            return (magSS + magDS) / 2.
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # Strike slip
            return 3.88 + (2.5 * log10(fwid))
        else:
            # Dip slip (thrust or normal)
            return 3.63 + (2.5 * log10(fwid))


