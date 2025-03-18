from math import log10

class Blaser2010_Interface():
    """
    Blaser, L., F. Kruger, M. Ohrnberger, F. Scherbaum 2010. Scaling relations of earthquake
    source parameter estimates with special focus on subduction environment.
    Bulletin of the Seismological Society of America, 100(6), pp 2914-2926.

    Implements magnitude-length and length-magnitude scaling relationships,
    plus magnitude-width and width-magnitude scaling relationships.
    """

    def get_median_length(self, mag, rake=None, setting=None):
        """
        Calculates median fault length from magnitude.
        """
        if rake is None:
            # Their 'ALL' case
            if setting == 'oceanic':
                return pow(10., (0.54 * mag) - 2.07)
            if setting == 'continental':
                return pow(10., (0.57 * mag) - 2.26)
            else:
                return pow(10., (0.57 * mag) - 2.31)
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            if setting == 'oceanic':
                return pow(10., (0.62 * mag) - 2.56)
            if setting == 'continental':
                return pow(10., (0.62 * mag) - 2.55)
            else:
                return pow(10., (0.64 * mag) - 2.69)
        elif rake > 0.:
            # thrust/reverse
            if setting == 'oceanic':
                return pow(10., (0.62 * mag) - 2.81)
            if setting == 'continental':
                return pow(10., (0.54 * mag) - 2.17)
            else:
                return pow(10., (0.57 * mag) - 2.37)
        else:
            # normal
            return pow(10., (0.52 * mag) - 1.91)
        return None

    def get_median_mag_from_length(self, flen, rake=None, setting=None):
        """
        Calculates median magnitude from fault length.
        """
        if rake is None:
            # Their 'ALL' case
            if setting == 'oceanic':
                return (log10(flen) + 2.07)/0.54
            if setting == 'continental':
                return (log10(flen) + 2.26)/0.57
            else:
                return (log10(flen) + 2.31)/0.57
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            if setting == 'oceanic':
                return (log10(flen) + 2.56)/0.62 
            if setting == 'continental':
                return (log10(flen) + 2.55)/0.62
            else:
                return (log10(flen) + 2.69)/0.64
        elif rake > 0.:
            # thrust/reverse
            if setting == 'oceanic':
                return (log10(flen) + 2.81)/0.62
            if setting == 'continental':
                return (log10(flen) + 2.17)/0.54
            else:
                return (log10(flen) + 2.37)/0.57
        else:
            # normal
            return (log10(flen) + 1.91)/0.52
        return None

    def get_median_width(self, mag, rake=None, setting=None):
        """
        Calculates median fault width from magnitude.
        """
        if rake is None:
            # Their 'ALL' case
            if setting == 'oceanic':
                return pow(10., (0.44 * mag) - 1.76)
            if setting == 'continental':
                return pow(10., (0.34 * mag) - 1.14)
            else:
                return pow(10., (0.41 * mag) - 1.56)
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            if setting == 'oceanic':
                return pow(10., (0.27 * mag) - 0.66)
            if setting == 'continental':
                return pow(10., (0.27 * mag) - 0.75)
            else:
                return pow(10., (0.33 * mag) - 1.12)
        elif rake > 0.:
            # thrust/reverse
            if setting == 'oceanic':
                return pow(10., (0.45 * mag) - 1.79)
            if setting == 'continental':
                return pow(10., (0.45 * mag) - 1.83)
            else:
                return pow(10., (0.46 * mag) - 1.86)
        else:
            # normal
            return pow(10., (0.36 * mag) - 1.20)
        return None

    def get_median_mag_from_width(self, fwid, rake=None, setting=None):
        """
        Calculates median magnitude from fault width.
        """
        if rake is None:
            # Their 'ALL' case
            if setting == 'oceanic':
                return (log10(fwid) + 1.76)/0.44
            if setting == 'continental':
                return (log10(fwid) + 1.14)/0.34
            else:
                return (log10(fwid) + 1.56)/0.41
        elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            if setting == 'oceanic':
                return (log10(fwid) + 0.66)/0.27 
            if setting == 'continental':
                return (log10(fwid) + 0.75)/0.27
            else:
                return (log10(fwid) + 1.12)/0.33
        elif rake > 0.:
            # thrust/reverse
            if setting == 'oceanic':
                return (log10(fwid) + 1.79)/0.45
            if setting == 'continental':
                return (log10(fwid) + 1.83)/0.45
            else:
                return (log10(fwid) + 1.86)/0.46
        else:
            # normal
            return (log10(fwid) + 1.20)/0.36
        return 

    def get_median_area(self, mag, rake=None, setting=None):
        """
        Calculates median fault area from magnitude.
        """
        length = self.get_median_length(mag, rake, setting)
        width = self.get_median_width(mag, rake, setting)
        return length * width
