from math import log10, sqrt

class NZNSHMcrustal_2022():
    """
    M. Stirling, M. Fitzgerald, B. Shaw, C. Ross 2022. New Magnitude–Area Scaling Relations for 
    the New Zealand National Seismic Hazard Model 2022
    Bulletin of the Seismological Society of America, 114(1), pp 137-149.

    Implements magnitude-area and area-magnitude scaling relationships,
    plus magnitude-length and length-magnitude scaling relationships and
    plus magnitude-width and width-magnitude scaling relationships using
    an aspect ratio relation.
    """
    def faultAR(self, mag):
        """
        # Chiou & Youngs 2006 for aspect ratio for strike-slip crustal only
        # log(aspect ratio) = 0.01752 * pow(mag - 4, 3.097)
        """
        if mag < 4.:
            ar = 1.
        else:
            #ar = pow(10., 0.01752 * pow(mag - 4., 3.097)) # ss
            ar = pow(10., (0.01752-0.003) * pow(mag - 4., 3.097)) # modified SS/rev
            #ar = pow(10., (0.01752-0.01099) * pow(mag - 4., 3.097)) # reverse
        return ar

    def get_median_area(self, mag, rake=None):
        """
        Calculates median fault area from magnitude.
        """
        return pow(10., mag - 4.2)


    def get_median_mag_from_area(self, area, rake=None):
        """
        Calculates median magnitude from fault area.
        """
        return log10(area) + 4.2


    def get_median_length(self, mag, rake=None, setting=None):
        """
        Calculates median fault length from magnitude.
        """
        return sqrt(self.get_median_area(mag) * self.faultAR(mag))


    def get_median_mag_from_length(self, flen, rake=None):
        """
        Calculates median magnitude from fault length.
        """
        return None


    def get_median_width(self, mag, rake=None):
        """
        Calculates median fault width from magnitude.
        """
        return sqrt(self.get_median_area(mag) / self.faultAR(mag))


    def get_median_mag_from_width(self, fwid, rake=None):
        """
        Calculates median magnitude from fault width.
        """
        return None

