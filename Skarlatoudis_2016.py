from math import log10

class Skarlatoudis2016_Interface():
    """
    A.A. Skarlatoudis, P.G. Somerville, H.K. Thio 2016. Source-scaling relations of interface
    subduction earthquakes for strong ground motion and tsunami simlation.
    Bulletin of the Seismological Society of America, 106(4), pp 1652-1662.

    Implements magnitude-length and length-magnitude scaling relationships,
    plus magnitude-width and width-magnitude scaling relationships.
    """

    def get_median_area(self, mag, rake=None):
        """
        Calculates median fault area from magnitude.
        These are from Table 4, non-self-similar!!
        """
        # log(area) = log(ca) + cb * M0(in Nm)
        # log(area) = log(ca) + cb * (1.5Mw + 9.05)
        # log(area) = (cb*1.5*Mw) + (log(ca)+(cb*9.05))
        return pow(10., (0.93 * mag) - 3.15347)


    def get_median_mag_from_area(self, area, rake=None):
        """
        Calculates median magnitude from fault area.
        """
        return (log10(area) - 3.15347)/0.93


    def get_median_length(self, mag, rake=None, setting=None):
        """
        Calculates median fault length from magnitude.
        """
        #area = self.get_median_area(mag)
        #width = self.get_median_width(mag)
        #return area/width
        return pow(10., (0.6285 * mag) - 2.801824)


    def get_median_mag_from_length(self, flen, rake=None):
        """
        Calculates median magnitude from fault length.
        """
        return log10((flen + 2.801824)/0.6285)


    def get_median_width(self, mag, rake=None):
        """
        Calculates median fault width from magnitude.
        """
        return pow(10., (0.3015 * mag) - 0.351646)


    def get_median_mag_from_width(self, fwid, rake=None):
        """
        Calculates median magnitude from fault width.
        """
        return (log10(fwid) + 0.351646)/0.3015

