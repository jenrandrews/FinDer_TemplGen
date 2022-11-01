from math import log10
from openquake.hazardlib.scalerel.strasser2010 import StrasserInterface

class Strasser2010_Interface_Ext(StrasserInterface):
    """
    Strasser, F.O., M.C. Arango, and J.J. Bommer 2010. Scaling of the source dimensions of 
    interface and intraslab subduction-zone earthquakes with moment magnitude.
    Seismological Research Letters, 81(6), pp 941-950.

    Extends the base class with magnitude-length and length-magnitude scaling relationships,
    plus magnitude-width and width-magnitude scaling relationships.
    """

    def get_median_length(self, mag, rake=None):
        """
        Calculates median fault length from magnitude.
        """
        return pow(10.0, (0.585 * mag) - 2.477)


    def get_median_mag_from_length(self, flen, rake=None):
        """
        Calculates median magnitude from fault length.
        """
        return 4.868 + (1.392 * log10(flen))


    def get_median_width(self, mag, rake=None):
        """
        Calculates median fault width from magnitude.
        """
        return pow(10.0, (0.351 * mag) - 0.882)


    def get_median_mag_from_width(self, fwid, rake=None):
        """
        Calculates median magnitude from fault width.
        """
        return 4.410 + (1.805 * log10(fwid))

