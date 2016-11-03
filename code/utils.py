
""" General utility functions. """

import logging
import os
from astropy.io import fits
from numpy import isfinite

logger = logging.getLogger("ges")

def safe_int(x, fill_value=-1):
    try:
        return int(x)
    except:
        return fill_value


def wg_as_int(wg):
    return int(str(wg).strip().lower().lstrip("wg"))


def parse_node_filename(filename):
    # GES_iDR5_WG10_NodeTemplate_QC_15052016
    _ = os.path.basename(filename).split("_")
    wg, node_name = _[2:4]
    if node_name.lower().endswith(".fits"):
        node_name = node_name[:-5]
    wg = wg_as_int(wg)
    return (wg, node_name)


def parse_wg_from_file(filename):
    """
    Parse the working group from a filename containing working group-level
    recommended files, and check the header to make sure that things are 
    consistent.

    :param filename:
        The local path of the file to check.        
    """

    basename = os.path.basename(filename)
    path_wg = basename.split("_")[2]

    # Check the WG entry in the image header.
    image = fits.open(filename)
    try:
        header_wg = image[0].header["NODE1"]

    except:
        logger.warn("Could not find NODE1 header in {}".format(filename))

    else:
        if path_wg != header_wg:
            raise ValueError("implicit WG from path is different from the explicit"\
                " WG in the image header ({} != {})".format(path_wg, header_wg))
    wg = int("{}".format(path_wg).lower().lstrip("wg"))
    return wg 



def mh_or_feh(table):

    feh = sum(isfinite(table["feh"]))
    mh = sum(isfinite(table["mh"]))

    return "feh" if feh > mh else "mh"
