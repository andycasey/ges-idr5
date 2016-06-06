
""" General utility functions. """

import os
from numpy import isfinite


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



def mh_or_feh(table):

    feh = sum(isfinite(table["feh"]))
    mh = sum(isfinite(table["mh"]))

    return "feh" if feh > mh else "mh"
