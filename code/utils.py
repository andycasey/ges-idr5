
""" General utility functions. """

import os

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
    wg = wg_as_int(wg)
    return (wg, node_name)
