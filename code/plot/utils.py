
import numpy as np
from astropy.table import Table

def parse_isochrone(filename):
    """
    Parse a PARSEC or Siess isochrone.

    :param filename:
        The filename of the isochrone.
    """

    is_parsec = "parsec" in filename.lower()
    kwds = {
        "format": "ascii"
    }
    if is_parsec:
        kwds.update({
            "data_start": 0,
            "header_start": 13
        })

    isochrone = Table.read(filename, **kwds)

    # Make columns common.
    if is_parsec:
        isochrone["logg"] = isochrone["logG"]
        isochrone["teff"] = 10**isochrone["logTe"]

    else: # Seiss
        isochrone["logg"] \
            = np.log10(isochrone["Mass"]) - np.log10(isochrone["R"]**2) + 4.44
        isochrone["teff"] = isochrone["Teff"]

    return isochrone



def wg_as_int(wg):
    return int(str(wg).strip().lower().lstrip("wg"))




def mh_or_feh(table):

    if table["mh"][0] is None:
        return "feh"
    elif table["feh"][0] is None:
        return "mh"
    feh = np.sum(np.isfinite(table["feh"]))
    mh = np.sum(np.isfinite(table["mh"]))

    return "feh" if feh > mh else "mh"
