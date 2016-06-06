
""" Produce summary tables. """


import numpy as np
from collections import Counter
from astropy.table import Table

import utils

def stellar_parameter_range(database, wg=None):
    """
    Produce a summary table outlining the range of stellar parameters reported.

    :param database:
        A database for transactions.

    :param wg: [optional]
        A specific working group.
    """

    if wg is None:
        nodes = database.retrieve("SELECT id, wg, name FROM nodes")
    else:
        nodes = database.retrieve("SELECT id, wg, name FROM nodes WHERE wg = %s",
            (utils.wg_as_int(wg), ))

    rows = []
    for node_id, node_wg, node_name in nodes:
        results = database.retrieve_table(
            """SELECT teff, e_teff, logg, e_logg, mh, e_mh, xi, e_xi
               FROM results WHERE node_id = %s""", (node_id, ))

        name = "WG{}/{}".format(node_wg, node_name) if wg is None else node_name
        if results is None or len(results) == 0:
            rows.append([name] + [np.nan] * 16)
            continue

        row = [name]
        for column in ("teff", "logg", "mh", "xi"):
            for column in [column, "e_{}".format(column)]:
                row.extend([
                    np.nanmin(results[column]),
                    np.nanmax(results[column])
                ])
        rows.append(row)

    return Table(rows=rows, names=("Name", 
        "Min. TEFF", "Max. TEFF", "Min. E_TEFF", "Max. E_TEFF",
        "Min. LOGG", "Max. LOGG", "Min. E_LOGG", "Max. E_LOGG",
        "Min. MH", "Max. MH", "Min. E_MH", "Max. E_MH",
        "Min. XI", "Max. XI", "Min. E_XI", "Max. E_XI"))


def stellar_parameter_summary(database, wg=None):
    """
    Produce a summary table outlining the number of valid results produced by
    different nodes.

    :param database:
        A database for transactions.

    :param wg: [optional]
        The working group.
    """

    # Get node ids.
    if wg is None:
        nodes = database.retrieve("""SELECT id, wg, name FROM nodes""")
    else:
        nodes = database.retrieve(
            """SELECT id, wg, name FROM nodes WHERE wg = %s""",
            (utils.wg_as_int(wg), ))

    rows = []
    for node_id, node_wg, node_name in nodes:
        results = database.retrieve_table(
            """SELECT teff, e_teff, logg, e_logg, mh, e_mh, xi, e_xi,
                      tech, remark, peculi
               FROM results
               WHERE node_id = %s""", (node_id, ))

        name = "WG{}/{}".format(node_wg, node_name) if wg is None else node_name
        if results is None or len(results) == 0:
            rows.append([name] + [0] * (4 * 2 + 3 + 1))
            continue

        N = len(results)
        valid = []
        for column in ("teff", "logg", "mh", "xi"):
            valid.extend([
                np.isfinite(results[column]).sum(),
                np.isfinite(results["e_{}".format(column)]).sum()
            ])

        num_flags = []
        for column in ("tech", "remark", "peculi"):
            num_flags.append(np.sum(results[column] != ""))

        row = [name, N]
        row.extend(valid)
        row.extend(num_flags)
        rows.append(row)

    return Table(rows=rows, names=("Node", "Number of results", "Valid TEFF",
        "Valid E_TEFF", "Valid LOGG", "Valid E_LOGG", "Valid MH", "Valid E_MH",
        "Valid XI", "Valid E_XI", "TECH entries", "REMARK entries", "PECULI entries"))



def tech_flags(database, wg, node_name, column="TECH"):
    """
    Produce a summary table outlining the number of times certain flags were
    used.

    :param database:
        A database for transactions.

    :param wg:
        The working group.

    :param node_name:
        The name of the node to summarize results for.

    """

    # FLAG / TOTAL_COUNT


    node_id = database.retrieve_node_id(wg, node_name)

    flags = database.retrieve_table(
        """SELECT {} FROM results WHERE node_id = %s""".format(column),
        (node_id, ))

    if flags is None: return None

    counts = Counter([each.strip() \
        for each in "|".join(flags["tech"]).split("|") if each.strip()])
    rows = sorted(counts.iteritems(), key=lambda (k, v): v)[::-1]
    
    return Table(rows=rows, names=("{} FLAG".format(column.upper()), "N"))
