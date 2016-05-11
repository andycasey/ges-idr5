
"""
Node level figures for quality control.
"""

import matplotlib.pyplot as plt
import numpy as np
from collections import Counter



def flags_used(database, node, flag_column):
    """
    Return a histogram figure showing the distribution of flags used by a given
    node.

    :param database:
        A database for transactions.

    :param node:
        The node identifier. This can be a node identifier (integer) or a string
        that describes the node (e.g., 'wg11.vilnius').

    :param flag_column:
        The flag column to examine. This must be 'PECULI', 'REMARK', or 'TECH'.
    """

    flag_column = str(flag_column).lower()
    if flag_column not in ("peculi", "remark", "tech"):
        raise ValueError("flag column must be 'peculi', 'remark', or 'tech'")

    node_id = database.node_id(node)

    # Retrieve the data as a table.
    flags = database.retrieve_table(
        "SELECT cname, tech, peculi, remark FROM result WHERE node_id = %s",
        (node_id, ))

    count_flags = Counter(flags[flag_column])

    if not count_flags:
        raise WHOOOOANOFLAGS


    # Create the figure.
    fig, ax = plt.subplots()
    ax.hist()


    return fig




def hrd(database, node, where=None, isochrones=None, show_other_wg_nodes=False):
    """
    Return a Hertszprung-Russell Diagram (effective temperature versus surface
    gravity) for all stars reported by a single node.

    :param database:
        A database for transactions.

    :param node:
        The node identifier. This can be a node identifier (integer) or a string
        that describes the node (e.g., 'wg11.vilnius').

    :param where: [optional]
        Additional SQL constraints, which would start from '... and `where`'.

    :param isochrones: [optional]
        A single path or list of paths of isochrones to display on the figure.

    :param show_other_wg_nodes: [optional]
        If set to `True`, the results from other nodes in the same working
        group will be shown in the background to provide context.
    """

    node_id = database.node_id(node)

    where_str = "" if where is None else " and {}".format(where)
    if isochrones is not None and not isinstance(isochrones, (list, tuple)):
        isochrones = [isochrones]

    # Collect the results for this node.

    # Show 

