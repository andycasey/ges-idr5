

import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from matplotlib.colors import LogNorm

import scipy.sparse.csgraph


_WG14_NODE_IDS = {
    "01": "Arcetri",
    "02": "CAUP",
    "03": "EPINARBO",
    "04": "IAC",
    "05": "Lumba",
    "06": "MaxPlanck",
    "07": "MyGIsFOS",
    "08": "Nice",
    "09": "OACT",
    "10": "OAPA",
    "11": "UCM",
    "12": "ULB",
    "13": "Vilnius",
    "14": "GSSP",
    "15": "IAC",
    "16": "Liege",
    "17": "MGNDU",
    "18": "Mntp",
    "19": "ON",
    "20": "ROB",
    "21": "ROBGrid",
    "22": "BIN",
    "23": "Halpha",
    "24": "NBfilters",
    "25": "TSNE"
}

def heatmap(database, wg, kind="tech", show_multiple_flags_per_node=True,
    group_by=None, use_cuthill_mckee=True, minimum_entry=1, figsize=(12, 12),
    **kwargs):
    """
    Plot a heat map of the flags used by a working group. A 2D heat map is
    created showing the frequency of flags used by different nodes. Flag entries 
    are only shown if they are used by at least two nodes.

    :param database:
        A database to connect to.

    :param wg:
        The working group (e.g., 10).

    :param kind: [optional]
        The kind of flag to show. Options include 'tech', 'remark', and 'peculi'

    :param show_multiple_flags_per_node: [optional]
        If set to `True`, then multiple flags used by the same node will be
        shown in the figure. If `False`, then only flags used by two different
        nodes will be shown.

    :param group_by: [optional]
        Group the grid. Options include: `None` for no grouping, `nodes` to
        group by nodes, and `issue` to group by issues.

    :param use_cuthill_mckee: [optional]
        Use the Cuthill-McKee algorithm to sort the flags to maximize symmetry
        and minimize the diagonal bandwidth of the matrix.

    :param minimum_entry: [optional]
        The minimum number of acceptable entries for any two dimensional grid
        bin. Any bins with values less than `minimum_entry` will not be shown.

    :param figsize: [optional]
        The size of the figure in inches `(xsize, ysize)`.
    """

    # Ensure to strip on | because some nodes give trailing | entries.
    _strip_chars = " |"

    kind = kind.lower()
    kind_available = ("tech", "remark", "peculi")
    if kind not in kind_available:
        raise ValueError("kind '{}' not recognized in: {}".format(
            kind, ", ".join(kind_available)))

    if group_by is not None:
        group_by = group_by.lower()
        group_by_available = ("issue", "node")
        if group_by not in group_by_available:
            raise ValueError("ordering by '{}' not available: {}".format(
                group_by, ", ".join(group_by_available)))

    # Select everything with a non-zero flag entry
    results = database.retrieve_table(
        """ WITH n as (
                SELECT id FROM nodes WHERE wg = {wg})
            SELECT r.node_id, r.cname, r.{kind}
            FROM n, results as r
            WHERE TRIM(r.{kind}) <> ''
              AND TRIM(r.{kind}) <> 'NaN'
              AND r.node_id = n.id
        """.format(wg=wg, kind=kind)) 

    flat_flags \
        = sum([_.strip(_strip_chars).split("|") for _ in results[kind]], [])

    # Unique issue id numbers.
    issue_ids = np.sort(np.unique([each.split("-")[0] for each in flat_flags]))

    # Node ids.
    node_ids = np.sort(np.unique([each.split("-")[2] for each in flat_flags]))

    L, M = len(issue_ids), len(node_ids)
    Z = np.zeros((L * M, L * M), dtype=int)

    if group_by == "node" or group_by is None:
        get_index = lambda issue, node: \
            (np.where(node_ids == node)[0][0] * len(issue_ids) \
                + np.where(issue_ids == issue)[0][0])
        labels = np.tile(issue_ids, M)

    elif group_by == "issue":
        get_index = lambda issue, node: \
            (np.where(issue_ids == issue)[0][0] * len(node_ids)) \
                + np.where(node_ids == node)[0][0]
        labels = np.repeat(issue_ids, M)

    else:
        raise ValueError("sorting by '{}' not available".format(group_by))

    # Group results by CNAME.
    for group in results.group_by("cname").groups:

        if isinstance(group[kind], (str, unicode)):
            group_flags = group[kind].strip(_strip_chars).split("|")
        else:
            group_flags \
                = sum([_.strip(_strip_chars).split("|") for _ in group[kind]], [])

        if len(group_flags) == 1:
            continue

        # Remove the WG and confidence entry.
        flags = np.unique(
            ["-".join([_.split("-")[0], _.split("-")[2]]) for _ in group_flags])

        for x, y in combinations(flags, 2):

            x_issue, x_node = x.split("-")
            y_issue, y_node = y.split("-")

            if x_node == y_node and not show_multiple_flags_per_node: continue

            xi = get_index(x_issue, x_node)
            yi = get_index(y_issue, y_node)

            Z[xi, yi] += 1
            Z[yi, xi] += 1


    before_reorder = np.sort(labels[np.where(Z == Z.max())[0]])

    # Show structure within each node?
    if use_cuthill_mckee:

        if group_by is None:

            matrix_indices = scipy.sparse.csgraph.reverse_cuthill_mckee(
                scipy.sparse.csr_matrix(Z))
            Z = Z[matrix_indices, :][:, matrix_indices]

            labels = np.array([
                ["-".join([iid, nid]) for iid in issue_ids] for nid in node_ids])
            labels = labels.flatten()[matrix_indices]


        elif group_by == "node":

            # Re-order the issue indices so that they show structure, but keep
            # the same issue order for each of the nodes.

            """
            matrix_indices = scipy.sparse.csgraph.reverse_cuthill_mckee(
                scipy.sparse.csr_matrix(Z))

            issue_indices = np.array(
                [_ for i, _ in enumerate(matrix_indices % L) \
                    if _ not in (matrix_indices % L)[:i]]).flatten().astype(int)

            for i in range(M):
                labels[L*i:L*(i + 1)] = labels[L*i + issue_indices]
                Z[L*i:L*(i + 1), :] = Z[L*i + issue_indices, :]
                Z[:, L*i:L*(i + 1)] = Z[:, L*i + issue_indices]
            """

            # Re-order the issue indices so that they show structure

            for i in range(M):
                indices = scipy.sparse.csgraph.reverse_cuthill_mckee(
                    scipy.sparse.csr_matrix(Z[L*i:L*(i+1), L*i:L*(i+1)]))

                labels[L*i:L*(i + 1)] = labels[L*i + indices]
                Z[L*i:L*(i + 1), :] = Z[L*i + indices, :]
                Z[:, L*i:L*(i + 1)] = Z[:, L*i + indices]


        else:
            raise NotImplementedError 


    after_reorder = np.sort(labels[np.where(Z == Z.max())[0]])
    assert np.all(before_reorder == after_reorder)

    keep = np.any(Z >= minimum_entry, axis=1)
    Z = Z[keep, :][:, keep]
    labels = labels[keep]
    
    gridlines = np.array([sum(keep[L*i:L*(i + 1)]) for i in range(M)])


    fig, ax = plt.subplots(figsize=figsize)
    kwds = dict(cmap="viridis", aspect="auto", interpolation="nearest", norm=LogNorm())
    kwds.update(**kwargs)

    ax.imshow(np.eye(*Z.shape), cmap="Greys", vmin=0, vmax=2, interpolation="nearest")
    image = ax.imshow(Z, **kwds)

    ax.set_xticks(np.arange(Z.shape[0]))
    ax.set_yticks(np.arange(Z.shape[0]))

    # Put gridlines
    if group_by is not None and gridlines is None:
        a, b = (L, M) if group_by == "issue" else (M, L)
        for _ in range(a):
            ax.axhline(_*b - 0.5, c="#000000", linewidth=2)
            ax.axvline(_*b - 0.5, c="#000000", linewidth=2)

    if gridlines is not None:
        for _ in np.hstack([0, np.cumsum(gridlines)]):
            ax.axhline(_ - 0.5, c="#000000", linewidth=2)
            ax.axvline(_ - 0.5, c="#000000", linewidth=2)

        
        ticks = np.hstack([0, np.cumsum(gridlines)])
        ticks = 0.5 * np.diff(ticks) + ticks[:-1] - 0.5
        ax_right = ax.twinx()
        ax_right.set_yticks(ticks)
        ax_right.set_yticklabels(
            [_WG14_NODE_IDS.get(node_id, node_id) for node_id in node_ids],
            verticalalignment="center")
        ax_right.set_ylim(ax.get_ylim())
        ax_right.tick_params(width=0)


    ax.set_xlim(-0.5, Z.shape[0] - 0.5)
    ax.set_ylim(Z.shape[0] - 0.5, -0.5)

    ax_right.set_ylim(ax.get_ylim())

    labels = labels if labels is not None else np.tile(issue_ids, M)
    ax.set_yticklabels(labels)
    ax.set_xticklabels(labels, rotation=90)

    ax.tick_params(width=0)

    fig.tight_layout()

    _height, _space = 0.025, 0.025
    cbar = plt.colorbar(image,
        cax=fig.add_axes([
            fig.subplotpars.left, 
            1 - _height - _space,
            fig.subplotpars.right - fig.subplotpars.left, 
            _height
            ]),
        orientation="horizontal")
    fig.subplots_adjust(top=1 - _height - 2*_space)

    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.xaxis.set_label_position("top")

    cbar.ax.xaxis.set_tick_params(width=0)
    # Empty metadata dict.
    meta = {}

    return (fig, meta)
    
