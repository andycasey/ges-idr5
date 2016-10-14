
"""
Plot things relevant to the ensemble model.
"""

import cPickle as pickle
import itertools
import logging
import numpy as np
import scipy.sparse

import matplotlib.pyplot as plt
from collections import OrderedDict
from matplotlib.ticker import MaxNLocator

import brewer2mpl # For pretty colors.

logger = logging.getLogger("ges")

_DEFAULT_SAVEFIG_KWDS = dict(dpi=300, bbox_inches="tight")


def systematic_uncertainty(model, ax=None, N_bins=50, xlabel=None, legend=True,
    **kwargs):
    """
    Plot the distribution of constant systematic uncertainties from all nodes.

    :param model:
        A fully-sampled ensemble model.

    """

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = ax.figure, ax

    N = model._chains["vs_c"].shape[1]

    colors = brewer2mpl.get_map("Set1", "qualitative", N).mpl_colors

    max_abs_sigma = np.sqrt(np.max(np.abs(model._chains["vs_c"])))
    bins = np.linspace(0, max_abs_sigma, N_bins)

    nodes = (model._metadata["node_ids"], model._metadata["node_names"])
    for i, (node_id, node_name) in enumerate(zip(*nodes)):

        ax.hist(np.sqrt(model._chains["vs_c"][:, i]),
            facecolor=colors[i], edgecolor=colors[i], bins=bins, 
            alpha=0.5, lw=2, normed=True, histtype="stepfilled",
            label=r"${{\rm {0}}}$".format(node_name.strip()))

    latex_labels = {
        "teff": r"${\rm Constant}$ ${\rm systematic}$ ${\rm uncertainty}$ "
                r"${\rm in}$ ${\rm effective}$ ${\rm temperature},$ "
                r"$\sigma_{c,T_{\rm eff}}$ $({\rm K})$",
        "logg": r"${\rm Constant}$ ${\rm systematic}$ ${\rm uncertainty}$ "
                r"${\rm in}$ ${\rm surface}$ ${\rm gravity},$ "
                r"$\sigma_{c,\log{g}}$ $({\rm dex})$",
        "feh":  r"${\rm Constant}$ ${\rm systematic}$ ${\rm uncertainty}$ "
                r"${\rm in}$ ${\rm metallicity},$ "
                r"$\sigma_{c,[{\rm Fe/H}]}$ $({\rm dex})$"
    }

    ax.set_xlabel(
        xlabel or latex_labels.get(model._parameter, model._parameter))
    ax.set_ylabel(r"${\rm Frequency}$")

    if legend:
        kwds = dict(frameon=False, ncol=2, loc="upper center")
        kwds.update(kwargs.get("legend_kwds", {}))
        ax.legend(**kwds)

    ax.set(
        adjustable="box-forced", 
        aspect=np.ptp(ax.get_xlim())/np.ptp(ax.get_ylim()))

    fig.tight_layout()
    # Monkey patch the savefig to incorporate the default keywords to ensure
    # the resulting figure is ready for publication.

    old_savefig = fig.savefig
    def new_savefig(self, *args, **kwargs):
        kwds = _DEFAULT_SAVEFIG_KWDS.copy()
        kwds.update(kwargs)
        return old_savefig(self, *args, **kwds)

    fig.savefig = new_savefig
    return fig



def node_relative_systematic_uncertainty(model, axes=None, quartiles=[16, 50, 84],
    ylims=(0.5, 10), **kwargs):
    """
    Plot the systematic uncertainty from all nodes as a function of the stellar
    parameters.
    """

    quartiles = np.hstack(quartiles)
    Q = len(quartiles)
    if Q not in (1, 3):
        raise ValueError("quartiles must be length 1 or 3")

    chains = model._chains
    assert chains is not None, "Has the model been sampled?"

    nodes = (model._metadata["node_ids"], model._metadata["node_names"])

    # TODO: Common colors.
    N = len(nodes[0])
    colors = brewer2mpl.get_map("Set1", "qualitative", N + 1).mpl_colors

    parameter_bounds = OrderedDict([
        ("teff", (3000, 8000)),
        ("logg", (0, 5)),
        ("feh", (-3.5, 0.5))
    ])

    if axes is None:
        fig, axes = plt.subplots(1, len(parameter_bounds))
    else:
        fig = axes[0].figure

    latex_xlabels = kwargs.get("latex_xlabels", dict(
        teff=r"${\rm Effective}$ ${\rm temperature},$ $T_{\rm eff}$ $({\rm K})$",
        logg=r"${\rm Surface}$ ${\rm gravity},$ $\log{g}$",
        feh=r"${\rm Metallicity},$ $[{\rm Fe/H}]$"))

    latex_ylabels = kwargs.get("latex_ylabels", dict(
        teff=r"$\sigma_{sys,T_{\rm eff}}/c_{sys}$",
        logg=r"$\sigma_{sys,\log{g}}/c_{sys}$",
        feh=r"$\sigma_{sys,[{\rm Fe/H}]}/c_{sys}$"))

    K = len(parameter_bounds)
    S = 500
    xs = np.linspace(0, 1, S)

    nticks = kwargs.get("nticks", 5)

    for i, (parameter, bounds) in enumerate(parameter_bounds.items()):

        x, ax = np.linspace(bounds[0], bounds[1], S), axes[i]

        for j, (node_id, node_name) in enumerate(zip(*nodes)):

            values = chains["vs_a"][:, i, j] * (1 - xs).reshape(-1, 1)**chains["vs_b"][:, i, j]

            sigma = np.sqrt(np.exp(values))
            q = np.percentile(sigma, quartiles, axis=1)

            ax.plot(x, q[Q / 2], lw=2, c=colors[j], zorder=10,
                label=r"${{\rm {0}}}$".format(node_name.strip()))
    
            # Show filled region.
            if Q > 1:
                ax.fill_between(x.flatten(), q[0], q[-1], facecolor=colors[j], 
                    alpha=0.5, edgecolor="none")


    # Set common ylimits.
    ylims = ylims or (0, np.max([ax.get_ylim() for ax in axes]))
    for i, (parameter, bounds) in enumerate(parameter_bounds.items()):

        ax = axes[i]

        ax.xaxis.set_major_locator(MaxNLocator(nticks))
        ax.yaxis.set_major_locator(MaxNLocator(nticks))

        ax.set_xlim(bounds)
        ax.set_ylim(ylims)
        if ax.is_first_col():
            ax.set_ylabel(latex_ylabels.get(model._parameter, model._parameter))

        else:
            ax.set_yticklabels([])

        if ax.is_last_row():
            ax.set_xlabel(latex_xlabels.get(parameter, parameter))

        else:
            ax.set_xticklabels([])

        # Force the axes to be square when saving it.
        ax.set(
            adjustable="box-forced", 
            aspect=np.ptp(ax.get_xlim())/np.ptp(ax.get_ylim()))

    fig.tight_layout()

    old_savefig = fig.savefig
    def new_savefig(self, *args, **kwargs):
        kwds = _DEFAULT_SAVEFIG_KWDS.copy()
        kwds.update(kwargs)
        return old_savefig(self, *args, **kwds)

    fig.savefig = new_savefig

    return fig


def node_uncertainty_with_snr(model, quartiles=[16, 50, 84], show_cr_bound=True, 
    xlims=(1, 100), ylims=None, Ns=100, **kwargs):
    """
    Plot the total node uncertainty as a function of S/N ratio.

    :param model:
        A trained ensemble model.

    :param xlims: [optional]
        A two-length tuple giving the lower and upper bounds to show in SNR.

    :param N: [optional]
        The number of draws to use when calculating the projected node uncertainty.
    """

    quartiles = sum([quartiles], [])
    Q = len(quartiles)
    if Q not in (1, 3):
        raise ValueError("quartiles must be length 1 or 3")

    chains = model._chains

    assert chains is not None, "Has the model been sampled?"

    x = np.linspace(xlims[0], xlims[1], 500).reshape(-1, 1)
    nodes = (model._metadata["node_ids"], model._metadata["node_names"])
    N = len(nodes[0])

    colors = brewer2mpl.get_map("Set1", "qualitative", N + 1).mpl_colors

    fig, ax = plt.subplots()

    for i, (node_id, node_name) in enumerate(zip(*nodes)):

        # Calculate total uncertainty as a function of SNR
        # (assuming no increase in systematic uncertainty due to different parts
        #  of parameter space)
        sigma = np.sqrt(
            chains["alpha_sq"][:, i].reshape(-1, 1) / x.T +
            chains["vs_c"][:, i].reshape(-1, 1) * np.ones_like(x.T))
        q = np.percentile(sigma, quartiles, axis=0)

        ax.plot(x, q[Q / 2], lw=2, c=colors[i], zorder=10,
            label=r"${{\rm {0}}}$".format(node_name.strip()))

        # Show filled region.
        if Q > 1:
            ax.fill_between(x.flatten(), q[0], q[-1], facecolor=colors[i], 
                alpha=0.5, edgecolor="none")


    if show_cr_bound:
        # Calculate the minimum variance as a function of SNR using all the
        # information from every node.
        
        C = model._chains["alpha_sq"].shape[0]
        sigma_wg = np.zeros((C, x.size))

        y = np.zeros((Ns, x.size))
            
        indices = np.random.choice(range(C), size=Ns, replace=False)
        for i, j in enumerate(indices):

            diag = np.sqrt(
                chains["alpha_sq"][j].reshape(-1, 1) / x.T +
                chains["vs_c"][j].reshape(-1, 1) * np.ones_like(x.T))

            I = np.eye(N)
            rho = np.dot(
                np.dot(I, chains["L_corr"][i]),
                np.dot(I, chains["L_corr"][i]).T)

            for k in range(x.size):
                
                Sigma = np.tile(diag[:, k], N).reshape(N, N) \
                      * np.repeat(diag[:, k], N).reshape(N, N) \
                      * rho

                W = np.ones((N, 1))
                Cinv = np.linalg.inv(Sigma)
                y[i, k] = 1.0/np.sqrt(np.dot(np.dot(W.T, Cinv), W))

        g = np.nanpercentile(y, quartiles, axis=0)

        ax.plot(x.flatten(), g[Q / 2], lw=2, c=colors[-1], linestyle="--",
            label=\
                r"${\rm Homogenised}$"
                "\n"
                r"$({\rm Cram\'er$-${\rm Rao}$ ${\rm bound)}$")
    
        if Q > 1:
            ax.fill_between(x.flatten(), g[0], g[-1], 
                facecolor=colors[-1], alpha=0.5, edgecolor="none", zorder=10)

    ax.set_xlim(*xlims)
    ax.set_xlabel(r"${\rm Signal}$-${\rm to}$-${\rm noise}$ ${\rm ratio},$ $S/N$ $({\rm pixel}^{-1})$")

    default_ylabels = {
        "teff": r"${\rm Uncertainty}$ ${\rm in}$ ${\rm effective}$ ${\rm temperature},$ $\sigma_{T_{\rm eff}}$ $({\rm K})$",
        "logg": r"${\rm Uncertainty}$ ${\rm in}$ ${\rm surface}$ ${\rm gravity},$ $\sigma_\log{g}$ $({\rm dex})$",
        "feh": r"${\rm Uncertainty}$ ${\rm in}$ ${\rm metallicity},$ $\sigma_{\rm [Fe/H]}$ $({\rm dex})$",
    }
    ylabel = kwargs.get("ylabel", default_ylabels.get(model._parameter,""))
    ax.set_ylabel(ylabel)

    if ylims is None:
        default_ylims = dict(teff=(0, 500), logg=(0, 1), feh=(0, 0.5))
        ylims = default_ylims.get(model._parameter, None)

    ax.set_ylim(ylims)
    
    legend_kwds = dict(fontsize=10, loc="upper center", ncol=2, frameon=False)
    legend_kwds.update(kwargs.get("legend_kwds", {}))
    plt.legend(**legend_kwds)
    
    ax.set(
        adjustable="box-forced", 
        aspect=np.ptp(ax.get_xlim())/np.ptp(ax.get_ylim()))

    fig.tight_layout()

    # Monkey patch the savefig to incorporate the default keywords to ensure
    # the resulting figure is ready for publication.

    old_savefig = fig.savefig
    def new_savefig(self, *args, **kwargs):
        kwds = _DEFAULT_SAVEFIG_KWDS.copy()
        kwds.update(kwargs)
        return old_savefig(self, *args, **kwds)

    fig.savefig = new_savefig

    return fig



def node_correlations(model, reorder=True, plot_edges=True, animate=False,
    **kwargs):
    """
    Show a lower-diagonal matrix coloured by the median correlation coefficient
    between each node.

    :param reorder: [optional]
        Re-order the correlation matrix to highlight structure.

    :param plot_edges: [optional]
        Plot edges surrounding the lower-diagonal matrix.

    :param animate: [optional]
        If `True`, then return a generator that produces `frames` of the 
        correlation matrix.
    """

    if not animate:
        return _node_correlations(
            model, reorder=reorder, plot_edges=plot_edges, **kwargs)
    else:
        return _node_correlations_animated(
            model, reorder=reorder, plot_edges=plot_edges, **kwargs)


def _node_correlations_animated(model, reorder=True, plot_edges=True, frames=100,
    reverse=True, **kwargs):

    node_names = np.array(model._metadata["node_names"])
    S, N, _ = model._chains["L_corr"].shape
    assert N == len(node_names)

    # Construct the correlation coefficients.
    I = np.eye(N)
    rho = np.nan * np.ones((S, N, N))

    for s in range(S):
        L = model._chains["L_corr"][s]
        rho[s] = np.dot(np.dot(I, L), np.dot(I, L).T)

    # Compress down to percentiles.
    if not reverse:
        q = np.linspace(0, 100, frames)
    
    else:
        q = np.hstack([
                np.linspace(0, 100, frames/2),
                np.linspace(100, 0, frames/2)
            ])

    F = q.size
    rho = np.percentile(rho, q, axis=0)

    if reorder:

        L_corr_median = np.median(model._chains["L_corr"], axis=0)
        rho_median = np.dot(np.dot(I, L_corr_median), np.dot(I, L_corr_median).T)

        permutations = list(itertools.permutations(range(N), N))
        score = np.nan * np.ones(len(permutations))
        for i, permutation in enumerate(permutations):

            _ = np.array(permutation)
            m = rho_median[:, _][_, :]

            score[i] = np.nansum(np.abs(np.diff(m, axis=0))) \
                     + np.nansum(np.abs(np.diff(m, axis=1)))

        matrix_indices = np.array(permutations[np.argmin(score)])

        # Check if we should flip this or not for visual aspects
        rho_ = rho_median[:, matrix_indices][matrix_indices, :]
        if np.sum(rho_[:, 0]) > np.sum(rho_[:, -1]):
            matrix_indices = matrix_indices[::-1]

        for f in range(F):
            rho[f] = rho[f][:, matrix_indices][matrix_indices, :]
        node_names = node_names[matrix_indices]

    # Common plotting parameters between frames.
    vrange = np.round(np.nanmax(np.abs(rho)), 1)

    assert vrange <= 1.0

    vmin = kwargs.get("vmin", None) or -vrange
    vmax = kwargs.get("vmax", None) or +vrange
    cmap = kwargs.get("cmap", "coolwarm")

    fig, ax = plt.subplots()

    # Monkey patch the savefig to incorporate the default keywords to ensure
    # the resulting figure is ready for publication.

    # Set NaN's for the upper triangle.
    rho = rho[:, 1:, :-1] 
    for f in range(F):
        for i in range(N - 1):
            for j in range(i + 1, N - 1):
                rho[f, i, j] = np.nan        

    im = ax.imshow(rho[0], vmin=vmin, vmax=vmax, cmap=cmap, interpolation="nearest")
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.xaxis.set_ticks_position("none") # Seriously matplotlib what the fuck
    ax.yaxis.set_ticks_position("none")

    ax.set_xticks(range(N - 1))
    ax.set_yticks(range(N - 1))
    ax.set_xticklabels(node_names[:-1], rotation=90)
    ax.set_yticklabels(node_names[1:])
    
    # Draw a line around everything.
    if plot_edges:
        kwds = dict(lw=2, c="k")

        edge = [-0.5, N - 1 - 0.5]
        ax.plot(edge, [N - 1 - 0.5, N - 1 -0.5], **kwds)
        ax.plot([-0.5, -0.5], edge, **kwds)

        x = []
        y = []
        for i in range(N - 1):
            x.extend([-0.5 + i, i + 0.5, i + 0.5])
            y.extend([i - 0.5, i - 0.5, i + 0.5])

        ax.plot(x, y, **kwds)
        tolerance = kwargs.get("__tolerance", 0.10)
        ax.set_xlim(-0.5 - tolerance, N - 1.5 + tolerance)
        ax.set_ylim(N - 1.5 + tolerance, -0.5 - tolerance)

    ax.set(adjustable="box-forced", aspect="equal")

    p = ax.get_position()
    cbar = plt.colorbar(ax=[ax], mappable=im, orientation="horizontal")
    cbar.set_ticks(np.linspace(-vrange, vrange, 5))
    cbar.ax.xaxis.set_ticks_position("none")

    default_labels = {
        "teff": r"${\rm Correlation}$ ${\rm coefficient},$ $\rho_{T_{\rm eff}}$",
        "logg": r"${\rm Correlation}$ ${\rm coefficient},$ $\rho_\log{g}$",
        "feh": r"${\rm Correlation}$ ${\rm coefficient},$ $\rho_{\rm [Fe/H]}$"
    }
    cbar.set_label(
        kwargs.get("label", default_labels.get(model._parameter, "rho")))

    fig.tight_layout()

    fig.subplots_adjust(bottom=0.35)
    plt.show()

    bottom, height = 0.05, 0.05
    co = kwargs.get("__cbar_offset", 0.01)
    cbar.ax.set_position([
        ax.get_position().x0 + co, 
        bottom, 
        ax.get_position().width - co*2, 
        bottom + height
    ])
    plt.show()

    """
    old_savefig = fig.savefig
    def new_savefig(self, *args, **kwargs):
        kwds = _DEFAULT_SAVEFIG_KWDS.copy()
        kwds.update(kwargs)
        return old_savefig(self, *args, **kwds)

    fig.savefig = new_savefig
    """
    
    for f in range(F):
        im.set_data(rho[f])
        #plt.show()
        yield fig


def _node_correlations(model, reorder=True, plot_edges=True, **kwargs):

    node_names = np.array(model._metadata["node_names"])
    S, N, _ = model._chains["L_corr"].shape
    assert N == len(node_names)

    # Construct the correlation coefficients from the Cholesky factors.
    I = np.eye(N)
    L = np.median(model._chains["L_corr"], axis=0)
    rho = np.dot(np.dot(I, L), np.dot(I, L).T)
    matrix_indices = np.arange(N, dtype=int)

    if reorder:
        permutations = list(itertools.permutations(range(N), N))
        score = np.nan * np.ones(len(permutations))
        for i, permutation in enumerate(permutations):

            _ = np.array(permutation)
            m = rho[:, _][_, :]

            score[i] = np.nansum(np.abs(np.diff(m, axis=0))) \
                     + np.nansum(np.abs(np.diff(m, axis=1)))

        matrix_indices = np.array(permutations[np.argmin(score)])

        # Check if we should flip this or not for visual aspects
        rho_ = rho[:, matrix_indices][matrix_indices, :]
        if np.sum(rho_[:, 0]) > np.sum(rho_[:, -1]):
            matrix_indices = matrix_indices[::-1]

        rho = rho[:, matrix_indices][matrix_indices, :]
        node_names = node_names[matrix_indices]

    # Get plotting parameters
    vrange = np.round(np.nanmax(np.abs(rho)), 1)

    vmin = kwargs.get("vmin", None) or -vrange
    vmax = kwargs.get("vmax", None) or +vrange
    cmap = kwargs.get("cmap", "coolwarm")

    fig, ax = plt.subplots()

    old_savefig = fig.savefig
    def new_savefig(self, *args, **kwargs):
        kwds = _DEFAULT_SAVEFIG_KWDS.copy()
        kwds.update(kwargs)
        return old_savefig(self, *args, **kwds)

    fig.savefig = new_savefig
    # Monkey patch the savefig to incorporate the default keywords to ensure
    # the resulting figure is ready for publication.

    # Set NaN's for the upper triangle.
    rho = rho[1:, :-1] 
    for i in range(N - 1):
        for j in range(i + 1, N - 1):
            rho[i, j] = np.nan        


    im = ax.imshow(rho, vmin=vmin, vmax=vmax, cmap=cmap, interpolation="nearest")
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.xaxis.set_ticks_position("none") # Seriously matplotlib what the fuck
    ax.yaxis.set_ticks_position("none")

    ax.set_xticks(range(N - 1))
    ax.set_yticks(range(N - 1))
    ax.set_xticklabels(node_names[:-1], rotation=90)
    ax.set_yticklabels(node_names[1:])
    
    # Draw a line around everything.
    if plot_edges:
        kwds = dict(lw=2, c="k")

        edge = [-0.5, N - 1 - 0.5]
        ax.plot(edge, [N - 1 - 0.5, N - 1 -0.5], **kwds)
        ax.plot([-0.5, -0.5], edge, **kwds)

        x = []
        y = []
        for i in range(N - 1):
            x.extend([-0.5 + i, i + 0.5, i + 0.5])
            y.extend([i - 0.5, i - 0.5, i + 0.5])

        ax.plot(x, y, **kwds)
        tolerance = kwargs.get("__tolerance", 0.10)
        ax.set_xlim(-0.5 - tolerance, N - 1.5 + tolerance)
        ax.set_ylim(N - 1.5 + tolerance, -0.5 - tolerance)

    ax.set(adjustable="box-forced", aspect="equal")

    p = ax.get_position()
    cbar = plt.colorbar(ax=[ax], mappable=im, orientation="horizontal")
    cbar.set_ticks(np.linspace(-vrange, vrange, 5))
    cbar.ax.xaxis.set_ticks_position("none")

    default_labels = {
        "teff": r"${\rm Correlation}$ ${\rm coefficient},$ $\rho_{T_{\rm eff}}$",
        "logg": r"${\rm Correlation}$ ${\rm coefficient},$ $\rho_\log{g}$",
        "feh": r"${\rm Correlation}$ ${\rm coefficient},$ $\rho_{\rm [Fe/H]}$"
    }
    cbar.set_label(
        kwargs.get("label", default_labels.get(model._parameter, "rho")))

    fig.tight_layout()

    fig.subplots_adjust(bottom=0.35)
    plt.show()

    bottom, height = 0.05, 0.05
    co = kwargs.get("__cbar_offset", 0.01)
    cbar.ax.set_position([
        ax.get_position().x0 + co, 
        bottom, 
        ax.get_position().width - co*2, 
        bottom + height
    ])
    plt.show()


    return fig


def biases(model, ax=None, N_bins=50, xlabel=None, legend=True, **kwargs):
    """
    Plot the distribution of biases from all nodes.

    :param model:
        A trained ensemble model.

    :param ax: [optional]
        The axes to plot the distributions.
    """

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = ax.figure, ax

    N = model._chains["biases"].shape[1]

    colors = brewer2mpl.get_map("Set1", "qualitative", N).mpl_colors

    max_abs_bias = np.max(np.abs(model._chains["biases"]))
    bins = np.linspace(-max_abs_bias, +max_abs_bias, N_bins)

    nodes = (model._metadata["node_ids"], model._metadata["node_names"])
    for i, (node_id, node_name) in enumerate(zip(*nodes)):

        ax.hist(model._chains["biases"][:, i], 
            facecolor=colors[i], edgecolor=colors[i], bins=bins, 
            alpha=0.5, lw=2, normed=True, histtype="stepfilled",
            label=r"${{\rm {0}}}$".format(node_name.strip()))

    latex_labels = {
        "teff": r"${\rm Bias}$ ${\rm in}$ ${\rm effective}$ ${\rm temperature},$ "
                r"$T_{\rm eff}$ $({\rm K})$",
        "logg": r"${\rm Bias}$ ${\rm in}$ ${\rm surface}$ ${\rm gravity},$ "
                r"$\log{g}$ $({\rm dex})$",
        "feh": r"${\rm Bias}$ ${\rm in}$ ${\rm metallicity},$ "
               r"$[{\rm Fe/H}]$ $({\rm dex})$"
    }

    ax.set_xlabel(
        xlabel or latex_labels.get(model._parameter, model._parameter))
    ax.set_ylabel(r"${\rm Frequency}$")

    if legend:
        kwds = dict(frameon=False, ncol=2, loc="upper center")
        kwds.update(kwargs.get("legend_kwds", {}))
        ax.legend(**kwds)

    ax.set(adjustable="box-forced")

    fig.tight_layout()
    # Monkey patch the savefig to incorporate the default keywords to ensure
    # the resulting figure is ready for publication.

    old_savefig = fig.savefig
    def new_savefig(self, *args, **kwargs):
        kwds = _DEFAULT_SAVEFIG_KWDS.copy()
        kwds.update(kwargs)
        return old_savefig(self, *args, **kwds)

    fig.savefig = new_savefig
    return fig



