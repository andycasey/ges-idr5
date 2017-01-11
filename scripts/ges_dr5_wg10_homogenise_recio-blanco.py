
"""
A script to homogenise the GES DR5 WG10 data using linear transformations to a
common node scale.
"""

import os
import numpy as np
import scipy.optimize as op
from astropy.table import Table

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


# Specify the WG10 node filenames to consider.
node_filenames = [
    "node-results/stellar-parameters/WG10/GES_iDR5_WG10_IAC.fits",
    "node-results/stellar-parameters/WG10/GES_iDR5_WG10_Lumba.fits",
    "node-results/stellar-parameters/WG10/GES_iDR5_WG10_Nice.fits"
]

# Specify a reference node filename.
reference_node_filename = "node-results/stellar-parameters/WG10/GES_iDR5_WG10_Lumba.fits"
assert reference_node_filename in node_filenames, \
    "The reference node filename must also be in the list of node filenames."
reference_node_index = node_filenames.index(reference_node_filename)

# Specify the output filename.
output_filename, overwrite = ("ges-dr5-wg10-recio-blanco.fits", True)

# ---------------------------------------------------------------------------- #

# Check that all node filenames are sorted in the same way, such that we can
# compare things on a row-by-row basis without worrying about mis-matches later.
node_submissions = [Table.read(node_filename) for node_filename in node_filenames]
for node_submission in node_submissions:
    node_submission.sort(["CNAME", "SETUP", "FILENAME"])

first_submission = node_submissions[0]
for node_submission in node_submissions[1:]:
    assert np.all(node_submission["CNAME"] == first_submission["CNAME"])
    assert np.all(node_submission["SETUP"] == first_submission["SETUP"])

# Some nodes use MH, some nodes use FEH. Make them all use FEH so that we can 
# just reference one column later on.
for node_filename, node_submission in zip(node_filenames, node_submissions):
    if not np.any(np.isfinite(node_submission["FEH"])) \
    and np.any(np.isfinite(node_submission["MH"])):
        print("Updating FEH to MH in {filename}".format(
            filename=os.path.basename(node_filename)))
        node_submission["FEH"] = node_submission["MH"]

# Perform homogenisation by different setups.
node_submissions = \
    [node_submission.group_by(["SETUP"]) for node_submission in node_submissions]

# Create a table that we will use to put in the homogenised values.
homogenised_table = first_submission.copy()

# Select a high-quality sub-sample of common stars between all nodes where the
# SNR > 15, and within reasonable ranges.

# Note that this can be different per SETUP. Any SETUPs not listed here will be
# ignored.

# TODO: In Alejandra's document she gives different requirements for HR10|HR21
#       compared to HR21 only.
all_transformation_subset_requirements = {
    "HR10|HR21": [
        # Column, lower range, upper range
        ("SNR", 15, np.inf),
        ("TEFF", 4000, 6500),
        ("LOGG", 0, 5),
        ("FEH", -2.0, +0.5)
    ],
    "HR21": [
        # Column, lower range, upper range
        ("SNR", 15, np.inf),
        ("TEFF", 4000, 6500),
        ("LOGG", 0, 5),
        ("FEH", -2.0, +0.5)
    ]
}

# Specify the parameters to homogenise.
parameters = ("TEFF", "LOGG", "FEH")

# We will reference the group indices from the first node submission.
group_indices = node_submissions[0].groups.indices
for index in range(group_indices.size - 1):
    gi, gj = group_indices[index:index + 2]

    setup = node_submissions[0]["SETUP"][gi]

    # Only consider the results from this setup.
    only_consider = np.zeros(len(first_submission), dtype=bool)
    only_consider[gi:gj] = True

    print("Running on SETUP = {}".format(setup))

    if setup not in all_transformation_subset_requirements:
        print("Skipping stars with SETUP = {0}".format(setup))

        # Put all the PARAMETER (and E_PARAMETER) columns as NaNs
        for parameter in parameters:
            homogenised_table[parameter][only_consider] = np.nan
            homogenised_table["E_{}".format(parameter)][only_consider] = np.nan

        # Continue to the next setup.
        continue

    transformation_subset = only_consider.copy()

    for column_name, lower, upper in all_transformation_subset_requirements[setup]:
        for node_filename, node_submission in zip(node_filenames, node_submissions):
            print("Checking {filename} for {column_name} values between {lower} and {upper}"\
                .format(
                    filename=os.path.basename(node_filename),
                    column_name=column_name,
                    lower=lower,
                    upper=upper))
            transformation_subset *= \
                  (upper > node_submission[column_name]) \
                * (node_submission[column_name] > lower) \
                * np.isfinite(node_submission[column_name])

    S = transformation_subset.sum()
    print("The transformation subset contains {0} rows".format(S))

    # OK, now calculate the linear transformations:
    #
    # \theta_N - \theta_{ref} = c_0 + c_1 * T_{{\rm eff},N} + c_2 * \log{g}_N
    #
    # where $\theta_{ref}$ stands for a parameter (teff, logg, or feh), as provided
    # by the reference node. The parameters with the $N$ subscript refers to those
    # of the node being translated into the reference scale.

    # The parameters to homogenise.
    parameters = ("TEFF", "LOGG", "FEH")

    # Get the reference node submission, because we will use it a lot later.
    reference_node_submission = node_submissions[reference_node_index]
    reference_node_name = reference_node_filename.split("_")[-1].split(".")[0]

    for i, (node_filename, node_submission) \
    in enumerate(zip(node_filenames, node_submissions)):

        for j, parameter in enumerate(parameters):

            if "CALIBRATED_{}".format(parameter) not in node_submission.dtype.names:
                node_submission["CALIBRATED_{}".format(parameter)] \
                    = np.nan * np.ones(len(node_submission))

            if i == reference_node_index:
                # This is the reference node.
                node_submission["CALIBRATED_{}".format(parameter)][only_consider] \
                    = node_submission[parameter][only_consider]
                continue

            # y = \theta_N - \theta_{ref}
            # x = [1, teff_N, logg_N]
            # c = [c0, c1, c2]

            # y = x \dot c

            x = np.array([
                np.ones(len(reference_node_submission)),
                node_submission["TEFF"],
                node_submission["LOGG"]
            ])
            y = node_submission[parameter] - reference_node_submission[parameter]

            # Only use the transformation subset for the fitting.
            x_subset = x[:, transformation_subset]
            y_subset = y[transformation_subset]

            # Specify the model, which is just the dot product.
            model = lambda x, *coeff: np.dot(coeff, x)

            # Find the coefficients and their uncertainties (the covariance matrix)
            coeff, cov = op.curve_fit(
                model, x_subset, y_subset, p0=np.zeros(x.shape[0]))

            # Apply the correction to all measurements by this node.
            node_submission["CALIBRATED_{}".format(parameter)][only_consider] \
                = (node_submission[parameter] - model(x, *coeff))[only_consider]


    # Calculate the node-to-node scatter (Section 2.3 of Alejandra's document)
    P = len(parameters)
    N = len(node_submissions)

    weights = {}

    for parameter in parameters:

        B = []
        A = []

        for i, ith_filename in enumerate(node_filenames):
            for j, jth_filename in enumerate(node_filenames[i + 1:]):

                node_i = node_submissions[i]
                node_j = node_submissions[i + 1 + j]

                delta = node_i[parameter] - node_j[parameter]

                i_idx = node_filenames.index(ith_filename)
                j_idx = node_filenames.index(jth_filename)

                _A = np.zeros(N)
                _A[i_idx] = 1
                _A[j_idx] = 1
                A.append(_A)
                B.append(np.nanstd(delta[only_consider]))

        # Solve for the node scatter terms.
        # (This will be of the same length as node_filenames)
        node_scatter = np.linalg.solve(np.array(A), np.array(B).T)

        # Calculate normalized weights.
        # w = \scatter/\sum{\scatter} (not \scatter^2 as naturally expected)
        weights[parameter] = node_scatter/node_scatter.sum()

    # Plot the node results before (black) and after (blue) calibration to make sure
    # that we are doing things correctly.

    limits = {
        "TEFF": (4000, 6500),
        "LOGG": (0, 5),
        "FEH": (-2.5, 0.5)
    }

    fig, axes = plt.subplots(P, N)
    axes = np.atleast_2d(axes)

    for n, (node_filename, node_submission) \
    in enumerate(zip(node_filenames, node_submissions)):

        for p, parameter in enumerate(parameters):

            ax = axes[p, n]

            # Before calibration.
            before_x = reference_node_submission[parameter][only_consider]
            before_y = node_submission[parameter][only_consider]

            # After calibration.
            after_x = reference_node_submission[parameter][only_consider]
            after_y = node_submission["CALIBRATED_{}".format(parameter)][only_consider]

            ax.scatter(before_x, before_y, 
                edgecolor="none", facecolor="#000000", s=1, zorder=1)
            ax.scatter(after_x, after_y, 
                edgecolor="none", facecolor="blue", s=1, zorder=10)

            lims = limits[parameter]
            ax.plot(lims, lims, c="#666666", zorder=-1)
            ax.set_xlim(lims)
            ax.set_ylim(lims)

            # Get the node name from the filename.
            node_name = node_filename.split("_")[-1].split(".")[0]

            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            
            ax.set_xlabel("{reference_node_name} {parameter}".format(
                reference_node_name=reference_node_name, parameter=parameter))

            ax.set_ylabel("{node_name} {parameter}".format(
                node_name=node_name, parameter=parameter))

            if ax.is_first_row():
                ax.set_title(node_name)
       
            if not ax.is_first_col():
                ax.set_yticklabels([])

    fig.tight_layout()
    fig.savefig(
        "ges-dr5-wg10-{0}-recio-blanco-calibration-comparison.pdf".format(setup))

    # Produce a homogenised set of values using the pre-calculated weights.
    for parameter, w in weights.items():
        x = np.array([
            node_submission["CALIBRATED_{}".format(parameter)][only_consider] \
            for node_submission in node_submissions])

        homogenised_table[parameter][only_consider] = np.sum(x.T * w, axis=1)

        # As per the report, the errors are taken as the standard deviation of
        # the results.
        homogenised_table["E_{}".format(parameter)][only_consider] = np.std(x, axis=0)
        
# WRITE IT
homogenised_table.write(output_filename, overwrite=overwrite)
print("Homogenised table written to {}".format(output_filename))