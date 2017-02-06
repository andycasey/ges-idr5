
""" Tasks to perform the top-level homogenisation for the Gaia-ESO Survey. """

import logging

logger = logging.getLogger("ges")


# DEFINE YOUNG AND OLD CLUSTERS
YOUNG_CLUSTER_NAMES = [
    "NGC6530",
    "Trumpler14",
    "IC2391",
    "IC2602",
    "IC4665",
    "NGC2547",
    "gamma2_Vel",
    "Cha_I",
    "NGC2264",
    "NGC3293",
    "Rho_Oph",
    "NGC2451"
]

OLD_CLUSTER_NAMES = [
    "Br25",
    "Br44",
    "Br81",
    "NGC2243",
    "NGC2516",
    "NGC3532",
    "NGC4815",
    "NGC6005",
    "NGC6633",
    "NGC6705",
    "NGC6802",
    "Pismis18",
    "Trumpler20",
    "Trumpler23",
]

def _follow_wg_provenance_rules(database, cname, column_names):
    """
    Return a row containing homogenised results for a single CNAME,
    by following the provenance rules specified by WG15 in the document at:

    https://docs.google.com/document/d/1_h8EgLwbTPqI1_269Q-QWudMqppHuHMoButrKmlV8rI/edit?ts=581b7497

    :param database:
        A database connection.

    :param cname:
        The CNAME of the object in question.

    :param column_names:
        A list of column names to retrieve for each object.
    """

    # Get all results for this object
    results = database.retrieve_table(
        """ SELECT  DISTINCT ON (wg_recommended_results.wg)
                    {column_names}
            FROM    spectra, wg_recommended_results
            WHERE   spectra.cname = '{cname}'
              AND   spectra.cname = wg_recommended_results.cname""".format(
        cname=cname, column_names=", ".join(column_names)))

    assert results is not None
    
    # Are there results from multiple WGs for this CNAME?
    if len(results) == 1:
        # No: there are only results from one WG
        # Provenance: Use this single result
        return results[0]

    else:
        # Yes: multiple results.
        # Is the star a calibrator (GES_TYPE is AR_SD_BM, GE_SD_BM, GE_SD_BC,
        # GE_SD_BW, GE_SD_RV, or GD_SD_PC)?

        if results["GES_TYPE"][0].strip() in \
        ("AR_SD_BM", "GE_SD_BM", "GE_SD_BW", "GE_SD_BC", "GE_SD_RV", "GD_SD_PC"):
            # Yes: it's a calibrator.


            if results["OBJECT"][0].strip() \
            in ("GJ205", "GJ436", "GJ526", "GJ551", "GJ581", "GJ699", "GJ880"):
                # Yes: this is an M-type star.
                # Provenance: Use the WG12 result.

                raise NotImplementedError

            else:
                # No: It's not an M-type star.
                # Did WG11 analyse this star?

                if 11 in results["wg"]:
                    # Yes: WG11 analysed this star
                    # Provenance: Use the WG11 result.

                    match = results["wg"] == 11
                    return results[match]

                else:
                    # No: WG11 did not analyse this star.
                    # Provenance: Use the WG10 result.
                    assert 10 in results["wg"]

                    # Do we have GES spectra *and* archival spectra?
                    raise NotImplementedError
                    # Take GES spectra over archival spectra


        else:
            # No: this star is not an M-type star

            # Is the star in a young cluster?
            # (Does it match any of the YOUNG_CLUSTER_NAMES?)
            if results["ges_fld"][0] in YOUNG_CLUSTER_NAMES:
                # Yes: this star is in a young cluster.

                # Did WG13 analyse this star, and does WG13 say that this star
                # has an effective temperature greater than 7000 K?
                if 13 in results["wg"] \
                and results["teff"][results["wg"] == 13][0] > 7000:
                    # Yes: WG13 analysed this star and they say that the
                    # temperature is greater than 7000 K.

                    # Provenance: Use the WG13 result.

                    match = results["wg"] == 13
                    return results[match]

                else:
                    # No: Either WG13 didn't analyse this star, or they did and
                    # they said the temperature is less than 7000 K.

                    # Provenance: Use the WG12 result

                    match = results["wg"] == 12
                    return results[match]

            else:
                # No: the star is not in a young cluster

                # Is the star in an old cluster?
                # (Does it match any of the OLD_CLUSTER_NAMES?)
                if results["ges_fld"][0] in OLD_CLUSTER_NAMES:
                    # Yes: this star is in an old cluster.
                    
                    # Was this star analysed by WG11?
                    if 11 in results["wg"]:
                        # Yes: this star was analysed by WG11.

                        match = results["wg"] == 11
                        return results[match]

                    else:
                        # No: this star was not analysed by WG11.

                        # Was this star observed with the HR9B setting AND is 
                        # the effective temperature of this star greater than
                        # 7000 K?

                        if any(["HR9B" in setup for setup in results["setup"]]) \
                        and any(results["teff"] > 7000):
                            # Yes, this star was observed with the HR9B setting
                            # AND at least one measurement of effective temperature
                            # says it is >7000 K.

                            # Provenance: Take the WG13 result

                            assert 13 in results["wg"]
                            match = results["wg"] == 13
                            return results[match]

                        else:
                            # No, this star was either not observed with the
                            # HR9B setting, or none of the temperature measurements
                            # are greater than 7000 K.

                            # Provenance: Take the WG10 result.

                            assert 10 in results["wg"]
                            match = results["wg"] == 10
                            return results[match]


                else:
                    # No, this star is not in an old cluster.
                    # It must be in the field.

                    # Was this star analysed by WG11?
                    if 11 in results["wg"]:
                        # Yes: this star was analysed by WG11.

                        # Provenance: Take the WG11 result.
                        match = results["wg"] == 11
                        return results[match]

                    else:
                        # No: this star was not analysed by WG11.

                        # Provenance: Take the WG10 result

                        match = results["wg"] == 10
                        return results[match]

    # You should never reach this part of the program.
    assert not ProgrammingError




def homogenise_survey_results(database):
    """
    Create a set of WG15/Survey level results based on the recommended results
    from the working groups.

    :param database:
        A database connection.
    """

    # Remove any previous entries?

    # Ensure that each cluster name actually exists.
    for cluster_name in YOUNG_CLUSTER_NAMES + OLD_CLUSTER_NAMES:
        check = database.retrieve(
            """ SELECT count(*) 
                FROM spectra
                WHERE ges_fld like '{}%'""".format(cluster_name))
        assert check is not None and check[0][0] > 1, \
            "Does the cluster '{}' exist?".format(cluster_name)
   
    # Get all column names.
    # HERE we select by WG12 because they will not have any provenance_columns   
    faux_table = database.retrieve_table(
        """ SELECT  DISTINCT ON (wg_recommended_results.wg)
                    spectra.*, wg_recommended_results.*
            FROM    spectra, wg_recommended_results
            WHERE   spectra.cname = wg_recommended_results.cname
              AND   wg_recommended_results.wg = '12'
            LIMIT   1
        """, prefixes=("spectra", "wg_recommended_results"))

    # Exclude provenance_* column names
    column_names = [c for c in faux_table.dtype.names if not c.startswith("provenance_")] 
    
    # Get all the unique CNAMEs for which we have results.
    cnames = database.retrieve("SELECT DISTINCT ON (cname) cname FROM results")
    N = len(cnames)

    for i, cname in enumerate(cnames):
        logger.info("Homogenising CNAME {0}/{1}: {2}".format(i, N, cname))

        result = _follow_wg_provenance_rules(database, cname[0], column_names)
        



    # Get the provenance for each star.

    # Update the table.


