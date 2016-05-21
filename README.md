
To-do:
- [X] Specify database structure
- [X] Ingest new results from FITS files (overwrite any existing results)
- [ ] After ingestion, create new figures at the node-level.
- [ ] If all node results are in, produce wg-level comparisons.

Node-level comparisons:
- FLAG distributions:
  - How many stars had parameters reported?
  - Which flags used, and how many?
  - Which flags are used in union with each other? (hard to display)
  - When a given node assigns a flag, what are the distributions of what other nodes say about that star/spectrum?

- HR diagram (colored by FEH) atop isochrones
- HR diagram (colored by FEH) atop all other nodes
- Cluster HR diagrams w/ FEH distributions:
  - NGC ???? (each above an isochrone)
- Setup-to-setup TEFF/LOGG/FEH/XI differences (same CNAME, different SETUP) for all nodes


WG-level comparisons:
- TEFF/LOGG/FEH/XI distributions:
  - corner plot
  - corner plot, highlighted by tech flags
  - delta distributions between this node and all others

- Distributions of TEFF/LOGG/FEH/XI for the benchmarks (box and whisker plots)
- RMS as a function of SNR for TEFF/LOGG/FEH/XI
- Reported error distributions in TEFF/LOGG/FEH/XI

# General checks
#1: Compare iDR4 results with the new results.
#2: What did you get for the Sun?
#3: Compare benchmark parameters with reference values.
#4: Extract open cluster and globular cluster members
#5: Review parameter distribution of field stars
#6: Review range of values for each results column.
#7: Compare results to photometric TEFF
#8: WG-specific checks/special cases


List of figures:
- [ ] Scatter plot showing iDR5 TEFF for a given node against the iDR4 RECOMMENDED value (1-to-1 and residual axis) [#1]
- [ ] Scatter plot showing iDR5 LOGG for a given node against the iDR4 RECOMMENDED value (1-to-1 and residual axis) [#1]
- [ ] Scatter plot showing iDR5 FEH for a given node against the iDR4 RECOMMENDED value (1-to-1 and residual axis) [#1]
- [ ] Scatter plot showing iDR5 XI for a given node against the iDR4 RECOMMENDED value (1-to-1 and residual axis) [#1]
- [ ] Scatter plot showing TEFF vs LOGG (colored by FEH) for all Solar spectra, highlighting the accepted position. [#2]
- [ ] Compare benchmark parameters with reference values as a box and whisker plot. [#3]
- [ ] Show a number of globular and open clusters compared to isochrones & accepted literature metallicity [#4]
- [ ] Show TEFF vs FEH, and TEFF vs LOGG and TEFF vs XI for globular/open cluster stars. [#4]
- [ ] Show LOGG vs FEH and LOGG vs TEFF and LOGG vs XI for globular/open cluster stars. [#4]
- [X] Histogram of TEFF [#5]
- [X] Histogram of LOGG [#5]
- [X] Histogram of FEH [#5]
- [X] Histogram of XI [#5]
- [X] HR diagram with metallicity color map [#5]
- [X] Histogram of E_TEFF [#6]
- [X] Histogram of E_LOGG [#6]
- [X] Histogram of E_FEH [#6]
- [X] Histogram of XI [#6]
- [ ] Compare TEFF to photometric TEFF [#7]
- [ ] Compare HR10|HR21 field with HR21 only stars in TEFF/LOGG [#8]
- [ ] Compare node-to-node parameters for WG13 [#8]

List of tables:
- [ ] Number of stars with valid TEFF/LOGG/FEH/XI values, and number of TECH/similar flags shown. [#1]
- [ ] Range of parameters reported.
