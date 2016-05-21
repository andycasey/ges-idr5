
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
