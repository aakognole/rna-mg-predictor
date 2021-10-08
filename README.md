# rna-mg-predictor

This set of scripts can be used to setup a GCMC-MD simulation to identify the
Mg2+ and K+ binding sites in RNA molecules. It has been provided in the form to
be applied to native states of RNA, however it was originally designed to study
ion-atomsphere during the folding of RNA in following articles.

Kognole AA, MacKerell AD Jr., "Mg2+ Impacts the Twister Ribozyme Through Push-
Pull Stabilization of Non-Sequential Phosphate Pairs" Biophysical J. (2020)
118(6): 1424-1437 https://doi.org/10.1016/j.bpj.2020.01.021

Kognole AA,MacKerell AD Jr., "Contributions and Competition of Mg2+ and K+ in
Folding and Stabilization of the Twister Ribozyme" RNA (2020)
https://doi.org/10.1261/rna.076851.120


## You need to have

- GROMACS (https://github.com/gromacs/gromacs)
- OpenMM with CUDA (conda install -c conda-forge openmm)
- Plumed plugin for OpenMM (conda install -c conda-forge openmm-plumed)
- ParmEd to utilize the GROMACS formats in OpenMM (conda install -c omnia parmed)
- MMTSB Toolset to manupulate the PDB files (installed on fly from github)
- MDAnalysis to run clustering analysis ( conda install -c conda-forge mdanalysis)


