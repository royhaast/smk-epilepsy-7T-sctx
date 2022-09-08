# Snakemake workflow: `smk-epilepsy-7T-sctx`
Snakemake workflow for processing BIDSified 7T MRI data as described in "title here" (link here)

Inputs:
- participants.tsv with target subject IDs
- For each target subject:
  - BIDS dataset
  - MP2RAGE data

Singularity containers required:
 - Gradient-distortion correction
 - Connectome Workbench

## Authors

* Roy AM Haast @royhaast 