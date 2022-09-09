# Snakemake workflow: `smk-epilepsy-7T-sctx`
Snakemake workflow for processing BIDSified 7T MRI data as described in "title here" (link here)

:minidisc: Inputs:
- participants.tsv file with target subject IDs
- For each target subject:
  - BIDS dataset
  - MP2RAGE data

:computer: Singularity containers required:
 - Gradient-distortion correction
 - Connectome Workbench

## External scripts that are used

- PreSurfer is developed by Sriranga Kashyap (@srikash). You can download and cite from here: 
https://doi.org/10.5281/zenodo.4626841.

- More info on the 7T-AMI atlas can be found in the original paper:
"Automatic segmentation of deep grey nuclei using a high-resolution 7T magnetic resonance imaging atlas - Quantification of T1 values in healthy volunteers" (https://doi.org/10.1111/ejn.15575)

## Authors

* Roy AM Haast @royhaast 
* Hugo Dary