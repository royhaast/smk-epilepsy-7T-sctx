# Snakemake workflow: `smk-epilepsy-7T-sctx`
Snakemake workflow for processing BIDSified 7T MRI data as described in
"Multi-scale structural alterations of basal ganglia in focal epilepsy as demonstrated by 7T MRI"
(link here)

![Analysis pipeline](https://github.com/royhaast/smk-epilepsy-7T-sctx/blob/main/resources/pipeline.jpg?raw=true)

:minidisc: Inputs:
- participants.tsv file with target subject IDs for both controls (prefixed with C, e.g., sub-C001) and patients.
- For each subject:
  - BIDS dataset that includes
    - B1+ map
    - MP2RAGE data
  - Demographic and patient info

:computer: Singularity containers required:
 - Gradient-distortion correction (https://hub.docker.com/r/khanlab/gradcorrect)
 - fMRIPrep (https://hub.docker.com/r/nipreps/fmriprep)
 - Connectome Workbench (https://hub.docker.com/r/khanlab/connectome-workbench)
 - Surfmorph (https://hub.docker.com/r/khanlab/surfmorph)
 - FSL

:lock: Also required:
 - MATLAB

## Python packages that are used

- Numpy
- Pandas
- Nibabel
- PyVista

## External scripts that are used

- PreSurfer is developed by Sriranga Kashyap (@srikash). You can download and cite from here: 
https://doi.org/10.5281/zenodo.4626841.

- More info on the 7T-AMI atlas can be found in the original paper:
"Automatic segmentation of deep grey nuclei using a high-resolution 7T magnetic resonance imaging atlas - Quantification of T1 values in healthy volunteers" (https://doi.org/10.1111/ejn.15575)

- Surfmorph is developed by Ali Khan (@akhanf). You can find the original code here:
https://github.com/khanlab/surfmorph

## Authors

* Roy AM Haast @royhaast 
* Hugo Dary