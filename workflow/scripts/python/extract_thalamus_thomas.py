import nibabel as nib
import pandas as pd
import numpy as np
import os.path
from pathlib import Path

def load_nifti(in_file):
    nii = nib.load(in_file)
    return nii.get_fdata()

# Load ROI labels and names
lut = pd.read_csv(snakemake.params.lut)

# Add additional multi-label ROIs
lut = lut.append(
    [{'label': lut['label'].values, 'roi': 'SUBNUCLEI', 'group':'all', 'r': 0, 'g': 0, 'b':0, 'a': 0}],
    ignore_index=True
)

# Common path to derivatives
subject_path = '../deriv/{}/sub-{}'

# Load subjects info
subject_info = pd.read_csv(snakemake.params.subject_info)

# Filter subjects file
subject_info = subject_info.loc[subject_info['ID'].isin(['{}'.format(snakemake.wildcards.subject)])]

# Load aseg stats output for eTIV
aseg = pd.read_csv(snakemake.input.stats,delimiter="\t")
aseg = aseg.set_index('Measure:volume')
etiv = aseg.loc['sub-{}'.format(snakemake.wildcards.subject)]['EstimatedTotalIntraCranialVol']

# Load T1 data
t1_img = load_nifti(snakemake.input.t1)

# Load ventricles mask
mask_img = load_nifti(snakemake.input.mask)

# Mask out CSF voxels (i.e., voxels with T1 > 4 sec)
t1_threshold = snakemake.params.t1_threshold

# Load DBM data
jac_img = load_nifti(snakemake.input.dbm)

# Compute volume, mean T1 and DBM for each ROI
d = []
for h,hemi in enumerate(['Left','Right']):
    # Load thalamic parcellation
    seg_img = load_nifti(snakemake.input.thomas[h])
    seg_img[t1_img>t1_threshold] = 0
    seg_img[mask_img!=1] = 0

    # Load whole thalamus 
    thalamus_img = load_nifti(snakemake.input.thalamus[h])
    thalamus_img[t1_img>t1_threshold] = 0
    thalamus_img[mask_img!=1] = 0

    # Extract label image path
    label_path = os.path.dirname(snakemake.input.thomas[h])

    # Get voxel volume
    zooms = nib.load(snakemake.input.thomas[h]).header.get_zooms()
    voxel_mm3 = np.prod(zooms)    

    for l in range(0,len(lut)):
        volume_vox  = np.sum(np.isin(thalamus_img if l == 0 else seg_img, lut.iloc[l]['label'], assume_unique=True))

        if volume_vox != 0:
            volume_mm   = volume_vox * voxel_mm3
            volume_perc = (volume_mm / etiv) * 100
            t1          = np.nanmean(t1_img[np.isin(thalamus_img if l == 0 else seg_img, lut.iloc[l]['label'], assume_unique=True)])
            dbm         = np.nanmean(jac_img[np.isin(thalamus_img if l == 0 else seg_img, lut.iloc[l]['label'], assume_unique=True)])
        else:
            volume_mm   = -999
            volume_perc = -999
            t1          = -999
            dbm         = -999

        d.append(('sub-'+snakemake.wildcards.subject,
            subject_info.iloc[0]['AGE'],
            subject_info.iloc[0]['SEX'],
            etiv,
            subject_info.iloc[0]['TYPE'],
            subject_info.iloc[0]['SIDE'],
            lut.iloc[l]['label'],
            lut.iloc[l]['group'],
            lut.iloc[l]['roi'],
            hemi,
            volume_mm,
            volume_perc,
            t1,
            dbm)
        )

# Convert to Pandas dataframe
df = pd.DataFrame(d, columns=('Subject', 'Age', 'Sex', 'eTIV', 'Type', 'Side', 'Label', 'Label group', 'ROI', 'Hemisphere', 'Volume (mm)', 'Volume (%)', 'T1', 'DBM'))

# Define epilepsy side and add to dataframe
hemi_dict = {'Left':'L','Right':'R'}
epi_side  = []

for i, hemi in enumerate(df['Hemisphere']):
    contra_hemi = 'Right' if hemi == 'Left' else 'Left'

    if df.iloc[i]['Side'] == hemi_dict[hemi]:
        epi_tmp = 'Ipsilateral'
    elif df.iloc[i]['Side'] == 'Bilateral {}>{}'.format(hemi_dict[hemi],hemi_dict[contra_hemi]):
        epi_tmp = 'Ipsilateral'
    elif df.iloc[i]['Side'] == hemi_dict[contra_hemi]:
        epi_tmp = 'Contralateral'        
    elif df.iloc[i]['Side'] == 'Bilateral {}>{}'.format(hemi_dict[contra_hemi],hemi_dict[hemi]):
        epi_tmp = 'Contralateral'
    elif df.iloc[i]['Side'] == 'Bilateral R=L':
        epi_tmp = 'Bilateral'
    else:
        epi_tmp = 'NA'
    
    epi_side.append(epi_tmp)

df['Affected'] = epi_side

# Save to file
df.to_pickle(snakemake.output.pkl)
df.to_csv(snakemake.output.csv)