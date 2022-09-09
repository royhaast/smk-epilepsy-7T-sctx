import nibabel as nib
import pandas as pd
import numpy as np

from pathlib import Path

def load_nifti(in_file):
    nii = nib.load(in_file)
    return nii.get_fdata()

# Load ROI labels and names
lut = pd.read_csv(snakemake.params.lut, delim_whitespace=True)

# Split ROI column to extract hemisphere and add to dataframe
hemi_info  = pd.DataFrame(lut['roi'].str.split('-',1).tolist(), columns = ['hemi','roi'])
lut['roi'] = hemi_info['roi']
lut.insert(loc=2, column='hemi', value=hemi_info['hemi'])

# All lh and rh labels
lh_labels = np.unique(lut[lut['hemi']=='Left']['label'].values)
rh_labels = np.unique(lut[lut['hemi']=='Right']['label'].values)

# Add additional multi-label ROIs
lut = lut.append(
    [{'label': [17,19,20,29], 'roi': 'pulvinar_total','group': 'posterior', 'hemi': 'Left', 'r': 0, 'g': 0, 'b':0, 'a': 0},
     {'label': [67,69,70,79], 'roi': 'pulvinar_total','group': 'posterior', 'hemi': 'Right', 'r': 0, 'g': 0, 'b':0, 'a': 0},
     {'label': lh_labels, 'roi': 'Whole', 'group': 'all', 'hemi': 'Left', 'r': 0, 'g': 0, 'b':0, 'a': 0},
     {'label': rh_labels, 'roi': 'Whole', 'group': 'all', 'hemi': 'Right', 'r': 0, 'g': 0, 'b':0, 'a': 0}],
    ignore_index=True
)

# Common path to derivatives
subject_path = '../deriv/{}/sub-{}'

# Load subjects info
subject_info = pd.read_csv(snakemake.params.subject_info)

# Filter subjects file
subject_info = subject_info.drop_duplicates(subset=['dicom'], ignore_index=True)
subject_info['bids'] = [ '{:03d}'.format(i) for i in subject_info['bids'] ]
subject_info = subject_info.loc[subject_info['bids'].isin(['{}'.format(snakemake.wildcards.subject)])]

# Load aseg stats output for eTIV
aseg = pd.read_csv(snakemake.input.stats,delimiter="\t")
aseg = aseg.set_index('Measure:volume')
etiv = aseg.loc['sub-{}'.format(snakemake.wildcards.subject)]['EstimatedTotalIntraCranialVol']

# Load thalamic parcellation
seg_file = nib.load(snakemake.input.thalamus)
seg_img = seg_file.get_fdata()

# Load T1 data
t1_img = load_nifti(snakemake.input.t1) 

# Load ventricles mask
mask_img = load_nifti(snakemake.input.mask)

# Mask out CSF voxels (i.e., voxels with T1 > 4 sec and ventricles mask)
t1_threshold = snakemake.params.t1_threshold
seg_img[t1_img>t1_threshold] = 0
seg_img[mask_img!=1] = 0

# If present, load other data
if Path(snakemake.params.dbm).is_file():
    jac = True

    # Load R2s and QSM data
    jac_img = load_nifti(snakemake.params.dbm)
else:
    jac = False  

# Get voxel volume
zooms = seg_file.header.get_zooms()
voxel_mm3 = np.prod(zooms)    

# Compute volume, mean R2s and T1 for each ROI
d = []
for l in range(0,len(lut)):
    volume_mm   = np.sum(np.isin(seg_img, lut.iloc[l]['label'], assume_unique=True)) * voxel_mm3
    volume_perc = (volume_mm / etiv) * 100
    t1     = np.nanmean(t1_img[np.isin(seg_img, lut.iloc[l]['label'], assume_unique=True)])

    if jac == True:
        dbm = np.nanmean(jac_img[np.isin(seg_img, lut.iloc[l]['label'], assume_unique=True)])
    else:
        dbm = -999

    d.append(('sub-'+snakemake.wildcards.subject,
                subject_info.iloc[0]['age'],
                subject_info.iloc[0]['sex'],
                etiv,
                subject_info.iloc[0]['type'],
                subject_info.iloc[0]['side'],
                lut.iloc[l]['label'],
                lut.iloc[l]['group'],
                lut.iloc[l]['roi'],
                lut.iloc[l]['hemi'],
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