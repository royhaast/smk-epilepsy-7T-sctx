import numpy as np
import nibabel as nib
from scipy import stats

tmp = nib.load(snakemake.input[0]).darrays[0].data
nvertices = len(tmp)

atlas_data = np.zeros((nvertices,len(snakemake.input)))
for i, in_gii in enumerate(snakemake.input):
    atlas = nib.load(in_gii)                
    atlas_data[:,i] = atlas.darrays[0].data

# Replace non-ROI labels with NaN
atlas_data[~np.isin(atlas_data, snakemake.params.labels)] = np.nan

m = stats.mode(atlas_data, axis=1, nan_policy='omit')
parcellation = m[0][:,0]

gii = nib.gifti.GiftiImage()
gii.add_gifti_data_array(
    nib.gifti.GiftiDataArray(
        data=parcellation.astype(np.float32)
        )
) 

nib.save(gii, snakemake.output[0])
