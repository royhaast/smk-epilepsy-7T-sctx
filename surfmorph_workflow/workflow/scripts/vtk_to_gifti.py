import numpy as np
import pyvista as pv
import nibabel as nib

# IO files
fn_in = snakemake.input.vtk

# Read in surfaces again and find corresponding points 
S = pv.read(fn_in)

# Write metric to Gifti
gii = nib.gifti.GiftiImage()
gii.add_gifti_data_array(
    nib.gifti.GiftiDataArray(
        data=S.point_data['scalars'].astype(np.float32)
        )
) 

nib.save(gii, snakemake.output[0])