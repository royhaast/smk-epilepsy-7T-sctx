import nibabel as nib
import pyvista as pv
import numpy as np

atlas      = nib.load(snakemake.input.atlas)
atlas_data = atlas.get_fdata()
atlas_ind  = np.argwhere(atlas_data)

affine     = atlas.affine
inv_affine = np.linalg.inv(atlas.affine)

S = pv.read(snakemake.input.surf)
S = S.compute_normals(cell_normals=False, auto_orient_normals=True)
v = S.points
n = S['Normals']

steps     = np.arange(0, 2, .1)
atlas_map = np.zeros((len(v)))
for vertex in np.arange(0,len(v)):
    label = 0
    for step in steps:
        vv_new = v[vertex][:3]+(n[vertex,:]*step)
        vv_new = np.dot(inv_affine, np.transpose(np.append(vv_new,[1])))
        vv_new = [int(np.round(vv_new[0])), int(np.round(vv_new[1])),int(np.round(vv_new[2]))]
        label  = atlas_data[vv_new[0],vv_new[1],vv_new[2]]  
        
        if label == 0:
            vv_new = v[vertex][:3]-(n[vertex,:]*step)
            vv_new = np.dot(inv_affine, np.transpose(np.append(vv_new,[1])))
            vv_new = [int(np.round(vv_new[0])), int(np.round(vv_new[1])),int(np.round(vv_new[2]))]       
            label  = atlas_data[vv_new[0],vv_new[1],vv_new[2]]    
        else:
            break
    atlas_map[vertex] = label

gii = nib.gifti.GiftiImage()
gii.add_gifti_data_array(
    nib.gifti.GiftiDataArray(
        data=atlas_map.astype(np.float32)
        )
) 

nib.save(gii, snakemake.output[0])
