import os
import re
import numpy as np
import nibabel as nib

seg = nib.load(snakemake.input.seg)
seg_data = seg.get_fdata()

i = 0
for r, roi in enumerate(snakemake.params.rois):
    for h, (hemi, label) in enumerate(snakemake.params.structures[roi].items()):
        binarized = np.zeros(seg_data.shape)
        binarized[np.isin(seg_data,label)] = 1

        img = nib.Nifti1Image(binarized,affine=seg.affine,header=seg.header)
        fname = [ s for s in snakemake.output.seg if all(x in s for x in (roi, hemi)) ][0]
        nib.save(img,fname)
        os.system('antsApplyTransforms -d 3 -i {0} -r {1} -o {0} -t {2} -n MultiLabel -u int'.format(
            fname, snakemake.params.mni, snakemake.input.xfm
        ))        
        i += 1

t1w = nib.load(snakemake.input.t1w)
img = nib.Nifti1Image(t1w.get_fdata(), affine=t1w.affine, header=t1w.header) 
nib.save(img, snakemake.output.t1w)
os.system('antsApplyTransforms -d 3 -i {0} -r {1} -o {0} -t {2} -n BSpline -u int'.format(
    snakemake.output.t1w, snakemake.params.mni, snakemake.input.xfm
))