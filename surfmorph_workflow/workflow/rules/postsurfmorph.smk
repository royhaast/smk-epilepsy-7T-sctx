##TODO: fix this (again) to take into account all (i.e., rigid, affine) transforms..... buh
# Convert vtk to gifti
rule vtk_to_gifti_mni_space:
   input: 
       done = 'deriv/surfmorph_{atlas}/displacement/{roi}_participant1_iter3.done',
       vtk = 'deriv/surfmorph_{atlas}/displacement/sub-{subject}_{hemi}/anat/sub-{subject}_{hemi}_space-MNI152NLin2009cAsym_label-{roi}_surfmorphinout.vtk'
   output: 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-MNI152NLin2009cAsym_label-{roi}.surf.gii'
   params: 
       structure = lambda wildcards: 'CORTEX_LEFT' if ( wildcards.hemi == 'L' or wildcards.hemi == 'Lflip' ) else 'CORTEX_RIGHT' 
   group: 'postsurfmorph'  
   shell:
       """
       mris_convert {input.vtk} {output}
       wb_command -set-structure {output} {params.structure} -surface-type ANATOMICAL -surface-secondary-type GRAY_WHITE
       """ 

# Flip left surface back
rule unflip_gii_mni_space:
    input: 'deriv/postsurfmorph_{atlas}/sub-{subject}_Lflip/sub-{subject}_hemi-Lflip_space-MNI152NLin2009cAsym_label-{roi}.surf.gii'
    output: 'deriv/postsurfmorph_{atlas}/sub-{subject}_L/sub-{subject}_hemi-L_space-MNI152NLin2009cAsym_label-{roi}.surf.gii'
    params:
        structure = 'CORTEX_LEFT'
    group: 'postsurfmorph'
    shell:
        """
        wb_command -surface-flip-lr {input} {output}
        wb_command -set-structure {output} {params.structure} -surface-type ANATOMICAL -surface-secondary-type GRAY_WHITE
        """

rule vtk_to_gifti_subject_space:
   input: 
       done = 'deriv/surfmorph_{atlas}/displacement/{roi}_participant1_iter3.done',
       vtk = 'deriv/surfmorph_{atlas}/displacement/sub-{subject}_{hemi}/anat/sub-{subject}_{hemi}_space-T1w_label-{roi}_surfmorphinout.vtk'
   output: 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-aligned_label-{roi}.surf.gii'      
   params: 
       structure = lambda wildcards: 'CORTEX_LEFT' if wildcards.hemi == 'L' else 'CORTEX_RIGHT' 
   group: 'postsurfmorph'  
   shell:
       """
       mris_convert {input.vtk} {output}
       wb_command -set-structure {output} {params.structure} -surface-type ANATOMICAL -surface-secondary-type GRAY_WHITE
       """ 

# Flip left surface back
rule unflip_gii_subject_space:
    input: 'deriv/postsurfmorph_{atlas}/sub-{subject}_Lflip/sub-{subject}_hemi-Lflip_space-aligned_label-{roi}.surf.gii' 
    output: 'deriv/postsurfmorph_{atlas}/sub-{subject}_L/sub-{subject}_hemi-L_space-aligned_label-{roi}.surf.gii' 
    params:
        structure = 'CORTEX_LEFT'
    group: 'postsurfmorph'
    shell:
        """
        wb_command -surface-flip-lr {input} {output}
        wb_command -set-structure {output} {params.structure} -surface-type ANATOMICAL -surface-secondary-type GRAY_WHITE
        """

rule rigid_itk_inv_to_world:
    input: 'deriv/presurfmorph/sub-{subject}/MNI152NLin2009cAsym_to_sub-{subject}_T1w_type-itk.txt'
    output: 'deriv/presurfmorph/sub-{subject}/MNI152NLin2009cAsym_to_sub-{subject}_T1w_type-world.txt'
    group: 'postsurfmorph'   
    shell:
        """
        wb_command -convert-affine -from-itk {input} -to-world {output}
        """ 

# Bring surface back into subject space
rule apply_inv_rigid_gii_subject_space:
   input: 
       surf = 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-aligned_label-{roi}.surf.gii',
       xfm = 'deriv/presurfmorph/sub-{subject}/MNI152NLin2009cAsym_to_sub-{subject}_T1w_type-world.txt'
   output: 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-{roi}.surf.gii',
   group: 'postsurfmorph'      
   shell:
       """
       wb_command -surface-apply-affine {input.surf} {input.xfm} {output}  
       """

# Transfer surface displacement to gifti metric file
rule disp_to_gifti:
   input: 
       done = 'deriv/surfmorph_{atlas}/displacement/{roi}_participant1_iter3.done',
       vtk = lambda wildcards: 'deriv/surfmorph_{atlas}/displacement/sub-{subject}_{hemi}/anat/sub-{subject}_{hemi}_space-MNI152NLin2009cAsym_label-{roi}_surfmorphinout.vtk'.format(
            atlas=wildcards.atlas, subject=wildcards.subject, hemi='Lflip' if wildcards.hemi == 'L' else 'R', roi=wildcards.roi 
        ),
   output: 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_label-{roi}_displacement.shape.gii'
   group: 'postsurfmorph'    
   script: '../scripts/vtk_to_gifti.py'     

# Calculate distance between L-R surfaces
rule calculate_lr_distance:
   input: 
       surf_l = 'deriv/postsurfmorph_{atlas}/sub-{subject}_Lflip/sub-{subject}_hemi-Lflip_space-aligned_label-{roi}.surf.gii',
       surf_r = 'deriv/postsurfmorph_{atlas}/sub-{subject}_R/sub-{subject}_hemi-R_space-aligned_label-{roi}.surf.gii'
   output:
       distance = 'deriv/postsurfmorph_{atlas}/sub-{subject}/sub-{subject}_space-aligned_label-{roi}_lrdistance.shape.gii',
       vectors = 'deriv/postsurfmorph_{atlas}/sub-{subject}/sub-{subject}_space-aligned_label-{roi}_lrvector.shape.gii',
   group: 'postsurfmorph'      
   shell:
       """
       wb_command -surface-to-surface-3d-distance {input.surf_l} {input.surf_r} {output.distance} -vectors {output.vectors}  
       """

# Calculate displacement between L-R surfaces
rule calculate_lr_displacement:
   input: 
       surf_l = 'deriv/postsurfmorph_{atlas}/sub-{subject}_L/sub-{subject}_hemi-L_label-{roi}_displacement.shape.gii',
       surf_r = 'deriv/postsurfmorph_{atlas}/sub-{subject}_R/sub-{subject}_hemi-R_label-{roi}_displacement.shape.gii'
   output: 'deriv/postsurfmorph_{atlas}/sub-{subject}/sub-{subject}_space-aligned_label-{roi}_lrdisplacement.shape.gii'
   group: 'postsurfmorph'      
   shell:
       """
       wb_command -metric-math '(l-r)' {output} -var l {input.surf_l} -var r {input.surf_r}
       """

# Calculate curvature
rule calculate_curvature:
   input: 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-{roi}.surf.gii'
   output: 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_label-{roi}_curvature.shape.gii'
   group: 'postsurfmorph'
   shell:
       """
       wb_command -surface-curvature {input} -mean {output}
       wb_command -metric-math 'x' {output} -fixnan 0 -var x {output}
       """

# Calculate area
rule calculate_area:
   input: 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-{roi}.surf.gii'
   output: 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_label-{roi}_area.shape.gii'
   group: 'postsurfmorph'
   shell:
       """
       wb_command -surface-vertex-areas {input} {output}
       wb_command -metric-math 'x' {output} -fixnan 0 -var x {output}
       """

# Smooth surface maps
rule smooth_shape_maps:
   input:
       surf = 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-{roi}.surf.gii',
       metric = 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_label-{roi}_{metric}.shape.gii'
   output: 'deriv/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_label-{roi}_{metric}-smoothed.shape.gii'
   params:
       kernel = 2
   group: 'postsurfmorph'
   shell:
       """
       wb_command -metric-smoothing {input.surf} {input.metric} {params.kernel} {output} -method GEO_GAUSS_EQUAL -fix-zeros
       """
