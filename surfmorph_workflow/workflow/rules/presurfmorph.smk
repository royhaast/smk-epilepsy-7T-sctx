# Prep THOMAS data for surface-based morphometry
def define_thomas_input(wildcards, atlas=config['atlas']):
    subject = '{wildcards.subject}'.format(wildcards=wildcards)
    group   = 'controls' if subject[0] == 'C' else 'patients'

    if atlas == 'thomas':
        return {
            't1w': join(config['deriv'][group],'subcortical_atlas/sub-{subject}/THOMAS_T1MAP_MASKED/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz'.format(subject=subject)),
            'seg': [
                join(config['deriv'][group],'subcortical_atlas/sub-{subject}/THOMAS_T1MAP_MASKED/left/1-THALAMUSfull.nii.gz'.format(subject=subject)),
                join(config['deriv'][group],'subcortical_atlas/sub-{subject}/THOMAS_T1MAP_MASKED/right/1-THALAMUSfull.nii.gz'.format(subject=subject))
                ],
            'atlas': [
                join(config['deriv'][group],'subcortical_atlas/sub-{subject}/THOMAS_T1MAP_MASKED/left/thomasfull.nii.gz'.format(subject=subject)),
                join(config['deriv'][group],'subcortical_atlas/sub-{subject}/THOMAS_T1MAP_MASKED/right/thomasrfull.nii.gz'.format(subject=subject)),
                ]
        }
    elif atlas == 'ami':
        return {
            't1w': join(config['deriv'][group],'subcortical_atlas/sub-{subject}/7TAMI/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz'.format(subject=subject)),
            'seg': join(config['deriv'][group],'subcortical_atlas/sub-{subject}/7TAMI/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz'.format(subject=subject)),
            'atlas': [
                join(config['deriv'][group],'subcortical_atlas/sub-{subject}/7TAMI/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz'.format(subject=subject)),
                join(config['deriv'][group],'subcortical_atlas/sub-{subject}/7TAMI/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz'.format(subject=subject))
                ]
        }    

# Perform ACPC alignment
rule acpc_alignment:
    input: unpack(define_thomas_input)
    output: 'deriv/presurfmorph/sub-{subject}/sub-{subject}_T1w_to_MNI152NLin2009cAsym_type-fsl.mat'
    params:
        mni = 'resources/MNI152NLin2009cAsym_t1_brain.nii.gz'
    group: 'presurfmorph'
    container: config['containers']['fsl_604']    
    shell:
        """
        template=`realpath {params.mni}`
        pushd `dirname {output}`
            robustfov -i {input.t1w} -m roi2full.mat -r robustfov.nii.gz
            convert_xfm -omat full2roi.mat -inverse roi2full.mat
            flirt -interp spline -in robustfov.nii.gz -ref $template -omat roi2std.mat -out acpc_mni.nii.gz
            convert_xfm -omat full2std.mat -concat roi2std.mat full2roi.mat
            aff2rigid full2std.mat `basename {output}`
        popd
        """

rule acpc_fsl_to_itk:
    input:
        unpack(define_thomas_input), 
        xfm = 'deriv/presurfmorph/sub-{subject}/sub-{subject}_T1w_to_MNI152NLin2009cAsym_type-fsl.mat'
    output:
        fwd = 'deriv/presurfmorph/sub-{subject}/sub-{subject}_T1w_to_MNI152NLin2009cAsym_type-itk.txt',
        inv = 'deriv/presurfmorph/sub-{subject}/MNI152NLin2009cAsym_to_sub-{subject}_T1w_type-itk.txt'
    params:
        mni = 'resources/MNI152NLin2009cAsym_t1_brain.nii.gz'
    group: 'presurfmorph' 
    container: config['containers']['surfmorph']
    shell:
        """
        c3d_affine_tool -ref {params.mni} -src {input.t1w} {input.xfm} -fsl2ras -oitk {output.fwd}
        c3d_affine_tool -ref {params.mni} -src {input.t1w} {input.xfm} -fsl2ras -inv -oitk {output.inv}
        """

rule apply_acpc_alignment:
    input:
        unpack(define_thomas_input),
        xfm = 'deriv/presurfmorph/sub-{subject}/sub-{subject}_T1w_to_MNI152NLin2009cAsym_type-itk.txt'
    output: 'deriv/presurfmorph/sub-{subject}/sub-{subject}_T1w_to_MNI152NLin2009cAsym.nii.gz',
    params:
        mni = 'resources/MNI152NLin2009cAsym_t1_brain.nii.gz'
    group: 'presurfmorph'
    # container: config['containers']['surfmorph']    
    shell:
        """
        antsApplyTransforms -d 3 -i {input.t1w} -r {params.mni} -o {output} -t {input.xfm} -n BSpline -u int
        """

# Prepare FreeSurfer data for surface-based morphometry        
rule prep_labels_surfmorph:
    input: 
        unpack(define_thomas_input),
        xfm = 'deriv/presurfmorph/sub-{subject}/sub-{subject}_T1w_to_MNI152NLin2009cAsym_type-itk.txt'
    output:
        t1w = 'deriv/surfmorph_{atlas}/labels/sub-{subject}/anat/sub-{subject}_brain_T1w.nii.gz',
        seg = expand('deriv/surfmorph_{{atlas}}/labels/sub-{{subject}}/anat/sub-{{subject}}_{hemi}_{roi}.nii.gz', hemi=hemis, roi=rois),
    params:
        rois = rois,
        structures = config['rois'][config['atlas']],
        script = "../scripts/prep_thomas.py" if config['atlas'] == 'thomas' else "../scripts/prep_surfmorph.py",
        mni = 'resources/MNI152NLin2009cAsym_t1_brain.nii.gz'
    group: 'presurfmorph'  
    script: "{params.script}"

# Flip left hemisphere
rule lr_flip_roi:
    input: 'deriv/surfmorph_{atlas}/labels/sub-{subject}/anat/sub-{subject}_left_{roi}.nii.gz'
    output: 'deriv/surfmorph_{atlas}/labels/sub-{subject}_Lflip/anat/sub-{subject}_Lflip_{roi}.nii.gz'
    container: config['containers']['surfmorph']
    group: 'presurfmorph'   
    shell:
        "c3d {input} -flip x -o {output}"

rule lr_flip_t1w:
    input: 'deriv/surfmorph_{atlas}/labels/sub-{subject}/anat/sub-{subject}_brain_T1w.nii.gz'
    output: 'deriv/surfmorph_{atlas}/labels/sub-{subject}_Lflip/anat/sub-{subject}_Lflip_brain_T1w.nii.gz'
    container: config['containers']['surfmorph']
    group: 'presurfmorph' 
    shell:
        "c3d {input} -flip x -o {output}"

# Copy right hemisphere
rule copy_right_roi:
    input: 'deriv/surfmorph_{atlas}/labels/sub-{subject}/anat/sub-{subject}_right_{roi}.nii.gz'
    output: 'deriv/surfmorph_{atlas}/labels/sub-{subject}_R/anat/sub-{subject}_R_{roi}.nii.gz'
    group: 'presurfmorph' 
    shell:
        "cp {input} {output}"

rule copy_right_t1w:
    input: 'deriv/surfmorph_{atlas}/labels/sub-{subject}/anat/sub-{subject}_brain_T1w.nii.gz'
    output: 'deriv/surfmorph_{atlas}/labels/sub-{subject}_R/anat/sub-{subject}_R_brain_T1w.nii.gz'
    group: 'presurfmorph'   
    shell:
        "cp {input} {output}"
