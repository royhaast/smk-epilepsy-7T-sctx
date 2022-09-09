# The pipeline expects the data to be in BIDS format with a file
# tar2bids.done as checkpoint for each subject.

# Gradient-distortion correction
rule gradcorrect:
    input: 'bids/sub-{subject}/tar2bids.done'
    output:
        output = 'deriv/gradcorrect/sub-{subject}/sub-{subject}_scans.tsv'
    params:
        bids_folder = 'bids', #TODO specify bids folder in config
        deriv_folder = 'deriv', #TODO specify deriv folder in config
        gradcorrect = config['containers']['gradcorrect'],
        script = 'scripts/gradcorrect/run.sh',
        coeff = config['gradientfile']
    shell:
        """
        singularity exec {params.gradcorrect} {params.script} {params.bids_folder} {params.deriv_folder}/gradcorrect participant --grad_coeff_file {params.coeff} --participant_label {wildcards.subject}
        """   

# Function to define input paths
def collect_input(wildcards,run='01'):
    subject = f'{wildcards.subject}'
    run = t1w_dict[subject]

    return {
        'inv1'    : 'deriv/gradcorrect/sub-{subject}/anat/sub-{subject}_inv-1_run-{run}_part-mag_MP2RAGE.nii.gz'.format(subject=subject,run=run),
        'inv2'    : 'deriv/gradcorrect/sub-{subject}/anat/sub-{subject}_inv-2_run-{run}_part-mag_MP2RAGE.nii.gz'.format(subject=subject,run=run),
        't1map'   : 'deriv/gradcorrect/sub-{subject}/anat/sub-{subject}_acq-MP2RAGE_run-{run}_T1map.nii.gz'.format(subject=subject,run=run),
        'uni-den' : 'deriv/gradcorrect/sub-{subject}/anat/sub-{subject}_acq-MP2RAGE_run-{run}_T1w.nii.gz'.format(subject=subject,run=run),
        'uni'     : 'deriv/gradcorrect/sub-{subject}/anat/sub-{subject}_acq-UNI_run-{run}_MP2RAGE.nii.gz'.format(subject=subject,run=run),
        'b1map'   : 'deriv/gradcorrect/sub-{subject}/fmap/sub-{subject}_part-phase_run-01_TB1TFL.nii.gz'.format(subject=subject),
        'b1map_m' : 'deriv/gradcorrect/sub-{subject}/fmap/sub-{subject}_part-mag_run-01_TB1TFL.nii.gz'.format(subject=subject)
    }

# Copy data for B1 correction
rule prepare_b1correct:
    input: unpack(collect_input)
    output:
        uni = 'deriv/presurfer/sub-{subject}/ImageUni.nii',
        inv1 = 'deriv/presurfer/sub-{subject}/ImageTI1.nii',
        inv2 = 'deriv/presurfer/sub-{subject}/ImageTI2.nii',
        b1map = 'deriv/presurfer/sub-{subject}/ImageB1mapP.nii'
    run:
        import os
        os.system('rm deriv/presurfer/sub-{}/*.nii*'.format(wildcards.subject))

        for image in output.keys():
            os.system('mri_convert {} {} -nc'.format(input[image], output[image]))
        
        for json in ['uni','b1map']:
            file_in = input[json].split('.nii.gz')[0]
            file_out = output[json].split('.nii')[0]
            os.system('cp {}.json {}.json'.format(file_in, file_out))

rule resample_b1map:
    input:
        unpack(collect_input),
        'deriv/presurfer/sub-{subject}/ImageB1mapP.nii'
    output:
        b1map_p = 'deriv/presurfer/sub-{subject}/B1reslicedImageB1mapP.nii'
    shell:
        """
        mri_convert {input.b1map} {output.b1map_p} -rl {input.inv2} -nc -rt cubic
        """

# Correct MP2RAGE data for B1+ inhomogeneities
rule b1correct:
    input: 'deriv/presurfer/sub-{subject}/B1reslicedImageB1mapP.nii'
    output:
        uni = 'deriv/presurfer/sub-{subject}/ImageUni_B1corrected.nii',
        t1map = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii',
        flaws = 'deriv/presurfer/sub-{subject}/syn_FLAWS_ImageUni_den_B1corrected.nii'
    params:
        out_folder = 'deriv/presurfer/sub-{subject}',
        wrapper = 'scripts/b1correct/b1correct.sh',
        b1correct = 'scripts/b1correct/b1correct_cemerem.m'
    shell:
        """
        bash {params.wrapper} {params.b1correct} `realpath {params.out_folder}`
        """

# Run presurfer (i.e, bias correction and brain extraction)
rule presurfer:
    input:
        uni = 'deriv/presurfer/sub-{subject}/ImageUni_B1corrected.nii',
        inv2 = 'deriv/presurfer/sub-{subject}/ImageTI2.nii'
    output:
        t1w = 'deriv/presurfer/sub-{subject}/presurf_MPRAGEise/presurf_UNI/ImageUni_B1corrected_MPRAGEised_biascorrected.nii',
	    mask = 'deriv/presurfer/sub-{subject}/presurf_MPRAGEise/presurf_UNI/ImageUni_B1corrected_MPRAGEised_brainmask.nii',
        classes = expand('deriv/presurfer/sub-{{subject}}/presurf_INV2/ImageTI2_class{i}.nii', i=range(3,7))
    params:
        wrapper = 'scripts/presurfer/presurfer.sh',
        presurf = 'scripts/presurfer/presurfer.m'
    group: 'mp2rage'        
    shell:
        """
        bash {params.wrapper} {params.presurf} {input.uni} {input.inv2}
        """

# Denoise the presurfer output using SALNM filter and apply mask
rule ants_denoise:
    input: 
        t1w = 'deriv/presurfer/sub-{subject}/presurf_MPRAGEise/presurf_UNI/ImageUni_B1corrected_MPRAGEised_biascorrected.nii',
	    mask = 'deriv/presurfer/sub-{subject}/presurf_MPRAGEise/presurf_UNI/ImageUni_B1corrected_MPRAGEised_brainmask.nii',
    output:
        denoised = 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map_T1w.nii.gz'
    group: 'mp2rage'        
    container: config['containers']['fmriprep']
    threads: 8
    resources:
        time = 30,
        mem_mb = 32000    
    shell:
        """
        DenoiseImage -d 3 -i {input.t1w} -n Rician -s 1 -o {output.denoised} -v
        """

# Compile a stripmask based on PreSurfer output (i.e., INV2 SPM segmentations)
rule generate_stripmask:
    input:
        lambda wildcards : expand('deriv/presurfer/sub-{{subject}}/presurf_INV2/ImageTI2_class{i}.nii', subject=wildcards.subject, i=range(3,7))
    output: 'deriv/presurfer/sub-{subject}/presurf_INV2/sub-{subject}_inv-2_part-mag_MP2RAGE_stripmask.nii.gz'
    container: config['containers']['fmriprep']
    group: 'mp2rage'  
    shell:
        """
        wb_command -volume-math '(1-((x1+x2+x3+x4)>0.5))' {output} -fixnan 0 -var x1 {input[0]} -var x2 {input[1]} -var x3 {input[2]} -var x4 {input[3]}
        """

# Use stripmask to mask PreSurfer processed and denoised UNI image
rule mask_t1w:
    input:
        t1w = 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map_T1w.nii.gz',
        mask = 'deriv/presurfer/sub-{subject}/presurf_INV2/sub-{subject}_inv-2_part-mag_MP2RAGE_stripmask.nii.gz'
    output: 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz'
    group: 'mp2rage'  
    container: config['containers']['fmriprep']
    shell:
        """
        fslmaths {input.t1w} -mul {input.mask} {output}
        """

# Copy FreeSurfer input file to subject's directory (i.e., 'mri/orig/001.mgz')
rule copy_t1w:
    input: 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz'
    output: 'deriv/freesurfer/sub-{subject}/mri/orig/001.mgz'
    group: 'mp2rage'        
    container: config['containers']['freesurfer']
    shell:
        """
        mri_convert {input} {output} -nc
        """

# Process the denoised image using FreeSurfer
rule freesurfer:
    input: 'deriv/freesurfer/sub-{subject}/mri/orig/001.mgz'
    output: 'deriv/freesurfer/sub-{subject}/scripts/recon-all.done'
    params:
        sd = 'deriv/freesurfer'
    group: 'mp2rage'        
    container: config['containers']['freesurfer']
    threads: 16
    resources:
        time = 600,
        mem_mb = 64000
    shell:
        """
        recon-all -all -s sub-{wildcards.subject} -hires -no-wsgcaatlas -notal-check -threads 16 -sd {params.sd}
        """

# Create ventricles mask
rule ventricles:
    input:
        freesurfer = 'deriv/freesurfer/sub-{subject}/scripts/recon-all.done',
        aseg = 'deriv/freesurfer/sub-{subject}/mri/aparc+aseg.mgz',
        t1 = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii'
    output: 'deriv/freesurfer/sub-{subject}/mri/ventricles.nii.gz'
    group: 'mp2rage'        
    container: config['containers']['freesurfer']
    shell:
        """
        mri_binarize --i {input.aseg} --ventricles --o {output}
        mri_convert {output} {output} -rl {input.t1} -nc -rt nearest
        """    

# Reslice to FreeSurfer space
rule reslice_t1:
    input:
        t1 = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii',
        orig = 'deriv/freesurfer/sub-{subject}/mri/orig.mgz'
    output: 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_space-FS_proc-B1map+PreSurfer_T1map.nii.gz'
    group: 'mp2rage'
    container: config['containers']['freesurfer']
    shell:
        """
        export FS_LICENSE=/home/rhaast/00_SOFTWARE/resources/license.txt
        mri_convert {input.t1} {output} -rl {input.orig} -nc -rt nearest
        """

# Run asegstats2table to collect segmentation statistics
rule asegstats:
    input: expand('deriv/freesurfer/sub-{subject}/stats/aseg.stats',subject=subjects),
    params:
        subjects = ['sub-{}'.format(s) for s in subjects],
        sd = 'deriv/freesurfer'
    output: 'deriv/freesurfer/aseg_volume.txt'
    group: 'mp2rage'
    container: config['containers']['fmriprep']
    shell:
        """
        export SUBJECTS_DIR={params.sd}
        asegstats2table --subjects {params.subjects} -t {output}
        """

# Transform subject T1w to 7T-AMI
rule t1w_to_7TAMI:
    input:
        template = 'resources/7TAMI_T1w_bet.nii.gz',
        t1w = 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz'
    output:
        affine = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj0GenericAffine.mat',
        warp = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1Warp.nii.gz',
        invwarp = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1InverseWarp.nii.gz'
    group: 'mp2rage'
    container: config['containers']['fmriprep']
    shell:
        """
        antsRegistrationSyN.sh -d 3 -f {input.t1w} -m {input.template} -o deriv/atlas/sub-{wildcards.subject}/7TAMI/7TAmi_T1w_bet_2_subj -n 6 -t s
        """

# Combine affine and warpfield
rule composite_transform:
    input:
        template = 'resources/7TAMI_T1w_bet.nii.gz',
        t1w = 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz',    
        affine = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj0GenericAffine.mat',
        warp = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1Warp.nii.gz'
    output: 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1CompositeWarp.nii.gz'
    group: 'mp2rage'
    container: config['containers']['fmriprep']
    shell:
        """
        antsApplyTransforms -d 3 -i {input.template} -r {input.t1w} -o [{output},1] -t {input.warp} -t [{input.affine},1]
        """    

# Get Jacobian determinant from warpfield
rule composite_jacobian_determinant:
    input: 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1CompositeWarp.nii.gz'
    output: 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1CompositeWarp_JacL.nii.gz'
    group: 'mp2rage'
    container: config['containers']['fmriprep']
    shell:
        """
        CreateJacobianDeterminantImage 3 {input} {output} 1 1
        """ 

# Apply transform to 7T-AMI atlas
rule transform_7TAMI:
    input:
        atlas = 'resources/7TAMI_DGN_SBA_v9.nii.gz',
        t1w = 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz',
        affine = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj0GenericAffine.mat',
        warp = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1Warp.nii.gz'
    output: 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_DGN_SBA_v9_to_subj.nii.gz'
    group: 'mp2rage'
    container: config['containers']['fmriprep']
    shell:
        """
        antsApplyTransforms -d 3 -i {input.atlas} -r {input.t1w} -o {output} -n MultiLabel -t {input.warp} -t {input.affine}
        """    

# Generate snapshot of WMn T1 map with THOMAS result overlaid
rule snapshot_7TAMI_seg:
    input: 
        t1map = 'deriv/subcortical_atlas/sub-{subject}/THOMAS_T1MAP/sub-{subject}_acq-MP2RAGE_desc-WMnull+denoised_T1map.nii.gz',
        atlas = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_DGN_SBA_v9_to_subj.nii.gz'
    output: 'deriv/snapshots/7TAMI/sub-{subject}_qc.png'
    group: 'mp2rage'
    shell:
        """
        fsleyes render --scene lightbox --outfile {output} --crop 5 --zaxis 2 --zrange 120 180 -ss 1.2 \
        -nc 8 {input.t1map} {input.atlas} -ot label
        """