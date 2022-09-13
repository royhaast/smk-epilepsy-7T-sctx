# This set of rules requires the 7T-AMI MALF pipeline to be ran a priori.
# Currently this pipeline is not yet publicly available but hopefully soon.
# Instead single-step coregistration results will be used here.

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
        antsRegistrationSyN.sh -d 3 -f {input.t1w} -m {input.template} -o deriv/atlas/sub-{wildcards.subject}/7TAMI/7TAmi_T1w_bet_2_subj -n 8 -t s
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

# Make sure the 7T-AMI MALF results are in same space as the initial MP2RAGE data
rule resample_ami_malf:
    input:
        ami = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_DGN_SBA_v9_to_subj.nii.gz',
        t1 = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii'
    output:
        ami = 'deriv/summary/sub-{subject}/7TAMI/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
    shell:
        """
        mri_convert {input.ami} {output} -rl {input.t1} -rt nearest
        """

# Generate snapshot of WMn T1 map with THOMAS result overlaid
rule snapshot_ami_malf_seg:
    input: 
        t1w = 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz',
        atlas = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_DGN_SBA_v9_to_subj.nii.gz'
    output: 'deriv/snapshots/7TAMI/sub-{subject}_qc.png'
    group: 'mp2rage'
    shell:
        """
        fsleyes render --scene lightbox --outfile {output} --crop 5 --zaxis 2 --zrange 120 180 -ss 1.2 \
        -nc 8 {input.t1w} {input.atlas} -ot label
        """

# Extract statistics from 7T-AMI MALF thalamic segmentations
rule thalamic_ami_malf_stats:
    input:
        stats = 'deriv/freesurfer/aseg_volume.txt',
        thalamus = 'deriv/summary/sub-{subject}/7TAMI/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
        t1 = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii',
        dbm = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1CompositeWarp_JacL.nii.gz',
        mask = 'deriv/presurfer/sub-{subject}/presurf_INV2/sub-{subject}_inv-2_part-mag_MP2RAGE_stripmask.nii.gz'
    params:
        t1_threshold = 3000,  
        subject_info = 'resources/participants.tsv',
        lut = 'resources/amilut.csv'
    output:
        pkl = 'deriv/summary/sub-{subject}/7TAMI/sub-{subject}_ami_summary.pkl',
        csv = 'deriv/summary/sub-{subject}/7TAMI/sub-{subject}_ami_summary.csv'
    group: 'subj'
    script: '../scripts/python/extract_thalamus_ami.py'

# Extract statistics from 7T-AMI non-thalamic segmentations
rule subcortical_ami_malf_stats:
    input:
        stats = 'deriv/freesurfer/aseg_volume.txt',
        subcortex = 'deriv/summary/sub-{subject}/7TAMI/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
        t1 = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii',
        dbm = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1CompositeWarp_JacL.nii.gz', 
        mask = 'deriv/presurfer/sub-{subject}/presurf_INV2/sub-{subject}_inv-2_part-mag_MP2RAGE_stripmask.nii.gz'
    params:
        t1_threshold = 3000,   
        subject_info = 'resources/participants.tsv',
        lut = 'resources/amilut_subcortex.csv'
    output:
        pkl = 'deriv/summary/sub-{subject}/7TAMI/sub-{subject}_subcortex_summary.pkl',
        csv = 'deriv/summary/sub-{subject}/7TAMI/sub-{subject}_subcortex_summary.csv'
    group: 'subj'
    script: '../scripts/python/extract_subcortex.py'

# Concatenate stats across subjects
rule concatenate_ami_malf_stats:
    input: expand('deriv/summary/sub-{subject}/7TAMI/sub-{subject}_{{atlas}}_summary.pkl', subject=subjects)
    output:
        pkl = 'deriv/summary/sub-all/7TAMI/sub-all_{atlas}_summary.pkl',
        csv = 'deriv/summary/sub-all/7TAMI/sub-all_{atlas}_summary.csv'
    group: 'group'
    run:
        import pandas as pd

        df_concat = []
        for s in range(0,len(input)):
            df = pd.read_pickle(input[s])
            df_concat.append(df)

        df_concat = pd.concat(df_concat, ignore_index=True)

        # Save to file
        df_concat.to_pickle(output.pkl)
        df_concat.to_csv(output.csv)

# # Prepare results for surfmorph analyses
# rule create_symlinks_for_surfmorph:
# #TODO edit to create symlinks
#     input:
#         ami = 'deriv/summary/7TAMI/sub-{subject}/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
#         t1w = 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz',
#         mask = 'deriv/presurfer/sub-{subject}/presurf_INV2/sub-{subject}_inv-2_part-mag_MP2RAGE_stripmask.nii.gz'
#     output:
#         ami = 'deriv/copy_for_shape/sub-{subject}/7TAMI/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
#         t1w = 'deriv/copy_for_shape/sub-{subject}/7TAMI/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz'
#     shell:
#         """
#         fslmaths {input.mask} -bin -mul {input.ami} {output.ami}
#         cp {input.t1w} {output.t1w}
#         """        