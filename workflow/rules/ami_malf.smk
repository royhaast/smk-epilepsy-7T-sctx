# This set of rules requires the 7T-AMI MALF pipeline to be ran a priori.
# Currently this pipeline is not yet publicly available but hopefully soon.
# Instead single-step coregistration results will be used here.

# Make sure the 7T-AMI MALF results are in same space as the initial MP2RAGE data
rule resample_ami_malf:
    input:
        ami = 'deriv/tmp/ami_malf/sub-{subject}/sub-{subject}_7TAMI60qfullJF_Labels.nii.gz',
        t1 = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii'
    output:
        ami = 'deriv/tmp/summary_ami_malf/sub-{subject}/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
    shell:
        """
        mri_convert {input.ami} {output} -rl {input.t1} -rt nearest
        """

# Prepare results for surfmorph analyses
rule create_symlinks_for_surfmorph:
#TODO edit to create symlinks
    input:
        ami = 'deriv/tmp/summary_ami_malf/sub-{subject}/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
        t1w = 'deriv/presurfer/sub-{subject}/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz',
        mask = 'deriv/presurfer/sub-{subject}/presurf_INV2/sub-{subject}_inv-2_part-mag_MP2RAGE_stripmask.nii.gz'
    output:
        ami = 'deriv/tmp/copy_for_shape/sub-{subject}/7TAMI/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
        t1w = 'deriv/tmp/copy_for_shape/sub-{subject}/7TAMI/sub-{subject}_acq-MP2RAGE_proc-B1map+PreSurfer_T1w.nii.gz'
    shell:
        """
        fslmaths {input.mask} -bin -mul {input.ami} {output.ami}
        cp {input.t1w} {output.t1w}
        """

# Extract statistics from 7T-AMI MALF thalamic segmentations
rule thalamic_ami_malf_stats:
    input:
        stats = 'deriv/freesurfer/aseg_volume.txt',
        thalamus = 'deriv/tmp/summary_ami_malf/sub-{subject}/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
        t1 = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii',
        mask = 'deriv/presurfer/sub-{subject}/presurf_INV2/sub-{subject}_inv-2_part-mag_MP2RAGE_stripmask.nii.gz'
    params:
        dbm = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1CompositeWarp_JacL.nii.gz',
        r2s = 'deriv/qsm/sub-{subject}/sub-{subject}_space-MP2RAGE_proc-GDC_R2s.nii.gz',
        qsm = 'deriv/qsm/sub-{subject}/sub-{subject}_space-MP2RAGE_proc-GDC_QSM.nii.gz',
        fa  = 'deriv/tmp/dwi/sub-{subject}/metrics/sub-{subject}_space-anat_fa.nii.gz',
        adc = 'deriv/tmp/dwi/sub-{subject}/metrics/sub-{subject}_space-anat_adc.nii.gz',  
        t1_threshold = 3000,  
        subject_info = 'resources/bids_decoder.csv',
        lut = 'resources/amilut.csv'
    output:
        pkl = 'deriv/tmp/summary_ami_malf/sub-{subject}/sub-{subject}_ami_summary.pkl',
        csv = 'deriv/tmp/summary_ami_malf/sub-{subject}/sub-{subject}_ami_summary.csv'
    group: 'subj'
    script: '../../scripts/extract_thalamus_ami.py'

# Extract statistics from 7T-AMI non-thalamic segmentations
rule subcortical_ami_malf_stats:
    input:
        stats = 'deriv/freesurfer/aseg_volume.txt',
        subcortex = 'deriv/tmp/summary_ami_malf/sub-{subject}/sub-{subject}_space-native_7TAMI60_qfullJF_Labels.nii.gz',
        amygdala = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_Amygdala_to_subj.nii.gz',
        t1 = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii',
        mask = 'deriv/presurfer/sub-{subject}/presurf_INV2/sub-{subject}_inv-2_part-mag_MP2RAGE_stripmask.nii.gz'
    params:
        dbm = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1CompositeWarp_JacL.nii.gz',
        r2s = 'deriv/qsm/sub-{subject}/sub-{subject}_space-MP2RAGE_proc-GDC_R2s.nii.gz',
        qsm = 'deriv/qsm/sub-{subject}/sub-{subject}_space-MP2RAGE_proc-GDC_QSM.nii.gz',
        fa  = 'deriv/tmp/dwi/sub-{subject}/metrics/sub-{subject}_space-anat_fa.nii.gz',
        adc = 'deriv/tmp/dwi/sub-{subject}/metrics/sub-{subject}_space-anat_adc.nii.gz',   
        t1_threshold = 3000,   
        subject_info = 'resources/bids_decoder.csv',
        lut = 'resources/amilut_subcortex.csv'
    output:
        pkl = 'deriv/tmp/summary_ami_malf/sub-{subject}/sub-{subject}_subcortex_summary.pkl',
        csv = 'deriv/tmp/summary_ami_malf/sub-{subject}/sub-{subject}_subcortex_summary.csv'
    group: 'subj'
    script: '../../scripts/extract_subcortex.py'

# Concatenate stats across subjects
rule concatenate_ami_malf_qmri:
    input: expand('deriv/tmp/summary_ami_malf/sub-{subject}/sub-{subject}_{{atlas}}_summary.pkl', subject=subjects)
    output:
        pkl = 'deriv/tmp/summary_ami_malf/sub-patients/sub-patients_{atlas}_summary.pkl',
        csv = 'deriv/tmp/summary_ami_malf/sub-patients/sub-patients_{atlas}_summary.csv'
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