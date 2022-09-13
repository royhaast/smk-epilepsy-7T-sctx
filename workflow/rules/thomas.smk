# Generate WM nulled T1 image
rule null_t1:
    input:
        t1map = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii',
        mask = 'deriv/presurfer/sub-{subject}/presurf_MPRAGEise/presurf_UNI/ImageUni_B1corrected_MPRAGEised_brainmask.nii'
    output: 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/sub-{subject}_acq-MP2RAGE_desc-WMnull_T1map.nii.gz'
    params:
        threshold = 700
    group: 'mp2rage'
    shell:
        """
        fslmaths {input.t1map} -recip -mul -{params.threshold} -exp -mul -2.0 -add 1.0 -abs -mul {input.mask} -thr 0 -uthr 2 {output}
        """

# Denoise WM nulled image prior to THOMAS
rule null_t1_denoise:
    input: 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/sub-{subject}_acq-MP2RAGE_desc-WMnull_T1map.nii.gz'
    output: 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/sub-{subject}_acq-MP2RAGE_desc-WMnull+denoised_T1map.nii.gz'
    group: 'mp2rage'
    shell:
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8
        DenoiseImage -d 3 -i {input} -o {output} -n Rician -v
        """

# Segment WM nulled T1 map with THOMAS
rule null_t1_thomas:
    input:
        t1map = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/sub-{subject}_acq-MP2RAGE_desc-WMnull+denoised_T1map.nii.gz'
    output:
        thalamus_left = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/left/1-THALAMUS.nii.gz',
        thalamus_right = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/right/1-THALAMUS.nii.gz',
        nuclei_left = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/left/thomasfull.nii.gz',
        nuclei_right = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/right/thomasrfull.nii.gz'
    params:
        thomas_pipeline = config['thomas']
    group: 'mp2rage'
    shell:
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8
        export THOMAS_HOME={params.thomas_pipeline}
        pushd deriv/atlas/sub-{wildcards.subject}/THOMAS_T1MAP
        tcsh $THOMAS_HOME/thomas_csh `basename {input.t1map}`
        popd
        """

# Reslice whole thalamus label to native space
rule null_t1_thomas_reslice:
    input:       
        thalamus = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/{hemi}/1-THALAMUS.nii.gz',
        t1map = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/sub-{subject}_acq-MP2RAGE_desc-WMnull_T1map.nii.gz'
    output: 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/{hemi}/1-THALAMUSfull.nii.gz'
    group: 'mp2rage'
    shell:
        """
        mri_convert {input.thalamus} {output} -rl {input.t1map} -rt nearest -nc
        """

# Extract statistics from THOMAS thalamic segmentations
rule thalamus_thomas_stats:
    input:
        stats = 'deriv/freesurfer/aseg_volume.txt',
        thomas = ['deriv/atlas/sub-{subject}/THOMAS_T1MAP/left/thomasfull.nii.gz',
                    'deriv/atlas/sub-{subject}/THOMAS_T1MAP/right/thomasrfull.nii.gz'],
        thalamus = ['deriv/atlas/sub-{subject}/THOMAS_T1MAP/left/1-THALAMUSfull.nii.gz',
                    'deriv/atlas/sub-{subject}/THOMAS_T1MAP/right/1-THALAMUSfull.nii.gz'],                    
        t1 = 'deriv/presurfer/sub-{subject}/T1map_B1corrected.nii',
        dbm = 'deriv/atlas/sub-{subject}/7TAMI/7TAmi_T1w_bet_2_subj1CompositeWarp_JacL.nii.gz',
        mask = 'deriv/presurfer/sub-{subject}/presurf_INV2/sub-{subject}_inv-2_part-mag_MP2RAGE_stripmask.nii.gz'
    params:
        t1_threshold = 3000,         
        subject_info = 'resources/participants.tsv',
        lut = 'resources/thomaslut.csv'
    output:
        pkl = 'deriv/summary/sub-{subject}/THOMAS_T1MAP/sub-{subject}_thomas_summary.pkl',
        csv = 'deriv/summary/sub-{subject}/THOMAS_T1MAP/sub-{subject}_thomas_summary.csv'
    group: 'subj'
    script: '../scripts/python/extract_thalamus_thomas.py'

# Concatenate stats across subjects
rule concatenate_thomas_stats:
    input: expand('deriv/summary/sub-{subject}/THOMAS_T1MAP/sub-{subject}_{{atlas}}_summary.pkl', subject=subjects)
    output:
        pkl = 'deriv/summary/sub-all/THOMAS_T1MAP/sub-all_{atlas}_summary.pkl',
        csv = 'deriv/summary/sub-all/THOMAS_T1MAP/sub-all_{atlas}_summary.csv'
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

# Will generate snapshot of WMn T1 map with THOMAS result overlaid
rule snapshot_thomas_seg:
    input: 
        t1map = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/sub-{subject}_acq-MP2RAGE_desc-WMnull+denoised_T1map.nii.gz',
        left = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/left/1-THALAMUS.nii.gz',
        right = 'deriv/atlas/sub-{subject}/THOMAS_T1MAP/right/1-THALAMUS.nii.gz',
    output: 'deriv/snapshots/THOMAS_T1MAP/sub-{subject}_qc.png'
    group: 'mp2rage'
    shell:
        """
        fsleyes render --scene lightbox --outfile {output} --crop 5 --zaxis 2 --zrange 120 180 -ss 1.2 \
        -nc 8 {input.t1map} {input.left} -ot label {input.right} -ot label
        """
