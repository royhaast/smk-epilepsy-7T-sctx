rule calculate_area:
    input: 'deriv/tmp/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-thalamus.surf.gii'
    output: 'deriv/tmp/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_label-thalamus_area.shape.gii'
    shell:
        """
        wb_command -surface-vertex-areas {input} {output}
        wb_command -metric-math 'x' {output} -fixnan 0 -var x {output}
        """

rule calculate_curvature:
    input: 'deriv/tmp/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-thalamus.surf.gii'
    output: 'deriv/tmp/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_label-thalamus_curvature.shape.gii'
    shell:
        """
        wb_command -surface-curvature {input} -mean {output}
        wb_command -metric-math 'x' {output} -fixnan 0 -var x {output}
        """

rule smooth_shape_map:
    input:
        surf = 'deriv/tmp/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-thalamus.surf.gii',
        metric = 'deriv/tmp/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_label-thalamus_{metric}.shape.gii'
    output: 'deriv/tmp/postsurfmorph_{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_label-thalamus_{metric}-smoothed.shape.gii'
    shell:
        """
        wb_command -metric-smoothing {input.surf} {input.metric} 2 {output} -fix-zeros -method GEO_GAUSS_EQUAL
        """    

rule summarize_shape_per_nuclei_ami:
    input:
        atlas = 'deriv/tmp/atlas_ami/sub-all/sub-all_label-thalamus_atlas.shape.gii',
        metrics = expand(
            'deriv/tmp/postsurfmorph_ami/sub-{{subject}}_{hemi}/sub-{{subject}}_hemi-{hemi}_label-thalamus_{metrics}-smoothed.shape.gii',
            hemi=['L','R'], metrics=['displacement','area','curvature'])
    params:
        subject_info = 'resources/subjects_displacement.csv',
        lut = 'resources/amilut.csv'
    output:
        csv = 'deriv/tmp/summary_shape_ami/sub-{subject}/sub-{subject}_label-thalamus_shape.csv',
        pkl = 'deriv/tmp/summary_shape_ami/sub-{subject}/sub-{subject}_label-thalamus_shape.pkl'
    script: '../../scripts/extract_shape_ami.py'
  
rule summarize_shape_per_nuclei_thomas:
    input:
        atlas = 'deriv/tmp/atlas_thomas/sub-all/sub-all_label-thalamus_atlas.shape.gii',
        metrics = expand(
            'deriv/tmp/postsurfmorph_thomas/sub-{{subject}}_{hemi}/sub-{{subject}}_hemi-{hemi}_label-thalamus_{metrics}-smoothed.shape.gii',
            hemi=['L','R'], metrics=['displacement','area','curvature'])
    params:
        subject_info = 'resources/subjects_displacement.csv',
        lut = 'resources/thomaslut.csv'
    output:
        csv = 'deriv/tmp/summary_shape_thomas/sub-{subject}/sub-{subject}_label-thalamus_shape.csv',
        pkl = 'deriv/tmp/summary_shape_thomas/sub-{subject}/sub-{subject}_label-thalamus_shape.pkl'
    script: '../../scripts/extract_shape_thomas.py'

rule concatenate_shape_per_nuclei:
    input: expand('deriv/tmp/summary_shape_{{atlas}}/sub-{subject}/sub-{subject}_label-thalamus_shape.pkl', subject=subjects)
    output:
        pkl = 'deriv/tmp/summary_shape_{atlas}/sub-all/sub-{group}_{atlas}_shape.pkl',
        csv = 'deriv/tmp/summary_shape_{atlas}/sub-all/sub-{group}_{atlas}_shape.csv'
    group: 'group'
    run:
        import pandas as pd

        if wildcards.group == 'controls':
            in_pkl = [ p for p in input if 'sub-C' in p ]
        elif wildcards.group == 'patients':
            in_pkl = [ p for p in input if not 'sub-C' in p ]
        
        df_concat = []
        for s in range(0,len(in_pkl)):
            df = pd.read_pickle(in_pkl[s])
            df_concat.append(df)

        df_concat = pd.concat(df_concat, ignore_index=True)

        # Save to file
        df_concat.to_pickle(output.pkl)
        df_concat.to_csv(output.csv, index=False)