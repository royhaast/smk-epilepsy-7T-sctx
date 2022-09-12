# Convert subject space surface to vtk
rule subject_space_to_vtk:
   input: 'deriv/postsurfmorph/{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-{roi}.surf.gii'
   output: 'deriv/postsurfmorph/{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-{roi}.vtk'
   group: 'atlas_participant'        
   shell:
       """
       mris_convert {input} {output}
       """ 

rule left_atlas_to_nifti:
    input: unpack(define_input) 
    output: 'deriv/atlas/{atlas}/sub-{subject}/sub-{subject}_hemi-L_atlas.nii.gz'
    group: 'atlas_participant'
    shell:
        """
        mri_convert {input.atlas[0]} {output} -nc -rt nearest
        """  

rule right_atlas_to_nifti:
    input: unpack(define_input)
    output: 'deriv/atlas/{atlas}/sub-{subject}/sub-{subject}_hemi-R_atlas.nii.gz'
    group: 'atlas_participant'
    shell:
        """
        mri_convert {input.atlas[1]} {output} -nc -rt nearest
        """

rule atlas_dilate:
    input: 'deriv/atlas/{atlas}/sub-{subject}/sub-{subject}_hemi-{hemi}_atlas.nii.gz'
    output: 'deriv/atlas/{atlas}/sub-{subject}/sub-{subject}_hemi-{hemi}_atlas_dilated-5mm.nii.gz'
    group: 'atlas_participant'    
    container: config['containers']['surfmorph']    
    shell:
        """
        fslmaths {input} -kernel boxv 5 -dilD {output}
        """

rule atlas_to_surf:
    input: 
        surf = 'deriv/postsurfmorph/{atlas}/sub-{subject}_{hemi}/sub-{subject}_hemi-{hemi}_space-T1w_label-{roi}.vtk',
        atlas = 'deriv/atlas/{atlas}/sub-{subject}/sub-{subject}_hemi-{hemi}_atlas_dilated-5mm.nii.gz'
    output: 'deriv/atlas/{atlas}/sub-{subject}/sub-{subject}_hemi-{hemi}_label-{roi}_atlas.shape.gii'
    group: 'atlas_participant'    
    script: 'scripts/atlas_to_gifti.py'

rule gen_atlas_gii:
    input: expand('deriv/atlas/{atlas}/sub-{subject}/sub-{subject}_hemi-{hemi}_label-{{roi}}_atlas.shape.gii', atlas=config['atlas'], subject=subjects, hemi=['R'])
    output: 'deriv/atlas/{atlas}/sub-all/sub-all_label-{roi}_atlas.shape.gii'
    params:
        labels = lambda wildcards: config['rois'][config['atlas']][wildcards.roi]['right'],
    group: 'atlas_group'    
    script: 'scripts/gen_atlas.py'    

rule get_atlas_label:
    input:
        gii = 'deriv/atlas/{atlas}/sub-all/sub-all_label-{roi}_atlas.shape.gii',
    output: 'deriv/atlas/{atlas}/sub-all/sub-all_label-{roi}_atlas.label.gii'
    params:
        lut = lambda wildcards: config['lut'][wildcards.atlas]
    group: 'atlas_group'
    shell:
        """
        wb_command -metric-label-import {input.gii} {params.lut} {output}
        """

rule get_atlas_border:
    input:
        surf = expand('deriv/postsurfmorph/{atlas}/sub-{subject}_R/sub-{subject}_hemi-R_space-T1w_label-{{roi}}.surf.gii', atlas=config['atlas'], subject=subjects[0]),
        label = 'deriv/atlas/{atlas}/sub-all/sub-all_label-{roi}_atlas.label.gii',
    output: 'deriv/atlas/{atlas}/sub-all/sub-all_label-{roi}_atlas.border'
    group: 'atlas_group'
    shell:
        """
        wb_command -label-to-border {input.surf} {input.label} {output}
        """        

# rule summarize_shape_per_nuclei_ami:
#     input:
#         atlas = 'deriv/atlas/ami/sub-all/sub-all_label-{roi}_atlas.label.gii',
#         metrics = expand(
#             'deriv/postsurfmorph/ami/sub-{{subject}}_{hemi}/sub-{{subject}}_hemi-{hemi}_label-{{roi}}_{metrics}-smoothed.shape.gii',
#             hemi=['L','R'], metrics=['displacement','area','curvature'])
#     params:
#         subject_info = 'resources/subjects_displacement.csv',
#         lut = lambda wildcards: 'resources/amilut.csv' if wildcards.roi == 'thalamus' else 'resources/amilut_subcortex.csv',
#         labels = lambda wildcards: config['rois'][config['atlas']][wildcards.roi]['right']
#     group: 'atlas_participant'
#     output:
#         csv = 'deriv/summary_shape_ami/sub-{subject}/sub-{subject}_label-{roi}_shape.csv',
#         pkl = 'deriv/summary_shape_ami/sub-{subject}/sub-{subject}_label-{roi}_shape.pkl'
#     script: '../scripts/extract_shape_ami.py'
  
# rule summarize_shape_per_nuclei_thomas:
#     input:
#         atlas = 'deriv/atlas/thomas/sub-all/sub-all_label-thalamus_atlas.shape.gii',
#         metrics = expand(
#             'deriv/postsurfmorph/thomas/sub-{{subject}}_{hemi}/sub-{{subject}}_hemi-{hemi}_label-thalamus_{metrics}-smoothed.shape.gii',
#             hemi=['L','R'], metrics=['displacement','area','curvature'])
#     params:
#         subject_info = 'resources/subjects_displacement.csv',
#         lut = 'resources/thomaslut.csv'
#     group: 'atlas_participant'        
#     output:
#         csv = 'deriv/summary_shape_thomas/sub-{subject}/sub-{subject}_label-thalamus_shape.csv',
#         pkl = 'deriv/summary_shape_thomas/sub-{subject}/sub-{subject}_label-thalamus_shape.pkl'
#     script: '../scripts/extract_shape_thomas.py'

# rule concatenate_shape_per_nuclei:
#     input: expand('deriv/summary_shape_{{atlas}}/sub-{subject}/sub-{subject}_label-thalamus_shape.pkl', subject=subjects_displacement)
#     output:
#         pkl = 'deriv/summary_shape_{atlas}/sub-all/sub-{group}_atlas-{atlas}_label-thalamus_shape.pkl',
#         csv = 'deriv/summary_shape_{atlas}/sub-all/sub-{group}_atlas-{atlas}_label-thalamus_shape.csv'
#     group: 'group'
#     run:
#         import pandas as pd

#         if wildcards.group == 'controls':
#             in_pkl = [ p for p in input if 'sub-C' in p ]
#         elif wildcards.group == 'patients':
#             in_pkl = [ p for p in input if not 'sub-C' in p ]
        
#         df_concat = []
#         for s in range(0,len(in_pkl)):
#             df = pd.read_pickle(in_pkl[s])
#             df_concat.append(df)

#         df_concat = pd.concat(df_concat, ignore_index=True)

#         # Save to file
#         df_concat.to_pickle(output.pkl)
#         df_concat.to_csv(output.csv, index=False)

# rule concatenate_shape_for_subcortex:
#     input: expand('deriv/summary_shape_ami/sub-{subject}/sub-{subject}_label-{roi}_shape.pkl', subject=subjects_displacement, roi=rois)
#     output:
#         pkl = 'deriv/summary_shape_ami/sub-all/sub-{group}_atlas-subcortex_shape.pkl',
#         csv = 'deriv/summary_shape_ami/sub-all/sub-{group}_atlas-subcortex_shape.csv'
#     group: 'group'
#     run:
#         import pandas as pd

#         if wildcards.group == 'controls':
#             in_pkl = [ p for p in input if 'sub-C' in p ]
#         elif wildcards.group == 'patients':
#             in_pkl = [ p for p in input if not 'sub-C' in p ]
        
#         df_concat = []
#         for s in range(0,len(in_pkl)):
#             df = pd.read_pickle(in_pkl[s])
#             df_concat.append(df)

#         df_concat = pd.concat(df_concat, ignore_index=True)

#         # Save to file
#         df_concat.to_pickle(output.pkl)
#         df_concat.to_csv(output.csv, index=False)