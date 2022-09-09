# Generate participants file
rule gen_participants_tsv:
    input:
        t1w = expand('deriv/surfmorph/{{atlas}}/labels/sub-{subject}_{hemi}/anat/sub-{subject}_{hemi}_brain_T1w.nii.gz',subject=subjects, hemi=['Lflip','R']),
        seg = expand('deriv/surfmorph/{{atlas}}/labels/sub-{subject}_{hemi}/anat/sub-{subject}_{hemi}_{roi}.nii.gz',subject=subjects, hemi=['Lflip','R'], roi=rois)
    output: 'deriv/surfmorph/{atlas}/labels/participants.tsv'
    group: 'surfmorph_group'
    params:
        labels_dir = lambda wildcards: 'deriv/surfmorph/{atlas}/labels'.format(atlas=wildcards.atlas)
    shell:
        """
        echo participant_id > {output}
        for d in `ls -d {params.labels_dir}/*_*/` ; do
            echo `basename $d` >> {output} ; 
        done
        """

# Run the actual surface morphometry pipeline: first group level, then participant level. Repeat three times.
rule run_surfmorph_group:
    input: lambda wildcards: rules.gen_participants_tsv.output if int(wildcards.iter)-1 == 0 else 'deriv/surfmorph/{atlas}/displacement/{roi}_participant1_iter{iteration}.done'.format(iteration=int(wildcards.iter)-1, atlas=wildcards.atlas, roi=wildcards.roi)
    output: touch('deriv/surfmorph/{atlas}/displacement/{roi}_group1_iter{iter}.done')
    params:
        in_dir = lambda wildcards: directory('deriv/surfmorph/{atlas}/labels'.format(atlas=wildcards.atlas)), 
        out_dir = lambda wildcards: directory('deriv/surfmorph/{atlas}/displacement'.format(atlas=wildcards.atlas)),
        run_script = config['surfmorph_script']
    group: 'surfmorph_group'
    singularity: config['containers']['surfmorph']
    shell: 
        """
        bash {params.run_script} {params.in_dir} {params.out_dir} group1 --in_seg_dir {params.in_dir} --matching_seg {wildcards.roi} --seg_name {wildcards.roi} --matching_T1w brain --skip_bet
        """

rule run_surfmorph_participant:
    input: lambda wildcards: 'deriv/surfmorph/{atlas}/displacement/{roi}_group1_iter{iteration}.done'.format(iteration=int(wildcards.iter), atlas=wildcards.atlas, roi=wildcards.roi)
    output: 'deriv/surfmorph/{atlas}/displacement/work/{roi}_iter{iter}/sub-{subject}_{hemi}/{roi}.surf_inout.txt'
    params:
        in_dir = lambda wildcards: directory('deriv/surfmorph/{atlas}/labels'.format(atlas=wildcards.atlas)), 
        out_dir = lambda wildcards: directory('deriv/surfmorph/{atlas}/displacement'.format(atlas=wildcards.atlas)),
        run_script = config['surfmorph_script']
    group: 'surfmorph_participant'
    singularity: config['containers']['surfmorph']     
    shell: 
        """
        bash {params.run_script} {params.in_dir} {params.out_dir} participant1 --participant_label sub-{wildcards.subject}_{wildcards.hemi} --in_seg_dir {params.in_dir} --matching_seg {wildcards.roi} --seg_name {wildcards.roi} --matching_T1w brain --skip_bet
        """ 

rule surfmorph_participant_done:
    input: expand('deriv/surfmorph/{{atlas}}/displacement/work/{{roi}}_iter{{iter}}/sub-{subject}_{hemi}/{{roi}}.surf_inout.txt', subject=subjects, hemi=['Lflip','R'])
    output: touch('deriv/surfmorph/{atlas}/displacement/{roi}_participant1_iter{iter}.done')