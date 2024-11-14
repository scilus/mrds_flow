#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["use_isotropic":"$params.use_isotropic",
                "r_threshold":"$params.r_threshold",
                "mrds_processes":"$params.mrds_processes",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "SCIL MRDS pipeline"
log.info "=========================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[MRDS options]"
log.info "Use isotropic Compartment: $params.use_isotropic"
log.info "Number of processes used: $params.processes"
log.info "Number of MRDS processes: $params.mrds_processes"
log.info ""
log.info ""
log.info ""

log.info "Input: $params.input"
root = file(params.input)

if(params.mrds_processes < 0) {
    error "Error params.mrds_processes should be higher than 0"
}

if(params.mrds_processes > params.processes) {
    error "Error params.mrds_processes should be lower than the params.processes - Currently ${params.mrds_processes} > ${params.processes}"
}

/* Watch out, files are ordered alphabetically in channel */
Channel
    .fromFilePairs("$root/**/*{bval,bvec,dwi.nii.gz}",
                    size: 3,
                    maxDepth:1,
                    flat: true) {it.parent.name}
    .set{data}

(dwi, only_dwi, gradients_to_be_converted, count_subjects) = data
    .map{sid, bvals, bvecs, dwi -> [tuple(sid, bvals, bvecs, dwi),
                                    tuple(sid, dwi),
                                    tuple(sid, bvals, bvecs),
                                    sid]}
    .separate(4)

only_dwi.into{only_dwi_for_mrds; only_dwi_for_todi; only_dwi_for_modsel}

Channel
    .fromFilePairs("$root/**/{*tracking*.*,}",
                    size: -1, maxDepth:1) {it.parent.name}
    .into{count_tractograms; tractogram_for_todi; tractogram_for_mask}

Channel
    .fromFilePairs("$params.input/**/*mask.nii.gz",
        size: -1) { it.parent.name }
    .into{count_masks; create_masks; mask_for_mrds; mask_for_modsel; mask_for_metrics}

count_subjects.count()
    .concat(count_tractograms.count())
    .concat(count_masks.count())
    .toList()
    .subscribe{a, b, c-> if (a != b || a != c || b != c) {
        error "Error ~ The number of subjects, tractograms and masks should be the same."
    }}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

process Convert_Scheme {
    input:
    set sid, path(bval), path(bvec) from gradients_to_be_converted

    output:
    set sid, "${sid}__scheme.b" into scheme_for_mrds

    script:
    """
    scil_gradients_convert.py ${bval} ${bvec} ${sid}__scheme.b -f --input_fsl
    """
}

process Compute_Mask {
    input:
    set sid, path(tracking) from tractogram_for_mask

    output:
    set sid, "${sid}__mask.nii.gz" into computed_mask_for_mrds

    when:
    (create_masks.count() == 0)

    script:
    """
    scil_tractogram_compute_density_map.py ${tracking} ${sid}__mask.nii.gz --binary -f
    """
}

only_dwi_for_mrds
    .combine(scheme_for_mrds, by:0)
    .combine(computed_mask_for_mrds
          .mix(mask_for_mrds), by:0)
    .set{dwi_scheme_mask_for_mrds}

process Fit_MRDS {
    input:
    set sid, path(dwi), path(scheme), path(mask) from dwi_scheme_mask_for_mrds

    output:
    set sid, "${sid}__D1_signal_fraction.nii.gz", \
             "${sid}__D1_evals.nii.gz", \
             "${sid}__D1_isotropic.nii.gz", \
             "${sid}__D1_num_tensors.nii.gz", \
             "${sid}__D1_evecs.nii.gz", \
             "${sid}__D2_signal_fraction.nii.gz", \
             "${sid}__D2_evals.nii.gz", \
             "${sid}__D2_isotropic.nii.gz", \
             "${sid}__D2_num_tensors.nii.gz", \
             "${sid}__D2_evecs.nii.gz", \
             "${sid}__D3_signal_fraction.nii.gz", \
             "${sid}__D3_evals.nii.gz", \
             "${sid}__D3_isotropic.nii.gz", \
             "${sid}__D3_num_tensors.nii.gz", \
             "${sid}__D3_evecs.nii.gz" into mrds_for_modsel optional true
    
    path("${sid}__DTInolin_ResponseAnisotropic.txt") optional true
    path("${sid}_d_perp.txt") optional true
    path("${sid}_d_par.txt") optional true
    path("*nii.gz") optional true

    script:
    iso=""
    if (params.use_isotropic)
        iso="-iso"
    """
    export OMP_NUM_THREADS=$params.mrds_processes
    
    dti $dwi $scheme $sid -mask $mask -response 0 -correction 0
    rm -rf *_DTInolin*.nii.gz *_DTInolin_ResponseIsotropic.txt
    cp ${sid}_DTInolin_ResponseAnisotropic.txt DTInolin_ResponseAnisotropic.txt

    # Extract diffisivities
    line=\$(head -n 1 DTInolin_ResponseAnisotropic.txt)
    values=(\${line// / })
    d_par=\${values[0]}
    d_perp=\$(echo "\${values[1]}+\${values[2]}/2" | bc -l | awk '{printf "%f", \$0}')

    echo \${d_perp} > ${sid}_d_perp.txt
    echo \${d_par} >  ${sid}_d_par.txt

    mdtmrds $dwi \
        $scheme \
        $sid \
        -correction 0 \
        -response \${d_par},\${d_perp},0.003 \
        -mask $mask \
        -modsel bic \
        -each -intermediate \
        $iso \
        -mse \
        -method Diff
    
    for v in V1 V2 V3; do
        mv ${sid}__MRDS_Diff_\${v}_COMP_SIZE.nii.gz ${sid}__\${v/V/D}_signal_fraction.nii.gz
        mv ${sid}__MRDS_Diff_\${v}_PDDs_CARTESIAN.nii.gz ${sid}__\${v/V/D}_evecs.nii.gz
        mv ${sid}__MRDS_Diff_\${v}_COMP_SIZE.nii.gz ${sid}__\${v/V/D}_signal_fraction.nii.gz
        mv ${sid}__MRDS_Diff_\${v}_EIGENVALUES.nii.gz ${sid}__\${v/V/D}_evals.nii.gz
        mv ${sid}__MRDS_Diff_\${v}_ISOTROPIC.nii.gz ${sid}__\${v/V/D}_isotropic.nii.gz
        mv ${sid}__MRDS_Diff_\${v}_NUM_COMP.nii.gz ${sid}__\${v/V/D}_num_tensors.nii.gz
    done
    """
}

only_dwi_for_todi
    .combine(tractogram_for_todi, by: 0)
    .set{dwi_tractogram_for_todi}

process Compute_TODI {
    input:
    set sid, path(dwi), path(tractogram) from dwi_tractogram_for_todi

    output:
    set sid, "${sid}__NUFO.nii.gz" into nufo_for_modsel
    path("${sid}__TODI_SH.nii.gz")

    script:
    """
    scil_tractogram_compute_TODI.py ${tractogram} \
        --out_todi_sh ${sid}__TODI_SH.nii.gz \
        --reference ${dwi} \
        --sh_basis descoteaux07 -f

    scil_fodf_metrics.py ${sid}__TODI_SH.nii.gz \
        --nufo ${sid}__NUFO.nii.gz \
        --not_all \
        --sh_basis descoteaux07 \
        --rt ${params.r_threshold} -f
    """
}

nufo_for_modsel
    .combine(only_dwi_for_modsel, by: 0)
    .combine(mask_for_modsel, by: 0)
    .combine(mrds_for_modsel, by: 0)
    .set{dwi_nufo_mrds_for_modsel}


process Modsel {
    input:
    set sid, path(nufo), path(dwi), path(mask), \
        path(n1_signal_fraction), path(n1_eigen), path(n1_iso), path(n1_num_tensors), path(n1_evecs), \
        path(n2_signal_fraction), path(n2_eigen), path(n2_iso), path(n2_num_tensors), path(n2_evecs), \
        path(n3_signal_fraction), path(n3_eigen), path(n3_iso), path(n3_num_tensors), path(n3_evecs) from dwi_nufo_mrds_for_modsel

    output:
    set sid, "${sid}__evals.nii.gz" into eigenvalues_for_metrics
    path("${sid}__signal_fraction.nii.gz")
    path("${sid}__isotropic.nii.gz")
    path("${sid}__num_tensors.nii.gz")
    path("${sid}__evecs.nii.gz")

    script:
    """
    scil_mrds_select_number_of_tensors.py ${nufo} ${dwi} \
        --N1 ${n1_signal_fraction} ${n1_eigen} ${n1_iso} ${n1_num_tensors} ${n1_evecs} \
        --N2 ${n2_signal_fraction} ${n2_eigen} ${n2_iso} ${n2_num_tensors} ${n2_evecs} \
        --N3 ${n3_signal_fraction} ${n3_eigen} ${n3_iso} ${n3_num_tensors} ${n3_evecs} \
        --prefix ${sid}_ \
        --mask ${mask}
    
    mv ${sid}__PDDs_CARTESIAN.nii.gz ${sid}__evecs.nii.gz
    mv ${sid}__COMP_SIZE.nii.gz ${sid}__signal_fraction.nii.gz
    mv ${sid}__EIGENVALUES.nii.gz ${sid}__evals.nii.gz
    mv ${sid}__ISOTROPIC.nii.gz ${sid}__isotropic.nii.gz
    mv ${sid}__NUM_COMP.nii.gz ${sid}__num_tensors.nii.gz
    """
}

eigenvalues_for_metrics
    .combine(mask_for_metrics, by: 0)
    .set{eigenvalues_mask_for_metrics}

process MRDS_Metrics {
    input:
    set sid, path(eigenvalues), path(mask) from eigenvalues_mask_for_metrics

    output:
    path("${sid}__MRDS_AD.nii.gz")
    path("${sid}__MRDS_RD.nii.gz")
    path("${sid}__MRDS_MD.nii.gz")
    path("${sid}__MRDS_FA.nii.gz")

    script:
    """
    scil_mrds_metrics.py ${eigenvalues} \
        --not_all \
        --fa ${sid}__MRDS_FA.nii.gz \
        --ad ${sid}__MRDS_AD.nii.gz \
        --rd ${sid}__MRDS_RD.nii.gz \
        --md ${sid}__MRDS_MD.nii.gz \
        --mask ${mask} 
    """
}