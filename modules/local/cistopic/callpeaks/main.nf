process CISTOPIC_CALLPEAKS {
    tag "Performing peak calling for ${celltype_name}"
    input:
    tuple val(meta), path(fragments)

    output:
    tuple val(meta), path("$task.ext.outputdir/*.narrowPeak"), env(large_peaks_num), env(all_peaks_num), emit: narrowPeak
    path 'versions.yml', emit: versions
    
    script:
    def args = task.ext.args ?: ''
    """
    # perform peak calling
    peak_calling.py \\
            $args \\
            --bed_path $fragments \\
            --output_dir "${task.ext.outputdir}" \\
            --cpus $task.cpus \\
            --skip_empty_peaks
    
    # count peaks
    large_peaks_num=\$(cut -f4 ${task.ext.outputdir}/*.narrowPeak | sed 's/[a-z]\$//' | sort -u | wc -l)
    all_peaks_num=\$(wc -l < ${task.ext.outputdir}/*.narrowPeak)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pycisTopic: \$( python -c "import pycisTopic; print(pycisTopic.__version__)" )
    END_VERSIONS
    """
}