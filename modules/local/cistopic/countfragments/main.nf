process CISTOPIC_COUNTFRAGMENTS {
    tag "Counting fragments for ${meta.id}"
    container 'docker://quay.io/cellgeni/pycistopic:2.0a0'

    input:
    tuple val(meta), path(fragments)

    output:
    tuple val(meta), path('fragment_counts.csv'), emit: csv
    path "versions.yml", emit: versions
    
    script:
    """
    zcat $fragments | calculate_fragments.awk > fragment_counts.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}' | tr -d ',')
    END_VERSIONS
    """
}