process CISTOPIC_COUNTFRAGMENTS {
    tag "Counting fragments per barcode"
    input:
    tuple val(meta), path(fragments)

    output:
    tuple val(meta), path('fragment_counts.csv')
    
    script:
    """
    zcat $fragments | calculate_fragments.awk > fragment_counts.csv
    """
}