#!/usr/bin/awk -f

# Skip comment lines that start with #
/^#/ { next }

# Process data lines
NF >= 5 {
    # Extract barcode (4th column) and count (5th column)
    barcode = $4
    count = $5
    
    # Sum counts for each barcode
    barcode_counts[barcode] += count
}

END {
    # Print results
    print "barcode\tTotal_Fragments"
    for (barcode in barcode_counts) {
        print barcode "," barcode_counts[barcode]
    }
}