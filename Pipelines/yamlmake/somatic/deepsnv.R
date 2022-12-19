# Keogh only used deepSNV to exclude (ie. to improve specificity) rather than to find new variants
suppressPackageStartupMessages(library(deepSNV))

args = commandArgs(trailingOnly=TRUE)
#args <- c('processed/keogh2018/SRP159015/SRR7762430/SRR7762430_bqsr.bam',
#'processed/keogh2018/SRP159015/SRR7762431/SRR7762431_bqsr.bam',
#'processed/keogh2018/SRP159015/SRR7762430/SRR7762430_SRR7762431_somatic_varscan_JESROI.tsv')
tumour_bam = args[1]
normal_bam = args[2]
potential_variants_file = args[3]
rm(args)

regions_of_interest <- read.table(
                        file=potential_variants_file,
                        header=FALSE,
                        col.names=c('CHROM', 'POS', 'ID', 'REF', 'ALT'))

bad_rownums = which(regions_of_interest$CHR == 'chr10' & regions_of_interest$POS == '13125860')
if (length(bad_rownums) > 0) {
    warning('Omitting bad data known to cause segfault in R:\n',
        paste0(capture.output(regions_of_interest[bad_rownums,]), collapse='\n'))
    regions_of_interest = regions_of_interest[-bad_rownums,]
}

# The start and end are inclusive,
# so identical start and stop params gives us single nucleotide regions of interest
regions <- data.frame(chr   = regions_of_interest$CHR,
                      start = regions_of_interest$POS,
                      stop  = regions_of_interest$POS)

# run the actual analysis
ds <- deepSNV(
    # the corresponsing .bai files must exist too
    test    = tumour_bam,
    control = normal_bam,
    regions = regions,
    q       = 10)

# What's a good cutoff?

# what format does bcftools annotate need?
combined_df <- cbind(regions_of_interest, ds@p.val)
combined_df$QUAL <- apply(combined_df, 1, function(row) {'.'})
combined_df$FILTER <- apply(combined_df, 1,
    function(row) {
        p_value_cutoff = 0.001
        potential_alt <- row['ALT']
        if (! is.na(row[potential_alt]) && row[potential_alt] <= p_value_cutoff) {
            # should you explicitly say PASS, or should it be NULL so it
            # can be merged with pre-existing failures (if present)?
            'PASS'
        } else {
            # made-up string for now, is there an expected convention here?
            'deepSNVFail'
        }
    })
combined_df$INFO <- apply(combined_df, 1, function(row) {'.'})

vcf_cols = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
# do this right
# quick and dirty vcf
# we want no quotes and to start with a #
write(vcf_cols, sep=)
# Also need to have:
#
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">          
##FILTER=<ID=deepSNVFail,Description="something something">  

# And then it actually works!!!

write.table(
    combined_df[,vcf_cols],
    sep='\t',
    col.names=FALSE
    row.names=FALSE)
