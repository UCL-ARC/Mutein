# Keogh only used deepSNV to exclude (ie. to improve specificity) rather than to find new variants
library("optparse")
library("tidyverse")

option_list = list(
  make_option(c("--tumour-bam"), type="character", dest='tbam',
              help="input tumour BAM", metavar="bam file name"),
  make_option(c("--normal-bam"), type="character", dest='nbam',
              help="input normal BAM", metavar="bam file name"),
  make_option(c("--potential-variants"), type="character", dest='potential_variants',
              help="variants to investigate", metavar="variants of interest TSV file name"),
  make_option(c("--out"), type="character",
              help="output VCF file name", metavar="VCF file")
)
 
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

#args <- c('processed/keogh2018/SRP159015/SRR7762430/SRR7762430_bqsr.bam',
#'processed/keogh2018/SRP159015/SRR7762431/SRR7762431_bqsr.bam',
#'processed/keogh2018/SRP159015/SRR7762430/SRR7762430_SRR7762431_somatic_varscan_JESROI.tsv')

# the corresponsing .bai files must exist too
tumour_bam = opts$tbam
normal_bam = opts$nbam
potential_variants_file = opts$potential_variants
out_file_name = opts$out
is_one_by_one = opts$one_by_one
if (is.null(tumour_bam)) {
    stop('must specify tumour bam')
}
if (is.null(normal_bam)) {
    stop('must specify normal bam')
}
if (is.null(potential_variants_file)) {
    stop('must specify potential variants')
}
if (is.null(out_file_name)) {
    stop('must specify outfile')
}

regions_of_interest <- read.table(
                        file=potential_variants_file,
                        header=FALSE,
                        col.names =c('CHROM', 'POS', 'ID', 'REF', 'ALT'),
                        colClasses=c('character', 'integer', 'character', 'character', 'character'))
write(sprintf('Potential variants file %s contains %d entries',
                potential_variants_file,
                nrow(regions_of_interest)),
      stderr())

OUT_VCF_COLS = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
# What's a good cutoff?
P_VALUE_CUTOFF = 0.001

# The start and end are inclusive,
# so identical start and stop params gives us single nucleotide regions of interest
regions <- data.frame(chr   = regions_of_interest$CHR,
                      start = regions_of_interest$POS,
                      stop  = regions_of_interest$POS)

suppressPackageStartupMessages(library(deepSNV))
# run the actual analysis
calldeepSNV <- function(regions) {
    # loading the library is slow so only do it just before it's needed,
    # so argument error checks etc can fail early
    write(sprintf("Calling on %d regions of interest", nrow(regions)), stderr())
    # deepSNV fails with only 0 or 1 ROIs, so work around
    if (nrow(regions) == 0) {
        empty_df = data.frame(matrix(nrow=0, ncol=length(OUT_VCF_COLS)))
        colnames(empty_df) = OUT_VCF_COLS
        return(empty_df)
    }
    single_roi_hack = (nrow(regions) == 1)
    if (single_roi_hack) {
        write("Enabling single-ROI deepSNV workaround",  stderr())
        regions = rbind(regions, regions)
    }
    
    ds <- deepSNV(
        test    = tumour_bam,
        control = normal_bam,
        regions = regions,
        q       = 10)
    ds_pvals <- as.data.frame(ds@p.val)
    if (single_roi_hack) {
        ds_pvals <- ds_pvals[1,]
    }

    # The returned dataframe is assumed to refer 1:1 to the input data in the same order.
    # There is no joining on values going on, so it had better be right!
    # So if they're different lengths something very bad has happened.
    if (nrow(regions_of_interest) != nrow(ds_pvals)) {
        stop(sprintf('table of results (%d) has different number of rows to input (%d)',
                        nrow(ds_pvals), nrow(regions_of_interest)))
    }

    combined_df <- cbind(regions_of_interest, ds_pvals)
    combined_df$QUAL <- apply(combined_df, 1, function(row) {'.'})
    return(combined_df)
}

combined_df = calldeepSNV(regions)

for (i in seq_len(nrow(combined_df))) {
    potential_alt <- combined_df[i, 'ALT']
    p_value <- combined_df[i, potential_alt]
    write(sprintf("i = %d", i), stderr())
    write(sprintf("potential_alt = %s", potential_alt), stderr())
    write(sprintf("p_value = %g", p_value), stderr())
    if (nchar(potential_alt) != 1) {
        # should perhaps have a different error for this case
        write("ALT not a single nucleotide, so fail this filter", stderr())
        p_value = NA
    }
    if (is.na(p_value)) {
        # Will get treated as a fail
        combined_df[i, 'INFO'] = '.'
    } else {
        combined_df[i, 'INFO'] = sprintf('deepSNVpValue=%g', p_value)
    }
    if (! is.na(p_value) && p_value <= P_VALUE_CUTOFF) {
        # should you explicitly say PASS, or should it be NULL so it
        # can be merged with pre-existing failures (if present)?
        combined_df[i, 'FILTER'] = 'PASS'
    } else {
        # made-up string for now, is there an expected convention here?
        # Treating an NA pvalue as a fail
        combined_df[i, 'FILTER'] = 'deepSNVFail'
    }
    # write(sprintf("Full row: %s", combined_df[i,]))
}
# Could experiment with summary(deepSNV, value="VCF"), but am doing a DIY VCF for now
# We want no quotes and to start with a #
out_file = file(out_file_name, open = "wt")
cat(sprintf('##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=deepSNVFail,Description="DeepSNV, failed to meet p-value cutoff of %g">
##INFO=<ID=deepSNVpValue,Number=1,Type=Float,Description="The actual p-value issued by deepSNV">
#', P_VALUE_CUTOFF), file=out_file)
cat(OUT_VCF_COLS, sep='\t', file = out_file)
cat('\n', file = out_file)
write.table(
    combined_df[,OUT_VCF_COLS],
    sep='\t',
    file=out_file,
    quote=FALSE,
    col.names=FALSE,
    row.names=FALSE)

close(out_file)