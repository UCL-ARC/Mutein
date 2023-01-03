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
  make_option(c("--one-by-one"), type="logical", action="store_true", default=FALSE,  dest='one_by_one',
              help="For debugging only, run areas of interest one by one to identify regions that cause a crash"),
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
                        col.names=c('CHROM', 'POS', 'ID', 'REF', 'ALT'))

# grep -h -B 2 'caught segfault'  yamlmake_logs/mut.20221221.172252.736840.deepsnv.*.err | grep Calling | sort | uniq

# grep -B4 'caught segfault' yamlmake_logs/mut.20221221.154521.512993.deepsnv.*.err | view -
bad_data = rbind(
    data.frame(CHR='chr1',  POS=20633892),
    data.frame(CHR='chr1',  POS=20645555),
    data.frame(CHR='chr1',  POS=20645618),
    data.frame(CHR='chr1',  POS=226883824),
    data.frame(CHR='chr2',  POS=201711098),
    data.frame(CHR='chr2',  POS=29320797),
    data.frame(CHR='chr2',  POS=74370400),
    data.frame(CHR='chr4',  POS=95206718),
    data.frame(CHR='chr5',  POS=70067122),
    data.frame(CHR='chr6',  POS=117301070),
    data.frame(CHR='chr6',  POS=117383318),
    data.frame(CHR='chr6',  POS=151879295),
    data.frame(CHR='chr6',  POS=152061247),
    data.frame(CHR='chr7',  POS=154052934),
    data.frame(CHR='chr7',  POS=154880936),
    data.frame(CHR='chr7',  POS=154889340),
    data.frame(CHR='chr9',  POS=136513000),
    data.frame(CHR='chr10', POS=13125860),
    data.frame(CHR='chr11', POS=121605089),
    data.frame(CHR='chr11', POS=32435148),
    data.frame(CHR='chr11', POS=64805130),
    data.frame(CHR='chr12', POS=40309044),
    data.frame(CHR='chr12', POS=40320071),
    data.frame(CHR='chr12', POS=40364850),
    data.frame(CHR='chr14', POS=36520594),
    data.frame(CHR='chr15', POS=34348017),
    data.frame(CHR='chr15', POS=34348177),
    data.frame(CHR='chr16', POS=31182621),
    data.frame(CHR='chr16', POS=46662372),
    data.frame(CHR='chr16', POS=68737447),
    data.frame(CHR='chr17', POS=39723335),
    data.frame(CHR='chr17', POS=45993928),
    data.frame(CHR='chr19', POS=15160960),
    data.frame(CHR='chr19', POS=15174241),
    data.frame(CHR='chr19', POS=15184323),
    data.frame(CHR='chr21', POS=34887027),
    data.frame(CHR='chr21', POS=41473456),
    data.frame(CHR='chr21', POS=41507982),
    data.frame(CHR='chr22', POS=32479203)

# Next go:
#  chr10 13125860
#  chr11 32435148
#  chr11 64805130
#  chr1 20633892
#  chr1 20645555
#  chr1 20645618
#  chr12 40309044
#  chr12 40320071
#  chr12 40364850
#  chr14 36520594
#  chr15 34348017
#  chr15 34348177
#  chr16 31182621
#  chr16 46662372
#  chr16 68737447
#  chr17 39723335
#  chr17 45993928
#  chr19 15160960
#  chr19 15174241
#  chr19 15184323
#  chr21 34887027
#  chr21 38403636
#  chr21 41473456
#  chr21 41507982
#  chr22 23767587
#  chr22 32479203
#  chr2 29320797
#  chr2 74370400
#  chr4 95206718
#  chr5 70067122



)
bad_rownums = which(apply(regions_of_interest, 1,
function(roi_row) {
    return(any(bad_data$CHR == roi_row['CHROM'] & bad_data$POS == roi_row['POS']))
}))

# bad_rownums = which(regions_of_interest$CHR == 'chr10' & regions_of_interest$POS == '13125860')
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

# loading the library is slow so only do it just before it's needed,
# so argument error checks etc can fail early
suppressPackageStartupMessages(library(deepSNV))

# run the actual analysis
if (is_one_by_one) {
    ds_pvals <- NULL
    for (rn in rownames(regions)) {
        single_region = regions[rn,]
        roi_str = sprintf("region: %s %s", single_region['chr'], single_region['start'])
        write(sprintf("Calling on %s", roi_str), stderr())
#         Error in (function (classes, fdef, mtable)  :
#   unable to find an inherited method for function ‘consensusSequence’ for signature ‘"numeric"’
# Calls: deepSNV ... .local -> .deepSNV -> consensusSequence -> <Anonymous>
# Execution halted

        dsrow <- deepSNV(
            test    = tumour_bam,
            control = normal_bam,
            # deepSNV fails mysteriously if a single row is given here,
            # so double up the row we care about, then remove the double results later
            regions = rbind(single_region, single_region),
            q       = 10)
        write(sprintf("Done calling on %s", roi_str), stderr())
        # init or append
        if (is.null(ds_pvals)) {
            ds_pvals <- dsrow@p.val[1,]
        } else {
            ds_pvals <- rbind(ds_pvals, dsrow@p.val[1,])
        }
    }
} else {
    ds <- deepSNV(
        test    = tumour_bam,
        control = normal_bam,
        regions = regions,
        q       = 10)
    ds_pvals <- ds@p.val
}

# What's a good cutoff?

P_VALUE_CUTOFF = 0.001
# what format does bcftools annotate need?
combined_df <- cbind(regions_of_interest, ds_pvals)
combined_df$QUAL <- apply(combined_df, 1, function(row) {'.'})
combined_df$FILTER <- apply(combined_df, 1,
    function(row) {
        p_value_cutoff = P_VALUE_CUTOFF
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
# quick and dirty vcf
# we want no quotes and to start with a #
out_file = file(out_file_name, open = "wt")
cat(sprintf('##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=deepSNVFail,Description="DeepSNV, failed to meet p-value cutoff of %g">
#', P_VALUE_CUTOFF), file=out_file)
cat(vcf_cols, sep='\t', file = out_file)
cat('\n', file = out_file)
write.table(
    combined_df[,vcf_cols],
    sep='\t',
    file=out_file,
    quote=FALSE,
    col.names=FALSE,
    row.names=FALSE)

close(out_file)