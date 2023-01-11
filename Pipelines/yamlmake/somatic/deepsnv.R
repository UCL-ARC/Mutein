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
  make_option(c("--exclude-bad"), type="logical", action="store_true", default=FALSE,  dest='exclude_bad',
              help="Exclude regions known to cause a crash"),
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
# fix the IDs in a column (row names get rewritten by semi_join)
regions_of_interest$roi_id = rownames(regions_of_interest)
write(sprintf('Potential variants file %s contains %d entries',
                potential_variants_file,
                nrow(regions_of_interest)),
      stderr())

OUT_VCF_COLS = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
P_VALUE_CUTOFF = 0.001

# This list is dataset specific, and determined experimentally using the --one-by-one option to this script.
# The crash logs can then be grepped as follows to identify the position(s) causing a crash.
# Only the first crashing position per sample will be detected, therefore after adding rows
# to bad_data, you will need to re-run everything, and you will get a (hopefully) smaller set of crashing rows to
# add, until you've caught every one.
# grep -h -B 2 'caught segfault'  yamlmake_logs/mut.20221221.172252.736840.deepsnv.*.err | grep Calling | sort | uniq

# Ideally I would find out the root cause of the crashing (maybe some property of
# the position is evident in the BAM file so bad locations could be automatically identified?)
bad_data = rbind(
    data.frame(CHR='chr1', POS=20633892),
    data.frame(CHR='chr1', POS=20645555),
    data.frame(CHR='chr1', POS=20645618),
    data.frame(CHR='chr1', POS=226883824),
    data.frame(CHR='chr2', POS=29320797),
    data.frame(CHR='chr2', POS=74370400),
    data.frame(CHR='chr2', POS=201711098),
    data.frame(CHR='chr4', POS=95206718),
    data.frame(CHR='chr5', POS=70067122),
    data.frame(CHR='chr5', POS=150120109),
    data.frame(CHR='chr6', POS=117301070),
    data.frame(CHR='chr6', POS=117383318),
    data.frame(CHR='chr6', POS=151879295),
    data.frame(CHR='chr6', POS=152061247),
    data.frame(CHR='chr7', POS=154052934),
    data.frame(CHR='chr7', POS=154880936),
    data.frame(CHR='chr7', POS=154889340),
    data.frame(CHR='chr9', POS=136513000),
    data.frame(CHR='chr10', POS=13125860),
    data.frame(CHR='chr10', POS=43120057),
    data.frame(CHR='chr11', POS=32435148),
    data.frame(CHR='chr11', POS=121605089),
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
    data.frame(CHR='chr21', POS=38403636),
    data.frame(CHR='chr21', POS=41473456),
    data.frame(CHR='chr21', POS=41507982),
    data.frame(CHR='chr22', POS=23766225),
    data.frame(CHR='chr22', POS=23767587),
    data.frame(CHR='chr22', POS=32479203)

# All files ever!
# $ grep -hr -B 2 'caught segfault' yamlmake_logs/mut.*.*.err  | grep 'region:' |sort|uniq -c
#       2 Calling on region: chr10 13125860
#       6 Calling on region: chr10 43120057
#      40 Calling on region: chr11 32435148
#       8 Calling on region: chr11 64805130
#       4 Calling on region: chr1 20633892
#       6 Calling on region: chr1 20645555
#       8 Calling on region: chr1 20645618
#       4 Calling on region: chr12 40309044
#      12 Calling on region: chr12 40320071
#       6 Calling on region: chr12 40364850
#       6 Calling on region: chr14 36520594
#       8 Calling on region: chr15 34348017
#       6 Calling on region: chr15 34348177
#      16 Calling on region: chr16 31182621
#      14 Calling on region: chr16 46662372
#       2 Calling on region: chr16 68737447
#      12 Calling on region: chr17 39723335
#       8 Calling on region: chr17 45993928
#       4 Calling on region: chr19 15160960
#      12 Calling on region: chr19 15174241
#       4 Calling on region: chr19 15184323
#      12 Calling on region: chr21 34887027
#       6 Calling on region: chr21 38403636
#       4 Calling on region: chr21 41473456
#       8 Calling on region: chr21 41507982
#       4 Calling on region: chr2 201711098
#       4 Calling on region: chr22 23766225
#       6 Calling on region: chr22 23767587
#      24 Calling on region: chr22 32479203
#       8 Calling on region: chr2 29320797
#       8 Calling on region: chr2 74370400
#      20 Calling on region: chr4 95206718
#       4 Calling on region: chr5 150120109
#       4 Calling on region: chr5 70067122
#       2 Calling on region: chr6 117301070
#       4 Calling on region: chr6 152061247
#      12 Calling on region: chr7 154889340

)

# use the roi_id because semi_join generates rows with new row names
bad_rownums = semi_join(regions_of_interest, bad_data, 
                                 by = c('CHROM' = 'CHR', 'POS' = 'POS'))$roi_id

if (length(bad_rownums) > 0) {
    if (opts$exclude_bad) {
        warning(sprintf('Omitting %d rows of bad data known to cause segfault in R:\n', length(bad_rownums)),
            paste0(capture.output(regions_of_interest[bad_rownums,]), collapse='\n'))
        regions_of_interest = regions_of_interest[!rownames(regions_of_interest) %in% bad_rownums,]
    } else {
        warning(sprintf('Identified %d rows of bad data known to cause segfault in R, ', length(bad_rownums)),
                'but not omitting them due to command line options. See help for more.')
    }
}

# finished with this
regions_of_interest$roi_id <- NULL

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
        write(sprintf("Calling on %d regions of interest", nrow(regions)), stderr())
        # deepSNV fails with only a single ROI, so work around
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
    }
    # What's a good cutoff?

    # what format does bcftools annotate need?
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
        write("ALT not a single nucleotide, will fail filter", stderr())
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
# Including VCF output through summary(deepSNV, value=”VCF”)


# combined_df$FILTER <- apply(combined_df, 1,
#     function(row) {
#         p_value_cutoff = P_VALUE_CUTOFF
#         potential_alt <- row['ALT']
#         if (! is.na(row[potential_alt]) && row[potential_alt] <= p_value_cutoff) {
#             # should you explicitly say PASS, or should it be NULL so it
#             # can be merged with pre-existing failures (if present)?
#             'PASS'
#         } else {
#             # made-up string for now, is there an expected convention here?
#             'deepSNVFail'
#         }
#     })
# combined_df$INFO <- apply(combined_df, 1, function(row) {'.'})

# diy vcf
# we want no quotes and to start with a #
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