# Install optparse package if it is not installed.
require('optparse', quietly=TRUE) || install.packages('optparse', repos='http://cran.us.r-project.org')

# Import logger and set the name of the logger.
source('logging.R')
set_logger(name='EBSeq')

# Parse options.
option_list = list(
  make_option(c('-i', '--input'), type='character', default=NULL, help='Input file name.'),
  make_option(c('-c', '--condition'), type='character', default=NULL, help='Input condition file name.'),
  make_option(c('-o', '--outdir'), type='character', default='result/ebseq', help='Output directory for EBSeq results. [default=%default]')
)
parser = OptionParser(option_list=option_list)
args = parse_args(parser)

# Install required packages if they are not installed.
source("https://bioconductor.org/biocLite.R")
require('EBSeq', quietly=TRUE) || biocLite('EBseq')

# Error handling for options.
if (is.null(args$input)) {
  print_help(parser)
  stop("Input file name[-i, --input] should be specified.", call.=FALSE)
}

if (is.null(args$condition)) {
  print_help(parser)
  stop("Condition file name[-c, --condition] should be specified.", call.=FALSE)
}

# Get input data and convert it into a matrix.
flog.info('Parsing input data...')
input_data = read.table(args$input, row.names=1, header=TRUE, check.names=FALSE)
data_mat = data.matrix(input_data)

# Parse conditions.
flog.info('Parsing condition table...')
conditions = read.table(args$condition, stringsAsFactors=FALSE, row.names=1, header=FALSE)

# Only take samples in intersection.
flog.info('Taking samples in intersection...')
common_samples = intersect(colnames(data_mat), rownames(conditions))
data_mat = data_mat[, common_samples]
conditions = conditions[common_samples, ]
flog.info('Total %d samples are going to be analyzed.', length(common_samples))

# Compute size factors.
flog.info('Computing size factors...')
size_factors = MedianNorm(data_mat)

# Discover DEGs.
flog.info('Discovering DEGS...')
flog.info('Running EBTest...')

iteration = 0
repeat {
  # Increase iteration numbers if the conditions are not met.
  # Hopefully most of the time, 10 iterations will be enough for convergence.
  iteration = iteration + 10
  EB_out = EBTest(Data=data_mat, Conditions=as.factor(conditions), sizeFactors=size_factors, maxround=5)
  flog.info('Running GetDEResults...')
  EB_DE_result = GetDEResults(EB_out, FDR=0.05)

  # Check convergences.
  # Each parameter should change less than 1e-3 between the last two iterations.
  l = length(EB_out$Alpha)
  if (!abs(EB_out$Alpha[l] - EB_out$Alpha[l-1]) < 1e-3) continue
  if (!abs(EB_out$Beta[l] - EB_out$Beta[l-1]) < 1e-3) continue
  if (!abs(EB_out$P[l] - EB_out$P[l-1]) < 1e-3) continue
  # If all conditions are met, break the loop.
  break
}

# Write DEGs to output file.
flog.info('Writing a DEG list to %s', file.path(args$outdir, 'ebseq_degs.txt'))
dir.create(args$outdir, recursive=TRUE, showWarnings=FALSE)
write.table(EB_DE_result$DEfound, file.path(args$outdir, 'ebseq_degs.txt'), quote=F, row.names=F, col.names=F)

# Write fold-change result to output file.
flog.info('Writing a fold-change list to %s', file.path(args$outdir, 'ebseq_fcs.txt'))
fold_changes = PostFC(EB_out)
write.table(fold_changes, file.path(args$outdir, 'ebseq_fcs.txt'), quote=F)

# Save RealFC vs PosteriorFC plot.
flog.info('Saving a RealFC vs PosteriorFC plot to %s', file.path(args$outdir, 'ebseq_realfc_vs_posteriorfc.png'))
png(file.path(args$outdir, 'ebseq_realfc_vs_posteriorfc.png'))
PlotPostVsRawFC(EB_out, fold_changes)
dev.off()

# Save Q-Q plot to check the model fit.
flog.info('Saving a Q-Q plot to %s', file.path(args$outdir, 'ebseq_qqplot.png'))
png(file.path(args$outdir, 'ebseq_qqplot.png'), width=960, height=480)
par(mfrow=c(1, 2))
QQP(EB_out)
dev.off()

# Save density fit plot to check the model fit.
flog.info('Saving a density fit plot to %s', file.path(args$outdir, 'ebseq_density_fit.png'))
png(file.path(args$outdir, 'ebseq_density_fit.png'), width=960, height=480)
par(mfrow=c(1, 2))
DenNHist(EB_out)
dev.off()
