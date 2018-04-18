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
EB_out = EBTest(Data=data_mat, Conditions=as.factor(conditions), sizeFactors=size_factors, maxround=5)
flog.info('Running GetDEResults...')
EB_DE_result = GetDEResults(EB_out, FDR=0.05)

flog.info('Writing to output file...')
dir.create(args$outdir, recursive=TRUE, showWarnings=FALSE)
write.table(EB_DE_result$DEfound, file.path(args$outdir, 'ebseq_degs.txt'), quote=F, row.names=F, col.names=F)
