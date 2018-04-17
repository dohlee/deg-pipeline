# Install required packages if they are not installed.
source("https://bioconductor.org/biocLite.R")
require('EBSeq', quietly=TRUE) || biocLite('EBseq')
require('optparse', quietly=TRUE) || install.packages('optparse', repos='http://cran.us.r-project.org')

# Parse options.
option_list = list(
  make_option(c('-i', '--input'), type='character', default=NULL, help='Input file name.'),
  make_option(c('-c', '--condition'), type='character', default=NULL, help='Input condition file name.'),
  make_option(c('-o', '--outdir'), type='character', default='result/ebseq', help='Output directory for EBSeq results. [default=%default]')
)
parser = OptionParser(option_list=option_list)
args = parse_args(parser)

# Error handling for options.
if (is.null(args$input)) {
  print_help(parser)
  stop("Input file name[-i, --input] should be specified.", call.=FALSE)
}

if (is.null(args$condition)) {
  print_help(parser)
  stop("Condition file name[-c, --condition] should be specified.", call.=FALSE)
}
