# Install futile.logger package if it is not installed.
require('futile.logger', quietly=TRUE) || install.packages('futile.logger', repos='http://cran.us.r-project.org')

set_logger = function(script_name) {
    # Layout example '[EBSeq|2018-04-18 11:04:44] Running EBSeq...'
    layout = layout.format(paste('[', script_name, '|~t] ~m'))
    flog.layout(layout)
}
