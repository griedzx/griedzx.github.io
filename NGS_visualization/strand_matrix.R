suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(pheatmap))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(forcats))
suppressMessages(library(EnrichedHeatmap))

parser <- ArgumentParser(prog = 'strand_matrix.R',
             description = "helping combind profile data from forward and reverse strand, and visualize it in R.",
             epilog = '')

parser$add_argument('--version', '-v', action = 'version', version = '%(prog)s 1.0.0')
parser$add_argument('--input', '-i', nargs = '+', help = 'the complexHeatmap output file, multiple files should be separated by spaced.', required = TRUE)
parser$add_argument('--output', '-o', help = 'the output file name', default = 'plot.pdf')
# parser$add_argument('--averageType', '-t', help = 'the type of stastics should be used for the profile, "mean" by default.', default = 'mean', choices = c('mean', 'max', 'min', 'median', 'sum'))
# parser$add_argument('--plotType', help = 'the plot type for profile, "line" by default.', default = 'line', choices = c('line', 'heatmap', 'both'))
parser$add_argument('--colors', nargs = '+', help = 'the colors used for plot lines, multiple colors should be separated by spaced and should be equal with group information size, "None" by default.', default = NULL, required = FALSE)

args <- parser$parse_args()
files <- args$input
output <- args$output

data <- lapply(files, FUN = function(file){
           tmp <- read.table(file = file, header = FALSE, sep = "\t", skip = 1)
           #if file_name contains "forward" or "positive", then it is forward strand
              if (grepl("forward|positive", file)){
                tmp$strand <- "forward"
                } else {
                tmp$strand <- "reverse"
                }
})
data <- purrr::reduce(data, rbind)

#line plot
y_label <- "signal(RPKM)"
y_limit <- c(round(min(data$V4)), round(max(data$V4))+1)
y_axis <- seq(y_limit[1], y_limit[2], 1)

x_label <- "Distance to TSS/TTS(bp)"
legend_label <- ""
fill_value <- c("#ea4790")
x_limit <- c(min(data$V1), max(data$V1))
x_axis <- seq(x_limit[1], x_limit[2], 1000)

p <- ggplot(data = data, aes(x = V1, y = V4, color = strand)) +
    geom_line(size = 1) +
    scale_y_continuous(breaks = y_axis, limits = y_limit, expand = c(0, 0)) +
    scale_x_continuous(breaks = x_axis, limits = x_limit, expand = c(0, 0)) +
    labs(x = x_label, y = y_label) +
    theme_classic(base_family = "Arial",
                    base_line_size = 0.3) +
    theme(legend.position = "none")

EnrichedHeatmap(data, )