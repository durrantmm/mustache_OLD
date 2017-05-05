library(readr); packageVersion("readr")
library(ggplot2); packageVersion("ggplot2")
library(tools)

args = commandArgs(trailingOnly=TRUE)

peak_depths_folder = args[1]
output_folder = args[2]

dir.create(output_folder)

peak_depths <- list.files(peak_depths_folder, full.names="TRUE")

for (filein in peak_depths){
  depths <- read_tsv(filein)
  gg <- qplot(data=depths, x=POS, y=VALID_DEPTH, geom='line', color='Total Read Depth') +
    geom_line(aes(x=POS, y=REVERSE_DEPTH, color='Reverse Read Depth'), linetype='dashed') +
    geom_line(aes(x=POS, y=FORWARD_DEPTH, color='Forward Read Depth'), linetype='dotted') +
    scale_color_manual(values=c("red", "blue", "black")) +
    xlab("Genomic Position") +
    ylab("Read Depth") +
    theme_bw()
  outfile <- file.path(output_folder, basename(filein))
  outpdf <- paste(file_path_sans_ext(outfile), '.pdf', sep='')
  ggsave(outpdf, width=5, height=3)
}

