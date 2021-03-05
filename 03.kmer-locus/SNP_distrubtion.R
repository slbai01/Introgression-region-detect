library(ggplot2)
library(dplyr)
library(stringr)
library('ggsci')

args=commandArgs(T)

input_file <- args[1]
out_file <- args[2]
y_lim <- as.numeric(args[3])

dat <- read.table(input_file, sep = '\t',
                  header = F, stringsAsFactors = F)

colnames(dat) = c("CHROM", "start", "end", "depth")

dat1 = dat %>% filter(str_detect(CHROM,pattern = "^[Cc]hr"))
dat1$CHROM = factor(dat1$CHROM)

png(out_file)
ggplot(dat1, aes(x = start, y = depth, colour = CHROM)) + geom_line() +
  scale_x_continuous(breaks = seq(0,650000000,100000000)) +
  scale_y_continuous(limits = c(0, y_lim)) +
  xlab("") +
  scale_color_npg() +
  facet_wrap(~ CHROM, ncol = 1) + theme_bw()
dev.off()


