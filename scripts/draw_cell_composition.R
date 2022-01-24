# Draw the deconvoluted cell composition of selected patients

library(tidyr)

load("../data/tcga_time_deconv.RData")

res.deconv.epic$cell_type <- rownames(res.deconv.epic)

p <- res.deconv.epic[, 401:406] %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_brewer(palette="Paired") +
  scale_x_discrete(limits = rev(levels(res.deconv.epic)))

pdf("../output/example_cell_composition.pdf", width = 8, height = 4)
print(p)
dev.off()
