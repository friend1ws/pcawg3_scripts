library(dplyr)
library(ggplot2)

# source("../../../conf/plot_config.R")

my_theme <- function() {
   theme_bw(base_family = "Helvetica") %+replace% 
   theme(title = element_text(size = 7),
         panel.border = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "grey20", size = 0.5), 
         axis.text = element_text(size = 7),
         axis.title = element_text(size = 7),
         legend.key = element_blank(), 
         legend.key.size = unit(0.25, "cm"),
         legend.margin = margin(0.5, 0.5, 0.5, 0.5),
         legend.text = element_text(size = 6),
         strip.text = element_text(size = 6),
         # strip.background = element_rect(fill = "white", colour = "black", size = 1),
         strip.background = element_blank(), 
         complete = TRUE)
} 

  
cacner_type_fdr <- read.table("../output/fdr/cancer_fdr.txt", sep = "\t", header = TRUE)

ggplot(cacner_type_fdr, aes(x = Cancer_Type, y = FDR)) +
  geom_bar(stat = "identity", fill = "#4daf4a") +
  my_theme() +
  labs(x = "Cancer type") +
  geom_abline(intercept = 0.10, slope = 0, colour = "#d73027", alpha = 0.8, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0))


ggsave("../output/fdr/cancer_type_fdr.pdf", width = 10, height = 5, dpi = 600, units = "cm")

