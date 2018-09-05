library(dplyr)
library(ggplot2)
library(cowplot)

my_theme <- function() {
  theme_minimal(base_family = "Helvetica") %+replace% 
    theme(title = element_text(size = 7),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.margin = unit(c(0.21, 0.21, 0.21, 0.21), "lines"),
          axis.line = element_line(colour = "grey20", size = 0.5), 
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 0),
          legend.key = element_blank(), 
          legend.key.size = unit(0.25, "cm"),
          legend.margin = margin(0.5, 0.5, 0.5, 0.5),
          legend.text = element_text(size = 6),
          strip.text = element_text(size = 6),
          # strip.background = element_rect(fill = "white", colour = "black", size = 1),
          strip.background = element_blank(), 
          complete = TRUE)
} 


# a <- read.table("../matome/omega.alt_junc.txt", sep = "\t", header = TRUE, quote = "")
a <- read.table("../output/intron_mut_annot/savnet.creatin_mutation.annot.txt", sep = "\t", header = FALSE, quote = "") %>%
  select(Is_Donor = V9, Dist = V11)

a$Is_Donor <- factor(a$Is_Donor, levels = c("donor", "acceptor"), label = c("Donor", "Acceptor"))

a$Log_Dist <- unlist(lapply(a$Dist, function(x) {ifelse(x, log10(x), 0)}))

p_d <- ggplot(a %>% filter(Is_Donor == "Donor"), aes(x = Log_Dist)) + 
  geom_histogram(binwidth = 0.2, fill = "#fdb462", colour = "grey30", size = 0.3) +
  ggtitle("Donor") +
  my_theme() +
  scale_y_continuous(limits = c(0, 135), expand = c(0, 0)) +
  labs(x = "", y = "") # +
  # labs(x = "Log10(distance from nearest exon)", y = "SAV Frequency")


p_a <- ggplot(a %>% filter(Is_Donor == "Acceptor"), aes(x = Log_Dist)) + 
  geom_histogram(binwidth = 0.2, fill = "#fdb462", colour = "grey30", size = 0.3) +
  ggtitle("Acceptor") +
  my_theme() +
  scale_y_continuous(limits = c(0, 135), expand = c(0, 0)) +
  labs(x = "", y = "") # +
  # labs(x = "Log10(distance from nearest exon)", y = "SAV Frequency")




p_d_a <- plot_grid(p_d, p_a, ncol = 2, align = "h")

xlabel <- ggdraw() + draw_label("Log10(distances from the nearest exons)", size = 7)
ylabel <- ggdraw() + draw_label("SAV Frequency", angle = 90, size = 7)

p_d_a_xl <- plot_grid(p_d_a, xlabel, ncol = 1, align = "v", rel_heights = c(1, 0.1))


plot_grid(ylabel, p_d_a_xl, ncol = 2, align = "h", rel_widths = c(0.05, 1))


# ggsave("../matome/alt_junc_pos.pdf", width = 8, height = 6)
ggsave("../output/intron_mut_annot/alt_junc_pos_log.pdf", width = 8, height = 4.5, units = "cm")



