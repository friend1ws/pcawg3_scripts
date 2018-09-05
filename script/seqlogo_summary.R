library(dplyr)
# library(seqLogo)
library(ggplot2)
library(ggseqlogo)
library(cowplot)
library(Cairo)

Cairo()

source("plot_config.R")

donor_pos_label <- c(add_emdash("3"), add_emdash("2"), add_emdash("1"), "+1", "+2", "+3", "+4", "+5", "+6")
acceptor_pos_label <- c("+6", "+5", "+4", "+3", "+2", "+1", add_emdash("1"))


splicing_mutation <- read.table("../output/motif/pcawg.savnet.motif_summary.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE) %>%
  filter(Mutation_Type %in% c("splicing donor creation", "splicing acceptor creation"))


theme_nonbottom <- function() {
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "lines"),
          axis.text = element_text(size = 7),
          # axis.title = element_text(size = 7),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title = element_blank())
}


theme_bottom <- function() {
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "lines"),
          panel.border = element_blank(),
          axis.text = element_text(size = 7),
          # axis.title = element_text(size = 7),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.line.x = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4)) # +
     # scale_x_continuous(breaks = 1:9,
     #                    labels = c("-3", "-2", "-1", "1", "2", "3", "4", "5", "6"))

}

get_base_count_mat <- function(motif_seqs, nbase) {

  base_count_mat <- matrix(0, nbase, 4)
  for(n in 1:nbase) {
    base_count <- table(substring(motif_seqs, n + 0, n + 0))
    base_count_mat[n, 1] <- base_count["A"]
    base_count_mat[n, 2] <- base_count["C"]
    base_count_mat[n, 3] <- base_count["G"]
    base_count_mat[n, 4] <- base_count["T"]
  }
  base_count_mat[is.na(base_count_mat)] <- 0  
  return(base_count_mat)
}

theme_bottom <- function() {
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "lines"),
          panel.border = element_blank(),
          axis.text = element_text(size = 7),
          # axis.title = element_text(size = 7),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.line.x = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4)) # +
}


motif_sub <- splicing_mutation %>% filter(Mutation_Type == "splicing donor creation")

motif_wt_d <- motif_sub$Motif_Seq

logo_wt_d <- ggseqlogo(t(get_base_count_mat(motif_wt_d, 9)), TRUE) +
    labs(x = "", y = "Bits") +
    scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0)) +
    theme_bottom() +
    scale_x_continuous(breaks = 1:9, labels = donor_pos_label)


motif_mt_d <- unlist(lapply(1:length(motif_wt_d),
                              function(x) {
                                paste(substring(motif_wt_d[x], 1, motif_sub$Rel_Pos[x] - 1),
                                      motif_sub$Alt_Base[x],
                                      substring(motif_wt_d[x], motif_sub$Rel_Pos[x] + 1, 9),
                                      sep ="")
                              }))


logo_mt_d <- ggseqlogo(t(get_base_count_mat(motif_mt_d, 9)), TRUE) +
    labs(x = "", y = "Bits") +
    scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0)) +
    theme_bottom() +
    scale_x_continuous(breaks = 1:9, labels = donor_pos_label)



motif_sub <- splicing_mutation %>% filter(Mutation_Type == "splicing acceptor creation")

motif_wt_a <- motif_sub$Motif_Seq

logo_wt_a <- ggseqlogo(t(get_base_count_mat(motif_wt_a, 7)), TRUE) +
    labs(x = "", y = "Bits") +
    scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0)) +
    theme_bottom() +
    scale_x_continuous(breaks = 1:7, labels = acceptor_pos_label)


motif_mt_a <- unlist(lapply(1:length(motif_wt_a),
                              function(x) {
                                paste(substring(motif_wt_a[x], 1, motif_sub$Rel_Pos[x] - 1),
                                      motif_sub$Alt_Base[x],
                                      substring(motif_wt_a[x], motif_sub$Rel_Pos[x] + 1, 9),
                                      sep ="")
                              }))


logo_mt_a <- ggseqlogo(t(get_base_count_mat(motif_mt_a, 7)), TRUE) +
    labs(x = "", y = "Bits") +
    scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0)) +
    theme_bottom() +
    scale_x_continuous(breaks = 1:7, labels = acceptor_pos_label)



plot_grid(ggdraw() + draw_label(""),
          ggdraw() + draw_label("Donor creation", size = 7),
          ggdraw() + draw_label("Acceptor creation", size = 7),
          ggdraw() + draw_label("WT", size = 7),
          logo_wt_d, logo_wt_a, 
          ggdraw() + draw_label("MT", size = 7),
          logo_mt_d, logo_mt_a, ncol = 3, rel_heights = c(0.2, 1, 1), rel_widths = c(0.2, 1, 0.9))

ggsave("../output/motif/motif_comp.pdf", width = 8, height = 5, units = "cm", dpi = 600) 


