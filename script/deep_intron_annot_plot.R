library(dplyr)
library(ggplot2)

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



D <- read.table("../output/intron_mut_annot/deep_mut_annot_summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

D2 <- D %>% group_by(Cancer_Type) %>% 
  summarize(Deep_Mut_Count = sum(Deep_Mut_Count), 
            CG1_Count = sum(CG1_Count),
            CG2_Count = sum(CG2_Count),
            CG3_Count = sum(CG3_Count),
            SINE_Count = sum(SINE_Count),
            SINE_Sense_Count = sum(SINE_Sense_Count),
            SINE_Antisense_Count = sum(SINE_Antisense_Count)) %>%
  mutate(CG1_Ratio = CG1_Count / Deep_Mut_Count,
         CG2_Ratio = CG2_Count / Deep_Mut_Count,
         CG3_Ratio = CG3_Count / Deep_Mut_Count,
         SINE_Ratio = SINE_Count / Deep_Mut_Count,
         SINE_Sense_Ratio = SINE_Sense_Count / Deep_Mut_Count,
         SINE_Antisense_Ratio = SINE_Antisense_Count / Deep_Mut_Count)


ggplot(D2, aes(x = Cancer_Type, y = CG1_Ratio)) + 
  geom_bar(stat = "identity") + 
  my_theme() +
  labs(x = "Cancer Type", y = "Lawrence et al. (2014)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../output/intron_mut_annot/LawrenceEtAl2014_ratio_pancan.pdf", width = 18, height = 6, units = "cm")


ggplot(D2, aes(x = Cancer_Type, y = CG2_Ratio)) + 
  geom_bar(stat = "identity") + 
  my_theme() +
  labs(x = "Cancer Type", y = "Cancer Genome Census") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../output/intron_mut_annot/CancerGenomeCensus_ratio_pancan.pdf", width = 18, height = 6, units = "cm")


ggplot(D2, aes(x = Cancer_Type, y = CG3_Ratio)) + 
  geom_bar(stat = "identity") + 
  my_theme() +
  labs(x = "Cancer Type", y = "Ye et al. (2016)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../output/intron_mut_annot/YeEtAl2016_pancan.pdf", width = 18, height = 6, units = "cm")


ggplot(D2, aes(x = Cancer_Type, y = SINE_Ratio)) + 
  geom_bar(stat = "identity") + 
  my_theme() +
  labs(x = "Cancer Type", y = "SINE") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../output/intron_mut_annot/SINE_pancan.pdf", width = 18, height = 6, units = "cm")


ggplot(D2, aes(x = Cancer_Type, y = SINE_Sense_Ratio)) + 
  geom_bar(stat = "identity") + 
  my_theme() +
  labs(x = "Cancer Type", y = "SINE (sense)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../output/intron_mut_annot/SINE_sense_pancan.pdf", width = 18, height = 6, units = "cm")



ggplot(D2, aes(x = Cancer_Type, y = SINE_Antisense_Ratio)) + 
  geom_bar(stat = "identity") + 
  my_theme() +
  labs(x = "Cancer Type", y = "SINE (antisense)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../output/intron_mut_annot/SINE_antisense_pancan.pdf", width = 18, height = 6, units = "cm")


write.table(S2 %>% filter(Annotation %in% c("SINE", "SINE sense", "SINE antisense")),
            "../result/intron_mut_annot/SINE_ratio_diff.txt", quote = FALSE, sep = "\t", row.names = FALSE)


D3 <- D2 %>% group_by() %>% 
  summarize(Deep_Mut_Count = sum(Deep_Mut_Count), 
            CG1_Count = sum(CG1_Count),
            CG2_Count = sum(CG2_Count),
            CG3_Count = sum(CG3_Count),
            SINE_Count = sum(SINE_Count),
            SINE_Sense_Count = sum(SINE_Sense_Count),
            SINE_Antisense_Count = sum(SINE_Antisense_Count)) %>%
  mutate(CG1_Ratio = CG1_Count / Deep_Mut_Count,
         CG2_Ratio = CG2_Count / Deep_Mut_Count,
         CG3_Ratio = CG3_Count / Deep_Mut_Count,
         SINE_Ratio = SINE_Count / Deep_Mut_Count,
         SINE_Sense_Ratio = SINE_Sense_Count / Deep_Mut_Count,
         SINE_Antisense_Ratio = SINE_Antisense_Count / Deep_Mut_Count)





S <- read.table("../output/intron_mut_annot/savnet.creatin_mutation.annot_summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(CG1_Ratio = CG1_Count / Deep_Mut_Count,
         CG2_Ratio = CG2_Count / Deep_Mut_Count,
         CG3_Ratio = CG3_Count / Deep_Mut_Count,
         SINE_Ratio = SINE_Count / Deep_Mut_Count,
         SINE_Sense_Ratio = SINE_Sense_Count / Deep_Mut_Count,
         SINE_Antisense_Ratio = SINE_Antisense_Count / Deep_Mut_Count)


Annotation <- factor(c("Lawrence et al., (2014)", "Lawrence et al., (2014)",
                "Cancer Gene Census", "Cancer Gene Census",
                "Ye et al., (2016)", "Ye et al., (2016)",
                "SINE", "SINE", "SINE sense", "SINE sense",
                "SINE antisense", "SINE antisense"),
                levels = c("Lawrence et al., (2014)", "Cancer Gene Census",
                           "Ye et al., (2016)", "SINE", 
                           "SINE sense", "SINE antisense"),
                labels = c("Lawrence (2014)", "Cancer Gene Census",
                           "Ye (2016)", "SINE", 
                           "SINE sense", "SINE antisense")
)

Ratio <- c(S$CG1_Ratio, D3$CG1_Ratio,
           S$CG2_Ratio, D3$CG2_Ratio,
           S$CG3_Ratio, D3$CG3_Ratio,
           S$SINE_Ratio, D3$SINE_Ratio,
           S$SINE_Sense_Ratio, D3$SINE_Sense_Ratio,
           S$SINE_Antisense_Ratio, D3$SINE_Antisense_Ratio
           )

IsSAVNET <- rep(c("SAV", "Backgroud"), 6)
          
                

S2 <- data.frame(Annotation = Annotation,
                 Ratio = Ratio, IsSAVNET = IsSAVNET)

ggplot(S2 %>% filter(Annotation %in% 
                       c("Lawrence (2014)", "Cancer Gene Census", "Ye (2016)")), 
       aes(x = IsSAVNET, y = Ratio, fill = IsSAVNET)) + 
  geom_bar(stat = "identity") +
  # theme_bw() +
  my_theme() +
  coord_flip() +
  facet_grid(Annotation~.) +
  ylim(c(0, 0.1)) +
  labs(x = "", y = "Cancer gene ratio") +
  guides(fill = FALSE) +
  scale_fill_manual(values = c("#666666", "#7fc97f"))

ggsave("../output/intron_mut_annot/CG_ratio_diff.pdf", width = 6, height = 9, units = "cm")


ggplot(S2 %>% filter(Annotation %in% 
                       c("SINE", "SINE sense", "SINE antisense")), 
       aes(x = IsSAVNET, y = Ratio, fill = IsSAVNET)) + 
  geom_bar(stat = "identity") +
  # theme_bw() +
  my_theme() +
  coord_flip() +
  facet_grid(Annotation~.) +
  ylim(c(0, 0.3)) +
  labs(x = "", y = "Fraction of overlap with Alu sequence") +
  guides(fill = FALSE) +
  scale_fill_manual(values = c("#666666", "#7fc97f"))


ggsave("../output/intron_mut_annot/SINE_ratio_diff.pdf", width = 6, height = 8, units = "cm")



1 - pbinom(S$CG1_Count, S$Deep_Mut_Count, D3$CG1_Ratio)
1 - pbinom(S$CG2_Count, S$Deep_Mut_Count, D3$CG2_Ratio)
1 - pbinom(S$CG3_Count, S$Deep_Mut_Count, D3$CG3_Ratio)

1 - pbinom(S$SINE_Count, S$Deep_Mut_Count, D3$SINE_Ratio)
1 - pbinom(S$SINE_Sense_Count, S$Deep_Mut_Count, D3$SINE_Sense_Ratio)
1 - pbinom(S$SINE_Antisense_Count, S$Deep_Mut_Count, D3$SINE_Antisense_Ratio)




 
 
