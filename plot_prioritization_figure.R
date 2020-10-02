library(data.table)
library(tidyverse)
library(ggplot2)
library(stringr)

dat <- fread("CELLECT-LDSC/results/prioritization.csv")

bmi <- dat %>% filter(gwas == "BMI_UKBB_Loh2018") %>% filter(specificity_id == "r_adan_lepR")

bmi %>% filter(annotation != "nan") %>%
  mutate(taxonomy = case_when(grepl("neurons",annotation) ~ "Neurons",
                              grepl("Tanycytes",annotation) ~ "Tanycytes",
                              grepl("Vlmcs",annotation) ~ "VLMCs",
                              grepl("endothelial",annotation) ~ "Endothelial",
                              T ~ annotation),
               bonf_significant = if_else(pvalue <= 0.05/nrow(bmi), true=T, false=F),
               p.value.mlog10 = -log10(pvalue)) %>%
  mutate(annotation = str_replace_all(annotation, "_", " "),
         taxonomy = str_replace_all(taxonomy, "_", " "),
         taxonomy = ifelse(str_detect(taxonomy, "[[:upper:]]"), taxonomy, str_to_title(taxonomy))) %>%
  mutate(taxonomy = factor(taxonomy, levels = rev(c("Neurons","Pitituary","Astrocytes","OPCs","Tanycytes","Oligodendrocytes",
                                              "VLMCs","Microglia","Endothelial","Nrn1 non neuronal")))) -> p.bmi 

p.bmi %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(taxonomy, as.character(annotation))]))) -> p.bmi

p.bmi$annotation %>% levels()
p.bmi$taxonomy %>% levels()


df.tax_text_position <- p.bmi %>% 
  group_by(taxonomy) %>% 
  summarize(n_annotations_in_tax=n()) %>% 
  arrange(taxonomy) %>%
  mutate(pos_mid=cumsum(n_annotations_in_tax)-n_annotations_in_tax/2,
         pos_start=cumsum(n_annotations_in_tax)-n_annotations_in_tax,
         pos_end=cumsum(n_annotations_in_tax),
         idx=1:n()) %>% # idx is used for identifying every n'th tax
  mutate(flag_draw_rect=if_else(idx %% 2 == 0, F, T)) %>% 
  arrange(n_annotations_in_tax) # Remove some text because they contain too few cell-types


plotting_func <- function(df.plot, df.tax_text_position){
  p.main <- ggplot() +
    ### add ggplot 'baselayer'. This makes our 'canvas' and needs to go first (for some reason I don't fully understand...)
    geom_blank(data=df.plot, aes(x=annotation, y=p.value.mlog10)) +
    ### add tax info and gray rects (add gray rect before the data points, so avoid them 'covering' the points)
    geom_text(data=df.tax_text_position, aes(x=pos_mid+0.5, y=-0.5, label=taxonomy), hjust="right", nudge_y = -2, size=rel(6)) +
    geom_rect(data=df.tax_text_position %>% filter(flag_draw_rect), aes(xmin=pos_start+0.5, xmax=pos_end+0.5, ymin=-5, ymax=Inf), color="gray", alpha=0.1) +
    ### cell-types
    geom_hline(yintercept=-log10(0.05/nrow(df.plot)), linetype="dashed", color="darkgray", size = 1) + 
    ### axes
    labs(x="", y=expression(-log[10](P[S-LDSC]))) +
    # coord
    coord_flip(ylim = c( 0, max(df.plot$p.value.mlog10) ), # This focuses the y-axis on the range of interest
               clip = 'off') +   # This keeps the labels from disappearing 
    theme_classic() + 
    theme(axis.text.y=element_text(size=rel(1.2)),
          axis.ticks.y=element_blank())
  p.main <- p.main + theme(plot.margin = unit(c(1,1,1,6), "cm")) # (t, r, b, l) widen left margin
  return(p.main)
  }


p.main <- plotting_func(p.bmi, df.tax_text_position)
p.main <- p.main + geom_point(data=p.bmi, aes(x=annotation, y=-log10(pvalue)), color="gray", size=5)
p.main <- p.main + geom_point(data=p.bmi %>% filter(bonf_significant), aes(x=annotation, y=-log10(pvalue), color=annotation), size=5)
p.main <- p.main + theme(legend.position="bottom", legend.text = element_text(size = 15),
                         axis.text.x = element_text(size = 15),
                         axis.title.x = element_text(size = 15))
p.main <- p.main + ggrepel::geom_text_repel(data=p.bmi %>% filter(bonf_significant), aes(x=annotation, y=-log10(pvalue), label=annotation, color=annotation), size = rel(7),
                                            hjust = 0, nudge_y = 1, show.legend=F, min.segment.length = 0, box.padding = 0.8, direction = "y") 

p.main <- p.main + ggtitle("CELLECT run on BMI gwas") +
  theme(plot.title = element_text(size=20, hjust = -0.85))
p.main

ggsave(p.main, filename="./CELLECT_bmi_loh_fig.pdf", width=12, height=9)



