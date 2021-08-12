#FIG 4
x <- Gene_Cell("BLCA", blcank, "PDGFD", c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD')) %>%
  rbind(Gene_Cell("BLCA", blcank, "PDGFRB", c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD'))) %>%
  mutate(cat = gsub("nk_resting", "ReNK", cat)) %>%
  mutate(cat = gsub("nk_primed_IL2", "IL2NK", cat)) %>%
  mutate(cat = gsub("IL2NK_PDGFD", "SPANK", cat)) 
x <- as.data.frame(x)
x$item <- factor(x$item, levels=c("1/1", "1/2", "2/1", "2/2"), ordered=TRUE)
x <- x %>% 
  #filter(total_living_days >= 0) %>%
  filter( grepl('PDGFD/', cat))
x$cat <- factor(x$cat, levels=c("PDGFD/ReNK", "PDGFD/IL2NK", "PDGFD/SPANK"), ordered=TRUE)



p4 =  ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("cat"),
  conf.int = T,
  legend.title = "Gene Abundance/NK",
  legend.labs = list("L/L", "L/H", "H/L", "H/H"),
  conf.int.alpha = 0.15,
  risk.table = F,
  pval = T,
  short.panel.labs = T,
  linetype = "item",
  palette = "jco"
) + theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) 

ggsave(plot = p4, "output/BLCA-NEW-p4.pdf",device = "pdf", height = 4, width = 10)
