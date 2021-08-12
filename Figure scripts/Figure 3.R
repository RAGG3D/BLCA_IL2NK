#FIG 3
x <- PDGFsurvival("BLCA", 2) %>%
  gather(cat, item, -sample) %>% 
  TCGA_clinical("BLCA")
x <- as.data.frame(x)
x$cat <- factor(x$cat, levels=c("PDGFD", "PDGFRB"), ordered=TRUE)
x$item <- factor(x$item, levels=c(1, 2), ordered=TRUE)
p3 =  ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("cat"),
  conf.int = T,
  risk.table = F,
  conf.int.alpha = 0.15,
  pval = T,
  linetype = "item",
  legend.title = "Abundance ",
  legend.labs = list("L ", "H "),
  short.panel.labs = T,
  palette = "jco"
) + theme_bw() + 
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) +
  theme(aspect.ratio=1)

ggsave(plot = p3, "output/BLCA-NEW-p3.pdf", device = "pdf", height = 4, width = 7)

