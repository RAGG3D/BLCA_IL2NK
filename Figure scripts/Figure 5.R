#FIG 5
x <- genelist_cancer(BLCA, "BLCA", list("KLRK1")) %>%
  filter(cat == "median")
x$symbol <- factor(x$symbol, levels=c("KLRK1"), ordered=TRUE)
p5a <- ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("symbol"),
  conf.int = T,
  risk.table = F,
  legend.labs = list("L", "H", "?", "?"),
  short.panel.labs = T,
  pval = T,
  nrow =2,
  conf.int.alpha = 0.15,
  legend.title = "Abundance",
  palette = "jco"
) + theme_bw() +
  labs(tag = "A") +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") 

cancer = "BLCA"
x <- foreach(i = "CD160", .combine = bind_rows) %do% {
  BLCA %>%
    filter(symbol == "TNFRSF14" | symbol == i) %>%
    dplyr::select(sample, symbol, raw_count_scaled) %>%
    spread(symbol, raw_count_scaled) %>%
    mutate(gene1 = factor(Hmisc::cut2(!!as.name("TNFRSF14"), g = 2), labels = c(1:2))) %>%
    mutate(gene2 = factor(Hmisc::cut2(!!as.name(i), g = 2), labels = c(1:2))) %>%
    unite("item", gene1, gene2, sep = "/", remove = T) %>%
    mutate(gene1 = "TNFRSF14", gene2 = i) %>%
    dplyr::select(sample, item, gene1, gene2)
} %>% 
  inner_join(read.csv(paste0("data/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
  mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
  mutate(na = is.na(total_living_days)) %>%
  mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
  dplyr::select(sample, gene1, gene2, item, vital_status, total_living_days, age) %>%
  mutate(vital_status = ifelse(vital_status == "Dead", 1, 0)) %>%
  mutate(item = ifelse(item == "2/2", "H/H", "other groups"))
x <- as.data.frame(x)
x$item = factor(x$item, levels=c("other groups", "H/H"))
p5b <- ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("gene1", "gene2"),
  conf.int = T,
  risk.table = F,
  conf.int.alpha = 0.15,
  pval = T,
  legend.title = "Row/Column",
  legend.labs = list("other groups", "H/H"),
  short.panel.labs = T,
  linetype = "item",
  palette = c("#01665E","#C71000FF")
) + theme_bw() + 
  labs(tag = "B") +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE)


Immunereceptor <- c("KLRK1", "TNFRSF14", "CD160", "TNFSF14" , "LTA", "BTLA")

x <- blcacell %>%
  spread(celltype, fraction) %>%
  dplyr::select(sample, contains("nk_"), contains("t_"), -mast_cell) %>%
  inner_join(foreach(i = Immunereceptor, .combine = bind_rows) %do%{
    BLCA %>% 
      filter(symbol == i) %>%
      mutate(item = log(raw_count_scaled + 1), cat = symbol) 
  } %>%
    dplyr::select(sample, cat, item) %>%
    spread(cat, item))

a <- reshape::melt(cor(x %>% dplyr::select(-sample))[-c(1:11), c(1:11)])
a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
a <- a %>% inner_join(as.data.frame(a.p$P[-c(1:11), c(1:11)]) %>% 
                        mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:11), c(1:11)]))) %>% 
                        gather(X2, p.value, -X1)) %>%
  inner_join(read_csv("data/cellname.csv"))
a$`Cell Type` <- factor(a$`Cell Type`, levels=c('ReNK', 'IL2NK', 'SPANK', "Helper T", "Naive CD8 T", "GD T", 'Memory CD4 T CTL', 
                                                'Memory CD4 T EFC', "Memory CD8 T CTL", "Memory CD8 T EFC", "Reg T"), ordered=TRUE)
a$X1 <- factor(a$X1, levels = c("LTA", "BTLA", "TNFSF14", "CD160", "TNFRSF14", "KLRK1"), ordered = T)
p5c <- as_tibble(a ) %>%
  tidybulk::rename("Receptor" = "X1") %>%
  tidyHeatmap::heatmap(
    Receptor,
    `Cell Type`,
    value,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    palette_value = circlize::colorRamp2(
      
      seq(2,-2, length.out = 11),
      
      RColorBrewer::brewer.pal(11, "RdBu")
      
    )
    
  ) 

ggsave("output/BLCA-NEW-p5-heat.pdf",device = "pdf", height =6, width = 15)

