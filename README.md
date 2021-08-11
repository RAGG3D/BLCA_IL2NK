# BLCA_IL2NK

## Packages and source required
```
options(connectionObserver = NULL) #This is to help load org.Hs.eg.db 
```
```{r,}
library(tidyverse)
library(tidybulk)
library(survminer)
library(survival)
library(gapminder)
library(foreach)
library(org.Hs.eg.db)
library(cowplot)
library(ggsci)
library(GGally)
library(gridExtra)
library(grid)
library(reshape)
library(Hmisc)
library(viridis)
library(furrr)

source("src/functions.R") #Source all the functions from src
```

## Prepare data
Firstly we need to run the signature of selected cell types with patient RNA-seq dataset in CIBERSORT:
```
res = res_run(readRDS("data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)])
blcank <- foreach(i = as.character(res$stratification_cellularity[[1]]$.cell_type), .combine = bind_rows) %do% {
  x <- res$stratification_cellularity[[1]] %>%
    filter(.cell_type == i)
  
  x$cell_type_proportions[[1]]  %>%
    mutate(celltype = i, fraction = .proportion, nkstate = celltype, sample = patient)%>%
    dplyr::select(sample, nkstate, fraction)
} %>%    filter(grepl("nk", nkstate)) 
```
Then we formatted the matrix of cell type fractions to make it easier to use in the following analysis:
```
blcank <- foreach(i = as.character(res$stratification_cellularity[[1]]$.cell_type), .combine = bind_rows) %do% {
  x <- res$stratification_cellularity[[1]] %>%
    filter(.cell_type == i)
  
  x$cell_type_proportions[[1]]  %>%
    mutate(celltype = i, fraction = .proportion, nkstate = celltype, sample = patient)%>%
    dplyr::select(sample, nkstate, fraction)
} %>%    
    filter(grepl("nk", nkstate)) 

blcacell <- foreach(i = as.character(res$stratification_cellularity[[1]]$.cell_type), .combine = bind_rows) %do% {
  x <- res$stratification_cellularity[[1]] %>%
    filter(.cell_type == i)
  
  x$cell_type_proportions[[1]]  %>%
    mutate(celltype = i, fraction = .proportion, sample = patient)%>%
    dplyr::select(sample, celltype, fraction)
}
```

After preparing the data we can start making figures.

## Figure 1
### FIG 1A

```
x <- Bar_TCGAcelltype_BLCA(blcank) %>%
  mutate(nkstate = gsub("nk_resting", "ReNK", nkstate)) %>%  
  mutate(nkstate = gsub("nk_primed_IL2_PDGFD", "SPANK", nkstate)) %>%
  mutate(nkstate = gsub("nk_primed_IL2", "IL2NK", nkstate)) 
  
x$nkstate = factor(x$nkstate, levels=c("ReNK", "SPANK", "IL2NK"))

p1a01 <- x %>%
  filter(cat == "Fraction") %>%
  ggplot(aes(x = reorder(order, x = sample), y = profile, fill = nkstate)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#e7e7e7",  "#cd534c", "#7aa6dc")) +
  coord_flip() + 
  theme_classic() +
  #theme(panel.background=element_rect(fill='transparent',color ="gray")) +
  labs( x = "Patients", y = "Fraction", tag = "A") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text=element_text(size=16, family="sans"),
        legend.position = "bottom")
        
p1a02 <- x %>%
  filter(cat == "Percentage") %>%
  ggplot(aes(x = reorder(order, x = sample), y = profile, fill = nkstate)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#e7e7e7",  "#cd534c", "#7aa6dc")) +
  coord_flip() + 
  theme_classic() +
  theme(panel.background=element_rect(fill='transparent',color ="gray")) +
  labs(x = "Patients", y = "Percentage", tag = "") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        text=element_text(size=16, family="sans"),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") + 
  guides(fill=guide_legend(title="NK Phenotype"))
```

### FIG 1B
```
x <- blcank %>%
  spread(nkstate, fraction) %>%
  mutate(nk_primed_IL2 = ifelse(
    nk_primed_IL2 > median(nk_primed_IL2), "H", "L"
  )) %>%
  mutate(nk_primed_IL2_PDGFD = ifelse(
    nk_primed_IL2_PDGFD > median(nk_primed_IL2_PDGFD), "H", "L"
  )) %>% 
  mutate(nk_resting = ifelse(
    nk_resting > median(nk_resting), "H", "L"
  )) %>%
  gather(celltype, fraction, -sample) %>%
  inner_join(read.csv("data/clinical_blca.csv"), by = c("sample" = "bcr_patient_barcode"))%>%
  mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
  mutate(na = is.na(total_living_days)) %>%
  mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
  dplyr::select(sample, celltype, fraction, vital_status, total_living_days, age) %>%
  mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))
  
x <- as.data.frame(x)

x$celltype = factor(x$celltype, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD'))

x$fraction = factor(x$fraction, levels=c("L", "H"))

p1b =  ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ fraction,
    data = x 
  ),
  data = x,
  facet.by = c("celltype"),
  conf.int = T,
  risk.table = F,
  conf.int.alpha = 0.15,
  pval = T,
  legend.title = "Cell Fraction",
  legend.labs = list("L ", "H "),
  short.panel.labs = T,
  panel.labs = list(celltype = c("ReNK", "IL2NK", "SPANK")),
  linetype = "fraction",
  palette = "jco"
) + theme_bw() + 
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) +
  labs(tag = "B")
```

Then we put the panels together:
```
plot_grid(plot_grid(p1a01, p1a02, nrow = 1), p1b, ncol = 1)
```
All the figure should be saved as a PDF for further edits.
(After modification)
![image](https://github.com/RAGG3D/BLCA_IL2NK/blob/master/figures/BLCA-NEW-p1.jpg)

## Figure 2
### FIG 2A
```
x <- Bar_TCGAcelltype_BLCA(blcank) %>%
  inner_join(read.csv("data/clinical_blca.csv"), by = c("sample" = "bcr_patient_barcode"))%>%
  mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
  mutate(na = is.na(total_living_days)) %>%
  mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
  dplyr::select(sample, nkstate, profile, cat, histologic_grade, vital_status, total_living_days, age) %>%
  mutate(vital_status = ifelse(vital_status == "Dead", 1, 0)) %>%
  filter(cat == "Fraction") %>%
  rbind(Bar_TCGAcelltype_BLCA(blcank) %>%
          inner_join(read.csv("data/clinical_blca.csv"), by = c("sample" = "bcr_patient_barcode"))%>%
          mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
          mutate(na = is.na(total_living_days)) %>%
          mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
          dplyr::select(sample, nkstate, profile, cat, histologic_grade, vital_status, total_living_days, age) %>%
          mutate(vital_status = ifelse(vital_status == "Dead", 1, 0)) %>%
          mutate(histologic_grade = "All") %>%
          filter(cat == "Fraction")) %>%
  filter(histologic_grade != "[Unknown]")
  
x <- as.data.frame(x)

x$nkstate = factor(x$nkstate, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD'))

x$histologic_grade = factor(x$histologic_grade, levels=c('All','Low Grade', 'High Grade'))

my_comparisons <- list( c("nk_resting", "nk_primed_IL2"), c("nk_primed_IL2", "nk_primed_IL2_PDGFD"), c("nk_resting", "nk_primed_IL2_PDGFD") )

pc <- ggboxplot(x, x = "nkstate", y = "profile", facet.by = "histologic_grade", fill = "nkstate", nrow = 1) + 
  scale_fill_brewer(name = "", palette="RdYlGn", labels = c("ReNK", "IL2NK", "SPANK"), direction = -1) +
  labs(x = "", y = "Fraction", tag = "A") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  
  geom_jitter(aes(color = nkstate), alpha = 0.6, size = 0.3) +
  scale_color_manual(name = "", values = c("#3b5629", "#e6b800", "#995c00"), labels = c("ReNK", "IL2NK", "SPANK")) +
  theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom")
 ```

### FIG 2B
```
lab = as.character(as.data.frame(x %>% group_by(histologic_grade) %>% summarise())[,1])
my_comparisons <- list( c(lab[1], lab[2]), c(lab[1], lab[3]),c(lab[2], lab[3]))
pd <- ggboxplot(x, x = "histologic_grade", y = "profile", facet.by = "nkstate", nrow = 1, fill = "histologic_grade",
                panel.labs = list(nkstate = c("ReNK", "IL2NK", "SPANK"))) + 
  scale_fill_brewer(name = "", palette="Set3", labels = lab, direction = 1) +
  labs(x = "", y = "Fraction") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  geom_jitter(aes(color = histologic_grade), alpha = 0.6, size = 0.3) +
  scale_color_manual(name = "", values = c("#328173", "#ff9900", "#5a5095"), labels = lab) +
  theme_bw() +
  labs(tag = "")+
  theme(text=element_text(size=16, family="sans"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom")
```
Then we need to format the 2 panels:
```
plot_grid(plot_grid(pc, pd, nrow = 1), 
plot_grid(NK_KM_Grade("Low Grade") + labs(tag = "B"), NK_KM_Grade("High Grade")+labs(tag = ""), ncol = 1), 
nrow = 2)
```
(After modification)
![image](https://github.com/RAGG3D/BLCA_IL2NK/blob/master/figures/BLCA-NEW-p2.jpg)

## Figure 3
### FIG 3
```
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
```
(After modification)
![image](https://github.com/RAGG3D/BLCA_IL2NK/blob/master/figures/BLCA-NEW-p3.jpg)

## Figure 4
### FIG  4
```
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
```
(After modification)
![image](https://github.com/RAGG3D/BLCA_IL2NK/blob/master/figures/BLCA-NEW-p4.jpg)

## Figure 5
### FIG 5A
```
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
```

### FIG 5B
```
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
```

### FIG 5C
```
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

a <- melt(cor(x %>% dplyr::select(-sample))[-c(1:11), c(1:11)])
a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
a <- a %>% inner_join(as.data.frame(a.p$P[-c(1:11), c(1:11)]) %>% 
                        mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:11), c(1:11)]))) %>% 
                        gather(X2, p.value, -X1)) 
a$stars <- cut(a$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
a$X2 <- factor(a$X2, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD', "t_helper", "t_CD8_naive", "t_gamma_delta", 't_CD4_memory_central', 
                              't_CD4_memory_effector', "t_CD8_memory_central", "t_CD8_memory_effector", "t_reg"), ordered=TRUE)
a$X1 <- factor(a$X1, levels = c("LTA", "BTLA", "TNFSF14", "CD160", "TNFRSF14", "KLRK1"), ordered = T)

p5c <- ggplot(data = a , aes(x=`X2`, y=`X1`, fill=value)) + 
  geom_tile() + 
  geom_text(aes(label=stars)) + 
  scale_fill_viridis() +
  labs(y = "", x = "") + #y = "Log Transformed Normalized Count", x = "Log Transformed NK Proportion"
  theme_bw() +
  scale_x_discrete(labels = c("ReNK", "IL2NK", "SPANK", "Helper\nT", "Naive\nCD8 T", "?æ? T", "Memory\nCD4 T CTL", 
                              "Memory\nCD4 T EFC", "Memory\nCD8 T CTL", "Memory\nCD8 T EFC", "Reg T")) +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -45),
        legend.key.width = unit(3, "line")) +
  labs(fill = "Correlation")
```

Then we need to format the 3 panels:
```
plot_grid(plot_grid(p5a, p5b, nrow = 1), p5c, ncol = 1)
```
(After modification)
![image](https://github.com/RAGG3D/BLCA_IL2NK/blob/master/figures/BLCA-NEW-p5.jpg)
