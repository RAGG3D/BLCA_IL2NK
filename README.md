# BLCA_IL2NK
This is Figure 1
#FIG 1
x <- Bar_TCGAcelltype_BLCA(blcank) %>%
  mutate(nkstate = gsub("nk_resting", "ReNK", nkstate)) %>%  
  mutate(nkstate = gsub("nk_primed_IL2_PDGFD", "SPANK", nkstate)) %>%
  mutate(nkstate = gsub("nk_primed_IL2", "IL2NK", nkstate)) 
x$nkstate = factor(x$nkstate, levels=c("ReNK", "SPANK", "IL2NK"))
p1a01 <- x %>%
  filter(cat == "Fraction") %>%
  ggplot(aes(x = reorder(order, x = sample), y = profile, fill = nkstate)) +
  geom_bar(stat = "identity") +
  #facet_grid(~cat, scales = "free") +
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
  #facet_grid(~cat, scales = "free") +
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
