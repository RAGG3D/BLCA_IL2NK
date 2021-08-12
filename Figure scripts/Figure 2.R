#FIG 2
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
p2 <- ggboxplot(x, x = "nkstate", y = "profile", facet.by = "histologic_grade", fill = "nkstate", nrow = 1) + 
  scale_fill_brewer(name = "", palette="RdYlGn", labels = c("ReNK", "IL2NK", "SPANK"), direction = -1) +
  labs(x = "", y = "Fraction") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  
  geom_jitter(aes(color = nkstate), alpha = 0.6, size = 0.3) +
  scale_color_manual(name = "", values = c("#3b5629", "#e6b800", "#995c00"), labels = c("ReNK", "IL2NK", "SPANK")) +
  theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom")

ggsave(plot = p2, "output/BLCA-NEW-p2.pdf", device = "pdf", height = 5, width = 8)





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
  theme(text=element_text(size=16, family="sans"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom")

ggsave(plot = pd, "output/BLCA-NEW-supp1.pdf", device = "pdf", height = 5, width = 8)
