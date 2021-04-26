#MHC-I with NK phenotypes
#Each gene a facet
geneblca <- BLCA_transcript()

x <- geneblca %>% filter(symbol == "HLA-A"|symbol == "HLA-B"|symbol == "HLA-C"|symbol == "HLA-E"|
                      symbol == "B2M"|symbol == "TAP1"|symbol == "TAP2") %>%
  inner_join(blcacell %>% 
               filter(grepl("nk_", celltype))) %>%
  clinical_combine("BLCA") %>%
  filter(vital_status == 1) %>%
  group_by(celltype, symbol) %>%
  mutate(midx = median(raw_count_scaled), midy = median(fraction)) %>%
  ungroup %>%
  mutate(celltype = gsub("nk_resting", "ReNK", celltype)) %>%
  mutate(celltype = gsub("nk_primed_IL2_PDGFD", "SPANK", celltype)) %>%
  mutate(celltype = gsub("nk_primed_IL2", "IL2NK", celltype)) 

x$celltype <- factor(x$celltype, levels = c("ReNK", "IL2NK", "SPANK"))

each_gene <-  x %>% ggplot(aes(x = raw_count_scaled, y = fraction, color = total_living_days)) + 
  geom_jitter() +
  facet_grid(celltype~symbol) +
  scale_color_distiller(direction = 1, palette = "BuGn") +
  geom_hline(aes(yintercept = midy, group = celltype), colour = 'red') +
  geom_vline(aes(xintercept = midx, group = celltype), colour = 'red') +
  theme_test() +
  labs(x = "MHC-I genes", y = "NK Phenotypes (fraction)", color = "", title = "(Only dead samples)") +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom")
  
  ggsave(plot = each_gene, "output/Extra-in-each-MHC-I.png", height = 8, width = 15)
  
#All MHC-I together
  x <- geneblca %>% filter(symbol == "HLA-A"|symbol == "HLA-B"|symbol == "HLA-C"|symbol == "HLA-E"|
                             symbol == "B2M") %>%
    group_by(sample) %>%
    summarise(sum(raw_count_scaled)) %>%
    tidybulk::rename(MHCI = `sum(raw_count_scaled)`) %>%
    ungroup %>%
    inner_join(blcacell %>% 
                 filter(grepl("nk_", celltype))) %>%
    clinical_combine("BLCA") %>%
    filter(vital_status == 1) %>%
    group_by(celltype) %>%
    mutate(midx = median(MHCI), midy = median(fraction)) %>%
    ungroup %>%
    mutate(celltype = gsub("nk_resting", "ReNK", celltype)) %>%
    mutate(celltype = gsub("nk_primed_IL2_PDGFD", "SPANK", celltype)) %>%
    mutate(celltype = gsub("nk_primed_IL2", "IL2NK", celltype)) 
  
  x$celltype <- factor(x$celltype, levels = c("ReNK", "IL2NK", "SPANK"))
  
  all_MHCI <-  x %>% ggplot(aes(x = MHCI, y = fraction, color = total_living_days)) + 
    geom_jitter() +
    facet_grid(celltype~.) +
    scale_color_distiller(direction = 1, palette = "BuGn") +
    geom_hline(aes(yintercept = midy, group = celltype), colour = 'red') +
    geom_vline(aes(xintercept = midx, group = celltype), colour = 'red') +
    theme_test() +
    labs(x = "MHC-I", y = "NK Phenotypes (fraction)", color = "", title = "(Only dead samples)") +
    theme(text=element_text(size=16, family="sans"),
          legend.position = "bottom")
  
  ggsave(plot = all_MHCI, "output/Extra-all-MHC-I.png", height = 10, width = 6)
