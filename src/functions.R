res_run <- function(x) {
  
  tibble(
    
    file =
      
      dir(
        
        pattern = "armet_",
        
        "D:/PhD Projects/Rprojects/BLCA_IL2NK/data/",
        
        full.names = T
        
      ) %>%
      
      grep("input.rds", ., value = T)
    
  ) %>%
    
    extract(file, "cancer_type", ".*armet_([A-Z]+)_input.*", remove = FALSE) %>%
    
    mutate(
      
      stratification_cellularity =
        
        furrr::future_map(
          
          file, ~ readRDS(.x) %>%
            
            test_stratification_cellularity(
              
              survival::Surv(PFI.time.2, PFI.2) ~ .,
              
              patient,
              
              transcript,
              
              count,
              
              cores = 1,
              
              reference = x 
              
            )
          
        )
      
    )
}

Bar_TCGAcelltype_BLCA <- function(a){
  x <- a %>% 
    group_by(sample) %>% 
    summarise(order = sum(fraction)) %>%
    right_join(a, by = c("sample")) %>%
    tidybulk::rename(profile = fraction) %>%
    mutate(cat = "Fraction") %>%
    rbind(a %>% 
            group_by(sample) %>% 
            summarise(sum = sum(fraction)) %>%
            right_join(a) %>%
            filter(sum > 0) %>%
            mutate(profile = fraction/sum) %>%
            left_join(a %>%
                        group_by(sample) %>% 
                        summarise(sum = sum(fraction)) %>%
                        right_join(a) %>%
                        filter(sum > 0) %>%
                        mutate(profile = fraction/sum) %>%
                        filter(nkstate == "nk_primed_IL2") %>%
                        mutate(order = profile) %>%
                        dplyr::select(sample, order)) %>%
            dplyr::select(-sum, -fraction) %>%
            mutate(cat = "Percentage"))
} 

Gene_Cell <- function(cancer, a, gene, cells){
  foreach(i = cells, .combine = bind_rows) %do%{
    BLCA %>%
      filter(symbol == gene) %>%
      dplyr::select(sample, symbol, raw_count_scaled) %>%
      spread(symbol, raw_count_scaled) %>%
      mutate(gene = factor(Hmisc::cut2(!!as.name(gene), g = 2), labels = c(1:2))) %>%
      inner_join(a %>%
                   spread(nkstate, fraction) %>%
                   mutate(cell = factor(Hmisc::cut2(!!as.name(i), g = 2), labels = c(1:2))) %>%
                   dplyr::select(sample, cell)) %>%
      unite("item", gene, cell, sep = "/", remove = T) %>%
      dplyr::select(sample, item) %>% 
      mutate(cat = paste0(gene, "/", i))
  } %>% TCGA_clinical(cancer)
}

NK_KM_Grade <- function(subtype) {
  
  x <- Bar_TCGAcelltype_BLCA(blcank) %>%
    inner_join(read.csv("data/clinical_blca.csv"), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    dplyr::select(sample, nkstate, profile, order, cat, histologic_grade, vital_status, total_living_days, age) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0)) %>%
    filter(histologic_grade == subtype) %>%
    mutate(nkstate = gsub("nk_resting", "ReNK", nkstate)) %>%  
    mutate(nkstate = gsub("nk_primed_IL2_PDGFD", "SPANK", nkstate)) %>%
    mutate(nkstate = gsub("nk_primed_IL2", "IL2NK", nkstate)) 
  
  x$nkstate = factor(x$nkstate, levels=c("ReNK", "SPANK","IL2NK"))
  p2a01 <- x %>%
    filter(cat == "Fraction") %>%
    ggplot(aes(x = reorder(order, x = sample), y = profile, fill = nkstate)) +
    geom_bar(stat = "identity") +
    #facet_grid(~cat, scales = "free") +
    scale_fill_manual(values = c("#deebf7", "#3182bd", "#9ecae1")) +
    coord_flip() + 
    theme_bw() +
    theme(panel.background=element_rect(fill='transparent',color ="gray")) +
    labs( x = subtype, y = "Fraction") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text=element_text(size=16, family="sans"),
          legend.position = "bottom")
  p2a02 <- x %>%
    filter(cat == "Percentage") %>%
    ggplot(aes(x = reorder(order, x = sample), y = profile, fill = nkstate)) +
    geom_bar(stat = "identity") +
    #facet_grid(~cat, scales = "free") +
    scale_fill_manual(values = c("#deebf7", "#3182bd", "#9ecae1")) +
    coord_flip() + 
    theme_bw() +
    theme(panel.background=element_rect(fill='transparent',color ="gray")) +
    labs(x = subtype, y = "Percentage") +
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
    dplyr::select(sample, celltype, fraction, histologic_grade, vital_status, total_living_days, age) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))%>%
    filter(histologic_grade == subtype)
  
  x <- as.data.frame(x)
  x$celltype = factor(x$celltype, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD'))
  x$fraction = factor(x$fraction, levels=c("L", "H"))
  p2b =  ggsurvplot(
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
    guides(linetype = FALSE) 
  
  plot_grid(plot_grid(p2a01, p2a02, nrow = 1), p2b, nrow = 1)}

PDGFsurvival <- function(cancer, n) {
  x <- BLCA %>% 
    filter(symbol == "PDGFD"|symbol == "PDGFRB") %>%
    dplyr::select(sample, symbol, raw_count_scaled) %>%
    spread(symbol, raw_count_scaled) %>%
    mutate(PDGFD = factor(Hmisc::cut2(PDGFD, g = n), labels = c(1:nlevels(Hmisc::cut2(PDGFD, g = n))))) %>%
    mutate(PDGFRB = factor(Hmisc::cut2(PDGFRB, g = n), labels = c(1:nlevels(Hmisc::cut2(PDGFRB, g = n))))) }

TCGA_clinical <- function(x, cancer) {
  x %>% 
    inner_join(read.csv(paste0("data/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    dplyr::select(sample, cat, item, vital_status, total_living_days, age) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))}

genelist_cancer <- function(y, cancer, genelist) {
  x <- foreach(i = genelist, .combine = bind_rows) %do% {
    y %>% filter(symbol == i) %>%
      mutate(median = factor(Hmisc::cut2(raw_count_scaled, g = 2), labels = c(1:2))) %>%
      mutate(quantile = factor(Hmisc::cut2(raw_count_scaled, g = 4), labels = c(1:nlevels(Hmisc::cut2(raw_count_scaled, g = 4))))) 
  } %>% 
    dplyr::select(sample, symbol, median, quantile) %>%
    gather(cat, item, -c(sample, symbol)) %>%
    inner_join(read.csv(paste0("data/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0)) %>%
    dplyr::select(sample, symbol, cat, item, vital_status, total_living_days, age)
  
  x <- as.data.frame(x)
  x$item = factor(x$item, levels=c(1,2,3,4))
  x$cat <- factor(x$cat, levels=c("median", "quantile"), ordered=TRUE)
  x
}

BLCA_transcript <- function(){
  foreach(i = list.files("data/Transcript/"), .combine = bind_rows) %do% {
    read_table2(paste0("data/Transcript/", i), col_names = F) %>%
      mutate(sample = i)
  }  %>% 
    tidybulk::rename(ensembl = `X1`) %>%
    separate(ensembl, c("ensembl_id", "c"), sep = "\\.") %>%
    inner_join(toTable(org.Hs.egENSEMBL)) %>%
    inner_join(toTable(org.Hs.egSYMBOL)) %>%
    dplyr::select(sample, symbol, raw_count) %>%
    mutate(sample = paste0(substr(basename(as.character(sample)), start = 1, stop = 49), ".gz")) %>%
    as.tibble() %>%
    inner_join(
      read_csv(paste0("data/gdc_sample_sheet.csv"))%>% mutate(sample = `File Name`)) %>%
    mutate(sample = `Case ID`) %>%
    dplyr::select(sample, symbol, raw_count) %>%
    tidybulk::aggregate_duplicates(.sample = sample, .abundance = raw_count, .transcript = symbol) %>%
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol) 
}
