
# packages ----------------------------------------------------------------

syn.client <- reticulate::import('synapseclient')
syn.utils <- reticulate::import('synapseutils')
syn <- syn.client$Synapse()

library(callr)
library(furrr)
library(clusterProfiler)
library(tidyverse)

# proteomics data ---------------------------------------------------------

syn$login()

# metadata synID
synid <- 'syn65986467'

# metadata
m.files = readxl::read_xlsx(syn$get(synid) %>% .$path,sheet = 2) %>% 
  mutate(sample = str_remove(SampleName, '-[1-9]$'),
         rep = str_extract(SampleName, '[1-9]$'),
         sample = if_else(sample == 'KD-sicon', 'KD-siCon', sample)
  ) %>% 
  mutate(name = str_remove(Name, '^ranjita_dia__') %>% 
           str_remove(., '.d.rar$'))

m.samples = 
  bind_rows(
    readxl::read_xlsx(
      syn$get(synid) %>% .$path,
      sheet = 1,
      range = 'L2:N32', 
      col_names = c('sample', 'target','passage_n')) %>% 
      mutate(cell = 'scramble'),
    readxl::read_xlsx(
      syn$get(synid) %>% .$path,
      sheet = 1,
      range = 'P2:R32', 
      col_names = c('sample', 'target','passage_n')) %>% 
      mutate(cell = 'psen2_kd')
  ) 

raw.m = left_join(m.files, m.samples, by = 'sample') 

# DIA-NN normalized data
raw.d <- read_tsv( syn$get('syn65986474') %>% .$path )

n.pep <- read_tsv( syn$get('syn65986477') %>% .$path ) %>% 
  summarise(n.pep = length(unique(Precursor.Id)), .by = c(Genes)) 


# align data & metadata  --------------------------------------------------

d <- raw.d %>% 
  select(-c( Knockdown_0_Scrambled_1, KD_SCR_Gene,Group) ) %>% 
  column_to_rownames(var = 'Sample.Name') %>% 
  mutate(across(everything(), ~ if_else( .x == 0, NA_real_, .x))) %>% 
  as.matrix() %>% t()

dups <- raw.m$Name[which(duplicated(raw.m$name))]

m <- raw.m[match( colnames(d), raw.m$name ),] 

# update metadata ---------------------------------------------------------

s.count <- map_dbl(1:ncol(d),  ~ sum(!is.na( d[,.x] )))
p.count <- map_dbl(1:nrow(d),  ~ sum(!is.na( d[.x,] )))
s.med <- map_dbl(1:ncol(d), ~ median(d[,.x], na.rm = T))

m <- bind_cols(m, n_prot = s.count, med_expr = s.med) 

# filtering & QC ----------------------------------------------------------

# samples to remove
gis_n <- which(m$SampleName == 'GIS')
outlier_n <- which(m$n_prot < 6100)

# data table
data <- bind_cols(m, t(d)) %>% slice(-c(gis_n,outlier_n))

# number of samples per group where protein is present
n.na = data %>% group_by(target, cell) %>% 
  summarise( across(Tap1:(ncol(data)-2), ~ sum( !is.na(.x) )) )

# which proteins are ID'd in at least 2 samples in at least 2 groups 
idx = map_lgl(3:ncol(n.na), ~ n.na[,.x] %>% .[.>=2] %>% length() >= 2)

# data table
data <- bind_cols(m, t(d[idx,])) %>% slice(-c(gis_n,outlier_n))

# differential expression -------------------------------------------------

fallbackIfSmallTukeyP <- TRUE

# remove metadata columns and log2 transform
dataclipped <- data[, which((names(data) %in% names(m))==F) ] %>% 
  mutate(across(everything(), ~log2(.x+1))) %>% 
  as.matrix()

# set up sample type and batch
SampleType <- str_c(data$target,':',data$cell) 
Batch <- as.character(data$rep)

# set up data for DE testing
aov <- aov(as.double(rnorm(n=nrow(data)))~SampleType+Batch, data=data)
tuk <- TukeyHSD(aov)
comparisonList <- rownames(tuk$SampleType)

# run DE analysis
start = Sys.time()
parallel_de <- function(dataclipped, SampleType, Batch, comparisonList) {
  
  # set up parallel environment
  library(furrr)
  library(tidyverse)
  plan('multisession', workers = parallel::detectCores() - 2)
  
  # calculate DE proteins using anova in parallel
  de.results <- furrr::future_map_dfr(
    
    # for all protein columns in data
    1:ncol(dataclipped), 
    
    ~ {
      # data for protein X
      x = as.double( dataclipped[,.x] )
      
      # perform ANOVA and Tukey HSD test
      aov <- aov( x ~ SampleType + Batch )
      anovaresult <- anova(aov)
      tuk <- TukeyHSD(aov)
      
      # collect comparisons
      local.comp <- rownames(tuk$SampleType)
      tukresult <- as_tibble(tuk$SampleType) %>% 
        rename_with(~str_replace_all(.x,' ','.')) %>% 
        mutate(comp = local.comp)
      
      # match up contrasts
      if( length(local.comp) != length(comparisonList) ) tukresult <-
        tukresult[match(comparisonList, local.comp), ]
      
      # count small p-values to re-compute with t-test
      zeroFallbacks=length(which(tukresult[,"p.adj"]<10^-8.5))
      
      # re-compute p-values using t-test for small p from Tukey HSD
      if (zeroFallbacks>0) for (comp in which(tukresult[,"p.adj"]<10^-8.5)) {
        this.comp=local.comp[comp]
        grp1=gsub("^(.*)\\-.*","\\1",this.comp)
        grp2=gsub("^.*\\-(.*)$","\\1",this.comp)
        
        g1.idx = which(SampleType == grp1)
        g2.idx = which(SampleType == grp2)
        
        lowTuk.p.estimate <- tryCatch( 
          p.adjust(
            t.test(
              x[g1.idx],
              x[g2.idx],
              alternative="two.sided",
              var.equal=FALSE)$p.value, 
            method="bonferroni",
            n=nrow(tukresult)), error=function(e) c(1) )
        
        tukresult[comp,"p.adj"] <- lowTuk.p.estimate
      }
      
      # join results into data frame row
      tibble(
        gene = colnames(dataclipped)[.x],
        anova.f = anovaresult$`F value`[1], 
        anova.p = anovaresult$`Pr(>F)`[1], 
        n.zero.pval = zeroFallbacks,
        tukey = list(tukresult)
      )
    }
  )
  
  # Save results to a file
  saveRDS(de.results, file = here::here('manuscript','results',"de_prot_results.rds"))
}

# Run the parallel computation in the background
bg_process <- callr::r_bg(func = parallel_de,
                          args = list(dataclipped, SampleType, Batch, comparisonList))

# Check if the background process is still running
bg_process$is_alive()

# read in DEP results -----------------------------------------------------

# read in the results:
de.results <- readRDS(here::here('manuscript','results',"de_prot_results.rds")) %>% 
  unnest(tukey) %>% 
  filter(!is.na(comp)) %>% 
  
  # - split the comparison column 
  mutate(
    diff.og = diff,
    tg.1 = str_split_fixed(comp, ':|-',4)[,1],
    cell.1 = str_split_fixed(comp, ':|-',4)[,2],
    tg.2 = str_split_fixed(comp, ':|-',4)[,3],
    cell.2 = str_split_fixed(comp, ':|-',4)[,4]
  ) %>% 
  
  # - reverse the directionality for comparisons that have the control first
  mutate(tg = if_else(grepl('siRNA',tg.1), tg.2,tg.1),
         diff = if_else(grepl('siRNA',tg.1), -1*diff, diff)) %>% 
  relocate(tg, cell.2, gene, diff, p.adj, comp ) %>% 
  
  # - join with the n.pep dataframe so we know how many peptides were measured
  left_join(., n.pep, by = c('gene' = 'Genes')) %>% 
  
  # remove proteins for which diff is NA, 
  # only retain comparisons that are within a cell line and include a comparison to control
  # or that are across cell lines and compare the two control samples
  filter(
    !is.na(diff), 
    (cell.1 == cell.2 & (grepl('siRNA', tg.1)|grepl('siRNA',tg.2))) 
    | (grepl('siRNA', tg.1) & grepl('siRNA',tg.2)) 
  ) 

write_tsv(de.results, here::here('manuscript','results','de_prot_results.tsv'))

# GSEA --------------------------------------------------------------------

# prepare the gene lists for GSEA
gsea = de.results %>% 
  
  # only keep the first gene symbol in the list
  mutate(g = str_remove_all(gene, '\\..*') %>% str_remove_all('|;.*')) %>% 
  
  # generate a gene list for each comparison
  group_by(tg, cell.2) %>% 
  summarise(
    all = diff %>% setNames(., g) %>% sort(decreasing = T) %>% list()
  ) %>% 
  pivot_longer(cols = c(all), values_to = 'gl') 

# run GSEA on the de results
gsea.results <- gsea %>%
  mutate(
    res = future_map(gl,
                     ~ clusterProfiler::gseGO(
                       geneList = .x, ont = 'ALL', keyType = 'SYMBOL',
                       OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                       eps = 0 , pvalueCutoff = 1
                     ))
  )

# Save results to a file
saveRDS(gsea.results, here::here('manuscript','results','gsea_results.rds'))

# clusterProfiler_4.12.6
# org.Mm.eg.db_3.19.1