library(tidyverse)

# Karen-Shul --------------------------------------------------------------

# Ido Amit DAM PMID: 28602351
url='https://ars.els-cdn.com/content/image/1-s2.0-S0092867417305780-mmc3.xlsx'
download.file(url, destfile =  here::here('data','dam_amit.xlsx'))

dam_amit = readxl::read_xlsx(here::here('data','dam_amit.xlsx'),
                             col_types = c('text','numeric','numeric','numeric','numeric')) %>% 
  filter(`DAM FDR p-value` > -log10(5e-2),
         `Fold-change (DAM to homeostatic microglia)` > 0) %>% 
  arrange(desc(`DAM FDR p-value`)) %>% 
  pull(`Gene name`)

homeo_amit = readxl::read_xlsx(here::here('data','dam_amit.xlsx'), 
                               col_types = c('text','numeric','numeric','numeric','numeric')) %>% 
  filter(`DAM FDR p-value` > -log10(5e-2),
         `Fold-change (DAM to homeostatic microglia)` < 0) %>% 
  arrange(desc(`DAM FDR p-value`)) %>% 
  pull(`Gene name`)
# 'Karen-Shaul_1-s2.0-S0092867417305780-mmc2.xlsx' 

cat('Amit DAM: ',length(dam_amit),
    '\nAmit Homeo: ', length(homeo_amit))


# Zhou --------------------------------------------------------------------

# Colona PMID: 31932797
url = 'https://pmc.ncbi.nlm.nih.gov/articles/instance/6980793/bin/NIHMS1542659-supplement-Sup_tab4.xlsx'
download.file(url, destfile =  here::here('data','colona_subsets.xlsx'))

# Micro1 is more abundant in Control
homeo.hs_colona = readxl::read_xlsx(here::here('data','colona_subsets.xlsx'), 13) %>%
  arrange(p_val_adj) %>%
  pull(gene)

# Micro0 is more abundant in AD
dam.hs_colona = readxl::read_xlsx(here::here('data','colona_subsets.xlsx'), 12) %>% 
  arrange(p_val_adj) %>% 
  pull(gene)

cat('Colonna DAM: ',length(dam.hs_colona),
    '\nColonna Homeo: ', length(homeo.hs_colona))


# Rangaraju ---------------------------------------------------------------

# Levey PMID: 29784049
url='https://pmc.ncbi.nlm.nih.gov/articles/instance/5963076/bin/13024_2018_254_MOESM1_ESM.xlsx'
download.file(url, destfile =  here::here('data','levey_subsets.xlsx'))

proinf.dam_levey <- readxl::read_xlsx(here::here('data','levey_subsets.xlsx'), 1, skip = 1) %>% 
  filter(Module == 'magenta'
         , kMEmagenta > quantile(kMEmagenta, .975)
  ) %>% 
  pull(GeneSymbol)

antiinf.dam_levey <- readxl::read_xlsx(here::here('data','levey_subsets.xlsx'), 1, skip = 1) %>% 
  filter(Module == 'yellow'
         , kMEyellow > quantile(kMEyellow, .975)
  ) %>% 
  pull(GeneSymbol)

homeo.3_levey <- readxl::read_xlsx(here::here('data','levey_subsets.xlsx'), 1, skip = 1) %>% 
  filter(Module == 'blue'
         , kMEblue > quantile(kMEblue, .975)
  ) %>% 
  pull(GeneSymbol)


cat('Levey Proinflam DAM: ',length(proinf.dam_levey),
    '\nLevey Antiinflam DAM: ',length(antiinf.dam_levey),
    '\nLevey Homeo: ', length(homeo.3_levey))


# Lloyd -------------------------------------------------------------------

url = 'https://ars.els-cdn.com/content/image/1-s2.0-S2211124724012592-mmc2.xlsx'
download.file(url, destfile =  here::here('data','lloyd_subsets.xlsx'))


uglia.modules <- readxl::read_xlsx(here::here('data','lloyd_subsets.xlsx'), sheet = 4) %>% 
  mutate(nm = str_c('[',`Protein Module Colour`,']: ',`Module GO Term`)) %>%
  group_by(nm) %>% summarise(gene = list(`Gene Name`))


# join & save -------------------------------------------------------------

cell.subsets = tibble(
  set = c('dam_amit',
          'homeo_amit',
          'homeo.hs_colona',
          'dam.hs_colona',
          'proinf.dam_levey',
          'antiinf.dam_levey',
          'homeo.3_levey'
  ),
  gene = list(dam_amit,
              homeo_amit,
              homeo.hs_colona,
              dam.hs_colona,
              proinf.dam_levey,
              antiinf.dam_levey,
              homeo.3_levey
  )) %>% 
  bind_rows(.,uglia.modules %>% select(set = nm, gene))

# call orthologs for Hs genes
cell.subsets$mm.gene = cell.subsets$gene
for(i in which( map_lgl(cell.subsets$gene, ~ all( head(.x) == str_to_upper(head(.x)) )) )  ){
  x = gprofiler2::gorth(cell.subsets$gene[[i]], 'hsapiens','mmusculus')
  cell.subsets$mm.gene[[i]] <- x$ortholog_name
}

saveRDS(cell.subsets, here::here('data','cellSubsets_DAM_etc.rds'))
  