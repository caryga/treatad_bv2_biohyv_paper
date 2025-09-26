# library(GEOquery)
# 
# Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 5)
# 
# id <- 'GSE181153' # Abeta Thp1 etc
# id <- 'GSE166500' # NHA astros
# gse <- getGEO(id)
# 
# length(gse)
# gse <- gse[[1]]

library(synapser)
library(cowplot)
library(tidyverse)

synLogin()

proj.dir <- here::here()

# pull score & omic data from synapse
tad.scores <- read_csv( synTableQuery("SELECT * FROM syn25575156")$filepath )
omics.scores <- read_csv( synTableQuery("SELECt * FROM syn22758536")$filepath )
gen.scores <- read_csv( synTableQuery('SELECT * FROM syn26844312')$filepath )

# pull biodom info from synapse
biodom <- readRDS( synGet('syn26642963')$path  )
biodom_genes <- readRDS(  synGet('syn26642963')$path  ) %>%
  select(Biodomain, GO_ID, symbol) %>% unnest_longer(symbol) %>% filter(!is.na(symbol), symbol != '')
dom.cols <- read_csv( synapser::synGet('syn26856828')$path )

theme_set(theme_bw())

pilot.tg.list <- tibble(
  hypothesis = c('mito_1a', 'mito_1b', 'immune_2a', 'immune_2b', 'syn_4a', 'syn_4b'),
  term = c('mitochondrion', 'mitochondrion', 'neutrophil degranulation', 'innate immune response', 'postsynaptic density', 'synapse'),
  GO_ID = c('GO:0005739','GO:0005739', 'GO:0043312','GO:0045087','GO:0014069','GO:0045202'),
  targets = list(
    c('ADH5','CAPN1','GABPA','KAT5','KIF1B','SQSTM1','STAT3','YY1'),
    c('ADH5', 'CAPN1', 'DLAT', 'DLD', 'PDHA1', 'PDHB', 'PDHX', 'PKM'),
    c('ANO6','AP2A2','CAPN1','IQGAP1','ITGAV','ITGB2','PKM','PSMC3','TMEM30A'),
    c('ARHGEF2', 'DDX1', 'DHX58', 'IFIH1', 'TBK1', 'TICAM1', 'TLR4', 'TRAF3'),
    c('ADD3','DLG4','PDLIM5','SLC6A9','STAT3'),
    c('CRKL','CYFIP2','DAG1','DLG4','GNA11','KIF1B','NDE1','NRXN1','TLN2')
  )) %>% unnest_longer(targets)

googlesheets4::gs4_auth('gregorycary@gmail.com')
gs = googlesheets4::gs4_get('1EHE9Kl7-ukxOIXEpQJB5mY-wGPRNw8-eMkEBTRKRvTc')

pilot.tg.list = googlesheets4::read_sheet(gs) %>% 
  select(target, BV2, `SH-SY5Y`, Both) %>% 
  distinct() %>% 
  left_join(pilot.tg.list, ., by = c('targets' = 'target'))

pilot.tg.list <- tad.scores %>% 
  select(targets = GeneName, ENSG, Overall, Overall_rank, GeneticsScore, OmicsScore) %>% 
  left_join(pilot.tg.list, ., by = 'targets')

pilot.tg.list <- omics.scores %>% 
  select(targets = GName, ENSG, RNA_Sig, RNA_Direction, RNA_TE, Pro_Sig, Pro_Direction, Pro_TE) %>% 
  left_join(pilot.tg.list, ., by = c('targets','ENSG'))

p = map(
  pilot.tg.list$hypothesis %>% unique(),
  ~ cowplot::plot_grid(
    omics.scores %>% 
      ggplot(aes(RNA_TE)) + 
      geom_histogram(
        data = subset(
          omics.scores,
          GName %in% (biodom_genes %>% filter(GO_ID %in% pilot.tg.list$GO_ID[which(pilot.tg.list$hypothesis == .x)]) %>% pull(symbol)) &
          RNA_Sig == 'NO'), alpha = .5, bins = 50) +
      geom_histogram(
        data = subset(
          omics.scores,
          GName %in% (biodom_genes %>% filter(GO_ID %in% pilot.tg.list$GO_ID[which(pilot.tg.list$hypothesis == .x)]) %>% pull(symbol)) &
          RNA_Sig == 'YES'), aes(fill = RNA_Direction), alpha = .5, bins = 50) +
      theme(legend.position = 'right')+labs(title = .x)
    ,
    omics.scores %>% 
      ggplot(aes(Pro_TE)) + 
      geom_histogram(
        data = subset(
          omics.scores,
          GName %in% (biodom_genes %>% filter(GO_ID %in% pilot.tg.list$GO_ID[which(pilot.tg.list$hypothesis == .x)]) %>% pull(symbol)) &
          Pro_Sig == 'NO'), alpha = .5, bins = 50) +
      geom_histogram(
        data = subset(
          omics.scores,
          GName %in% (biodom_genes %>% filter(GO_ID %in% pilot.tg.list$GO_ID[which(pilot.tg.list$hypothesis == .x)]) %>% pull(symbol)) &
          Pro_Sig == 'YES'), aes(fill = Pro_Direction), alpha = .5, bins = 50) +
      theme(legend.position = 'right') 
    , nrow = 1, align = 'hv')
)
cowplot::plot_grid(plotlist = p[c(1,3:6)], ncol = 1)


p = map(
  pilot.tg.list$hypothesis %>% unique(),
  ~ cowplot::plot_grid(
    omics.scores %>% 
      ggplot(aes(RNA_TE)) + 
      geom_histogram(
        data = subset(
          omics.scores,
          GName %in% (pilot.tg.list$targets[which(pilot.tg.list$hypothesis == .x)])  &
          RNA_Sig == 'NO'), alpha = .5, bins = 50) +
      geom_histogram(
        data = subset(
          omics.scores,
          GName %in% (pilot.tg.list$targets[which(pilot.tg.list$hypothesis == .x)])  &
          RNA_Sig == 'YES'), aes(fill = RNA_Direction), alpha = .5, bins = 50) +
      theme(legend.position = 'right')+labs(title = .x)
    ,
    omics.scores %>% 
      ggplot(aes(Pro_TE)) + 
      geom_histogram(
        data = subset(
          omics.scores,
          GName %in% (pilot.tg.list$targets[which(pilot.tg.list$hypothesis == .x)]) &
          Pro_Sig == 'NO'), alpha = .5, bins = 50) +
      geom_histogram(
        data = subset(
          omics.scores,
          GName %in% (pilot.tg.list$targets[which(pilot.tg.list$hypothesis == .x)]) &
          Pro_Sig == 'YES'), aes(fill = Pro_Direction), alpha = .5, bins = 50) +
      theme(legend.position = 'right') 
    , nrow = 1, align = 'hv')
)
cowplot::plot_grid(plotlist = p[c(1,3:6)], ncol = 1)


ranjita.list <- bind_rows(
  readxl::read_xlsx(paste0('~/treatAD_assayCellTypes', '/data/TREAT AD Targets and cell types Oct 2021.xlsx'), 
                    sheet = 1, skip = 2, col_names = c('Target','NA','Celltype')) %>% 
    mutate(round = 'one'),
  readxl::read_xlsx(paste0('~/treatAD_assayCellTypes', '/data/TREAT AD Targets and cell types Oct 2021.xlsx'), 
                    sheet = 2, skip = 2, col_names = c('Target','NA','Celltype')) %>% 
    mutate(round = 'two'),
  readxl::read_xlsx(paste0('~/treatAD_assayCellTypes', '/data/TREAT AD Targets and cell types Oct 2021.xlsx'), 
                    sheet = 1, skip = 2, col_names = c('Target','NA','Celltype')) %>% 
    mutate(round = 'smoc')
) %>% 
  select( -`NA`)

ranjita.list.og <- ranjita.list

tad.tg <- read_csv( synTableQuery('select * from syn25764512')$filepath )

ranjita.list <- 
  tibble( Target = setdiff(tad.tg$GENE, ranjita.list$Target),
          Celltype = NA,
          round = NA) %>% 
  bind_rows(ranjita.list, .)

m42.genes <- read_csv('~/treatAD_hypothesis/results/M42_genes.csv', col_names = F)

ranjita.list <- 
  tibble( Target = setdiff(m42.genes$X1, ranjita.list$Target),
          Celltype = NA,
          round = 'M42') %>% 
  bind_rows(ranjita.list, .)

rm(tad.tg, m42.genes)

all = tibble(geo = NULL, cell = NULL ,results = NULL)




# Proteomics --------------------------------------------------------------

# SH-SY5Y 
# PXD010776 # iTRAQ labeled SH-SY5Y, 4 replicates
system('wget http://ftp.ebi.ac.uk/pride-archive/2018/08/PXD010776/Rodriguez%20et%20al%202018%20iTRAQ%20note%20Tables%20S1-S6.xlsx')
readxl::excel_sheets('Rodriguez et al 2018 iTRAQ note Tables S1-S6.xlsx')
r = readxl::read_xlsx('Rodriguez et al 2018 iTRAQ note Tables S1-S6.xlsx', sheet = 1)
r$`Group Description` %>% str_extract(., 'GN=\\w+') %>% str_remove_all(., 'GN=') %>% unique() %>% length() # 3162 unique prot
intersect(r$`Group Description` %>% str_extract(., 'GN=\\w+') %>% str_remove_all(., 'GN=') %>% unique(),
          pilot.tg.list$targets %>% unlist()) %>% length() # 22 / 40 pilot proteins ID
# FOUND: ADH5 CAPN1 SQSTM1 STAT3 YY1 DLAT DLD PDHA1 PDHB PDHX PKM ANO6 AP2A2 IQGAP1 ITGAV PSMC3 ARHGEF2 DDX1 PDLIM5 CRKL GNA11 TLN2
# NOT FOUND: GABPA KAT5 KIF1B ITGB2 TMEM30A DHX58 IFIH1 TBK1 TICAM1 TLR4 TRAF3 ADD3 DLG4 SLC6A9 CYFIP2 DAG1 NDE1 NRXN1

# # PXD020807 # TMT of SH-SY5Y + ATRA (differentiation) 
# system('wget http://ftp.ebi.ac.uk/pride-archive/2020/10/PXD020807/2020_ATRA.xlsx')

# PXD031054 # SH-SY5Y , + RA (5d), +RA (3d) -> +RA+PMA (3d)
system('wget http://ftp.ebi.ac.uk/pride-archive/2022/04/PXD031054/artMS.zip') # 241 MB
x <- read_tsv( 'artMS/results/results-annotated.txt')
# FOUND: ADH5 CAPN1 GABPA SQSTM1 STAT3 YY1 DLAT DLD PDHA1 PDHB PDHX PKM ANO6 AP2A2 IQGAP1 ITGAV PSMC3 TMEM30A ARHGEF2 DDX1 TBK1 ADD3 PDLIM5 CRKL CYFIP2 GNA11 TLN2
# NOT FOUND: KAT5 KIF1B ITGB2 DHX58 IFIH1 TICAM1 TLR4 TRAF3 DLG4 SLC6A9 DAG1 NDE1 NRXN1

# PXD029012 # NGB overexpressing (NGB-FLAG) and control (CTRL)
system('wget https://ftp.ebi.ac.uk/pride-archive/2021/12/PXD029012/txt.zip') # 394 MB
# FOUND: ADH5 CAPN1 STAT3 DLAT DLD PDHA1 PDHB PDHX PKM ANO6 AP2A2 IQGAP1 ITGAV PSMC3 DDX1 ADD3 PDLIM5 CRKL CYFIP2 TLN2
# NOT FOUND: GABPA KAT5 KIF1B SQSTM1 YY1 ITGB2 TMEM30A ARHGEF2 DHX58 IFIH1 TBK1 TICAM1 TLR4 TRAF3 DLG4 SLC6A9 DAG1 GNA11 NDE1 NRXN1


# BV2
# PXD006753 # LFQ of BV2 micro + PBS | LPS (100ng/ml) | ShK-223 (100nM) or LPS+ShK-223; Seyfried/Dammer
system('wget http://ftp.ebi.ac.uk/pride-archive/2018/10/PXD006753/Rangaraju-BV2-MaxQuant_TXT_OUTPUT.zip') # 1143 MB

# PXD030447 # LPS-primed +/- baicalein
system('wget http://ftp.ebi.ac.uk/pride-archive/2022/02/PXD030447/PL%20pool%20ida2.group') # 1110 MB

# PXD006558 # BV2 cells + 1 ug/ml LPS | 10 ng/ml IFN-γ | 1ug/ml LPS plus 10 ng/ml IFR-γ; for 6 hr, 12 hr, 24 hr, and 48 hr
system('wget  http://ftp.ebi.ac.uk/pride-archive/2017/08/PXD006558/MAXQUANT_Result_Discovery.zip ')


# other
# PXD032733 # TMT of N2a with APPswe mutation


# Neurons -----------------------------------------------------------------


# GSE98834 - differentiated SH-SY5Y (RA+BDNF)
# groups: empty vector (pcDNA3), oe wt CLN1 (wtCLN1), oe mutatnt CLN1 (M57Nfs, L222P)
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98834/suppl/GSE98834_expression_table.txt.gz')
r = read_tsv(gzfile('GSE98834_expression_table.txt.gz'))
r <- r %>% 
  pivot_longer(cols = c(-Ensembl_geneid, -HGNC_gene_symbol, -HGNC_gene_name)) %>% 
  # filter(log2(value) > 2.5, !is.na(HGNC_gene_symbol)) %>%
  mutate(group = str_remove_all(name, '_[1-3]$'),
         group = case_when( grepl('pcDNA3', group) ~ 'empty_vector',
                            grepl('wtCLN1', group) ~ 'CLN1_wt',
                            grepl('L222P', group) ~ 'CLN1_L222P',
                            grepl('M57N', group) ~ 'CLN1_M57N')) %>% 
  rename(gene = HGNC_gene_symbol) %>% 
  group_by(gene, group) %>% 
  summarise(mn = mean(log2(value), na.rm = T),
            sd = sd(log2(value), na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, gene == tep.tg))+
#   geom_linerange(data = subset(r, gene == tep.tg),
#                  aes(ymin = mn-sd, ymax = mn+sd))+
#   coord_flip()

all = tibble(
  geo = 'GSE98834',
  cell = 'SH-SY5Y',
  data = list(r) ) %>% 
  bind_rows(all,.)

# GSE77383 - SH-SY5Y undiff & diff (RA-NBM)
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77383/suppl/GSE77383_expression_table.txt.gz')
r = read_tsv(gzfile('GSE77383_expression_table.txt.gz'))
r <- r %>% 
  pivot_longer(cols = c(-Ensembl_geneid, -HGNC_gene_symbol, -HGNC_gene_name)) %>% 
  # filter(log2(value) > 2.5, !is.na(HGNC_gene_symbol)) %>%
  mutate(group = str_remove_all(name, '_[1-3]$')) %>% 
  rename(gene = HGNC_gene_symbol) %>% 
  group_by(gene, group) %>% 
  summarise(mn = mean(log2(value), na.rm=T),
            sd = sd(log2(value), na.rm=T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, HGNC_gene_symbol == tep.tg))+
#   geom_linerange(data = subset(r, HGNC_gene_symbol == tep.tg), 
#                  aes(ymin = mn-sd, ymax = mn+sd))+
#   coord_flip()

all = tibble(
  geo = 'GSE77383',
  cell = 'SH-SY5Y',
  data = list(r) ) %>% 
  bind_rows(all,.)

# GSE74886 - SY-SY5Y undiff & diff (RA+BDNF)
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74886/suppl/GSE74886_processedFPKMS.xlsx')
r = readxl::read_xlsx('GSE74886_processedFPKMS.xlsx')
r <- r %>% 
  pivot_longer(cols = c(value_1, value_2) ) %>% 
  mutate(name = case_when( grepl('_1', name) ~ sample_1, grepl( '_2', name ) ~ sample_2)) %>% 
  select( gene = gene_id, group = name, mn = value ) %>% 
  mutate(sd = NA)
  # %>%   filter(log2(value) > 0, !is.na(gene_id)) 

# ggplot(r, aes(log2(value), group))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, gene == tep.tg))

all = tibble(
  geo = 'GSE74886',
  cell = 'SH-SY5Y',
  data = list(r) ) %>% 
  bind_rows(all,.)

# GSE150426 - HT22 (Mmus) +/- Nsmce1
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113436/suppl/GSE113436_gene_exp_with_symbol.txt.gz')

r = read_tsv(gzfile('GSE113436_gene_exp_with_symbol.txt.gz'))

ids = gprofiler2::gorth(r$Symbol, source_organism = 'mmusculus', target_organism = 'hsapiens')

r <- r %>% 
  left_join(., ids %>% select(input, gene = ortholog_name), by = c('Symbol'='input')) %>% 
  rename(Mmus_id = Symbol) %>% 
  pivot_longer(cols = 3:6) %>% 
  # filter( value > 2.5, !is.na(Symbol)) %>%
  mutate(group = str_remove_all(name, '[1-2]$'),
         group = case_when(group == 'N3' ~ 'Nsmce1_KO', group == 'PBIP' ~ 'WT')) %>% 
  group_by(Mmus_id, gene, group) %>% 
  summarise(mn = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, gene == tep.tg ))+
#   geom_linerange(data = subset(r, gene == tep.tg ),
#                  aes( ymin = mn -sd, ymax = mn +sd))+
#   coord_flip()

all = tibble(
  geo = 'GSE113436',
  cell = 'HT22',
  data = list(r) ) %>% 
  bind_rows(all,.)

# GSE127788 - HT22 +/- oe PGC1alpha +/- oxygen/glucose deprivation
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127788/suppl/GSE127788_matrix.txt.gz')

r = read_tsv(gzfile('GSE127788_matrix.txt.gz')) %>% 
  mutate(enst = str_remove_all(Track_id, '\\..*'))

ids = gprofiler2::gorth(r$enst, source_organism = 'mmusculus', target_organism = 'hsapiens')

r <- r %>% 
  left_join(., ids %>% select(input, gene = ortholog_name), by = c('enst'='input')) %>% 
  rename(Mmus_id = enst) %>% 
  pivot_longer(cols = 2:5, names_to = 'group') %>% 
  # filter( log2(value) > 2.5, !is.na(target)) %>% 
  group_by(Mmus_id, gene, group) %>% 
  summarise(mn = mean(log2(value), na.rm = T),
            sd = sd(log2(value), na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn ))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, hsap == tep.tg ))+
#   geom_linerange(data = subset(r, hsap == tep.tg ),
#                  aes(ymin = mn - sd, ymax = mn + sd))+
#   coord_flip()

all = tibble(
  geo = 'GSE127788',
  cell = 'HT22',
  data = list(r) ) %>% 
  bind_rows(all,.)


# GSE109699 - HT22 +/- Pja2
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109699/suppl/GSE109699_gene_exp_with_symbol.txt.gz')

r = read_tsv(gzfile('GSE109699_gene_exp_with_symbol.txt.gz')) 

ids = gprofiler2::gorth(r$Gene_ID, source_organism = 'mmusculus', target_organism = 'hsapiens')

r <- r %>% 
  left_join(., ids %>% select(input, gene = ortholog_name), by = c('Gene_ID'='input')) %>% 
  rename(Mmus_id = Gene_ID) %>% 
  pivot_longer(cols = 3:6) %>% 
  # filter( value > 2.5, !is.na(Symbol)) %>% 
  mutate(group = str_remove_all(name, '[1-2]$'),
         group = case_when(group == 'PBIP' ~ 'wt', group == 'P2' ~ 'Pja2_OE')) %>% 
  group_by(Mmus_id, gene, group) %>% 
  summarise(mn = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn ))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, gene == tep.tg))+
#   geom_linerange(data = subset(r, gene == tep.tg),
#                  aes(ymin = mn - 3*sd, ymax = mn + 3*sd))+
#   coord_flip()

all = tibble(
  geo = 'GSE109699',
  cell = 'HT22',
  data = list(r) ) %>% 
  bind_rows(all,.)


# Immune ------------------------------------------------------------------

# GSE181153 - Abeta stim of immune response in immune cell culture systems
# cell types: iPSC, iHPC, iMGL, HMC3, U87a, PBMC, THP1
# treatments: PBS (control), LPS and IFN-gamma, oligmerized Abeta, fibrilized Abeta
# ADAB_RNA_<CellType>_<Genotype>_<AbetaType>_<AbetaConc>_<OtherTreatments>_<Serum>_<Time>_S_quant
 
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181153/suppl/GSE181153_ADAB_Fig3_DESeqDesigns.xlsx')
# d = readxl::read_xlsx('GSE181153_ADAB_Fig3_DESeqDesigns.xlsx')

# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181153/suppl/GSE181153_ADAB_geneCounts.tsv.gz')

r = read_tsv(gzfile('GSE181153_ADAB_geneCounts.tsv.gz'))

r <- r %>% 
  pivot_longer(cols = c(-ensembl, -hgnc)) %>% 
  mutate(group = str_remove_all(name, 'ADAB_RNA_|\\.[1-2]\\.1$')) %>% 
  rename(gene = hgnc) %>% 
  # separate(name, c('cell','geno','abeta','conc','tx','serum','time','s','expt','rep')) %>% 
  # filter(log2(value) > 2.5, !is.na(hgnc)) %>%
  group_by(gene, group) %>% 
  summarise(mn = mean(log2(value), na.rm = T),
            sd = sd(log2(value), na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, gene == tep.tg)
#              # , position = position_jitter(seed = 123, width =0.2)
#              ) +
#   geom_linerange(data = subset(r, gene == tep.tg),
#                  aes(ymin = mn-sd, ymax = mn+sd)
#                  # , position = position_jitter(seed = 123, width =0.2)
#                  )+
#   coord_flip()


all = tibble(
  geo = 'GSE181153',
  cell = 'Immune_multi',
  data = list(r) ) %>% 
  bind_rows(all,.)


# GSE162526 - BV2 PU.1 knock-down vs overexpression
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162526/suppl/GSE162526_Spi1_KD.tsv.gz')
r = read_tsv(gzfile('GSE162526_Spi1_KD.tsv.gz'))
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162526/suppl/GSE162526_Spi1_OE.tsv.gz')
r = read_tsv(gzfile('GSE162526_Spi1_OE.tsv.gz')) %>% 
  full_join(r, ., by = 'gene')

ids = gprofiler2::gorth(r$gene, source_organism = 'mmusculus', target_organism = 'hsapiens')

r <- r %>% 
  rename( Mmus_id = gene ) %>% 
  left_join(., ids %>% select(input, gene = ortholog_name), by = c('Mmus_id'='input')) %>% 
  pivot_longer(cols = 2:13) %>% 
  # filter( log2(value) > 2.5, !is.na(gene)) %>% 
  mutate(group = str_remove_all(name, '[1-3]$'),
         group = case_when( grepl('KD_S', group) ~ 'Spi_KD',
                            grepl('KD_C', group) ~ 'Spi_KD_control',
                            grepl('OE_S', group) ~ 'Spi_OE',
                            grepl('OE_C', group) ~ 'Spi_OE_control') ) %>% 
  group_by(Mmus_id, gene, group) %>% 
  summarise(mn = mean(log2(value), na.rm = T),
            sd = sd(log2(value), na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn ))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, gene == tep.tg ))+
#   geom_linerange(data = subset(r, gene == tep.tg ),
#                  aes(ymin = mn - sd, ymax = mn + sd))+
#   coord_flip()


all = tibble(
  geo = 'GSE162526',
  cell = 'BV2',
  data = list(r) ) %>% 
  bind_rows(all,.)


# GSE132739 - BV2 + IFNg + Atg5KO
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132739/suppl/GSE132739_Orvedahl_et_al_BV2_IFNg_RNAseq_log2norm.xlsx')

r = readxl::read_xlsx('GSE132739_Orvedahl_et_al_BV2_IFNg_RNAseq_log2norm.xlsx') 

ids = gprofiler2::gorth(r$GeneID, source_organism = 'mmusculus', target_organism = 'hsapiens')

r <- r %>% 
  left_join(., ids %>% select(input, gene = ortholog_name), by = c('GeneID'='input')) %>% 
  rename(Mmus_id = Gene_symbol) %>% 
  pivot_longer(cols = -c(1:3,16)) %>% 
  # filter( value > 2.5, !is.na(Gene_symbol)) %>% 
  mutate(group = str_remove_all(name, '_Rep[A-C]$')) %>% 
  group_by(Mmus_id, gene, group) %>% 
  summarise(mn = mean(log2(value), na.rm =T),
            sd = sd(log2(value), na.rm =T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn ))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, hsap == tep.tg ))+
#   geom_linerange(data = subset(r, hsap == tep.tg ),
#                  aes(ymin = mn - sd, ymax = mn + sd))+
#   coord_flip()


all = tibble(
  geo = 'GSE132739',
  cell = 'BV2',
  data = list(r) ) %>% 
  bind_rows(all,.)


# GSE137741 - BV2 +LPS doses
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137741/suppl/GSE137741_count_RPMs_RPKMs.xls.gz')
# system('gunzip GSE137741_count_RPMs_RPKMs.xls.gz')
r = readxl::read_xls('GSE137741_count_RPMs_RPKMs.xls') 

ids = gprofiler2::gorth(r$...1, source_organism = 'mmusculus', target_organism = 'hsapiens')

r <- r %>% 
  left_join(., ids %>% select(input, gene = ortholog_name), by = c('...1'='input')) %>% 
  rename(Mmus_id = GeneName) %>% 
  pivot_longer(cols = -c(1:3,13)) %>% 
  # filter( log2(value) > 2.5, !is.na(GeneName)) %>% 
  mutate(group = str_remove_all(name, ' [4-6]$'),
         group = case_when(group == 'unprimed' ~ '1xLPS', 
                           group == 'high primed' ~ '2xLPS_highDose',
                           group == 'ultra-low primed' ~ '2xLPS_lowDose')) %>% 
  group_by(Mmus_id, gene, group) %>% 
  summarise(mn = mean(log2(value), na.rm = T),
            sd = sd(log2(value), na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn ))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, gene == tep.tg ))+
#   geom_linerange(data = subset(r, gene == tep.tg ),
#                  aes(ymin = mn - sd, ymax = mn + sd))+
#   coord_flip()

all = tibble(
  geo = 'GSE137741',
  cell = 'BV2',
  data = list(r) ) %>% 
  bind_rows(all,.)


# GSE173381 - BV2 +/- FL-APOE4


# Astrocytes --------------------------------------------------------------


# GSE166500 - normal human astrocytes (NHA)
# groups: Palmitic Acid, Tibolone, and Palmitic Acid + Tibolone; 
# controls: 
  # 1. DMEM, for the treatment with Tibolone, and 
  # 2. Vehicle, for the treatments with Palmitic Acid and Palmitic Acid + Tibolone.
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166500/suppl/GSE166500_salmon_merged_gene_counts.csv.gz')
# c = read_csv(gzfile('GSE166500_salmon_merged_gene_counts.csv.gz'))

# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166500/suppl/GSE166500_count_matrix_DESeq2_normalised.csv.gz')

r = read_csv(gzfile('GSE166500_count_matrix_DESeq2_normalised.csv.gz'))

r <- r %>% 
  pivot_longer(cols = c(-ensembl_gene_id, -external_gene_name)) %>%
  # filter(log2(value) > 2.5, !is.na(external_gene_name)) %>%
  mutate(group = str_remove_all(name, '[0-9]_[0-9]+$')) %>% 
  group_by(external_gene_name, group) %>% 
  summarise(mn = mean(log2(value), na.rm = T),
            sd = sd(log2(value), na.rm = T)) %>% 
  ungroup()

ggplot(r, aes(group, mn))+ theme_bw()+
  geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
  geom_point(data = subset(r, external_gene_name == tep.tg)) +
  geom_linerange(data = subset(r, external_gene_name == tep.tg),
                 aes(ymin = mn - sd, ymax = mn + sd)) +
  coord_flip()


all = tibble(
  geo = 'GSE166500',
  cell = 'NHA',
  data = list(r) ) %>% 
  bind_rows(all,.)


# GSE140471 - NHA WT & ATXN1 KO
# CIC knockout monoclonal lines are labeled A2 and H9. ATXN1L knockout monoclonal lines are labeled B82 and B16
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140471/suppl/GSE140471_Dwong_RNAseq_Raw-Counts.txt.gz')
r = read_tsv(gzfile('GSE140471_Dwong_RNAseq_Raw-Counts.txt.gz')) %>% 
  rename(gene = ...1)

r <- r %>% 
  pivot_longer(cols = c(-gene)) %>% 
  # filter(log2(value) > 2.5, !is.na(gene)) %>%
  mutate(group = str_remove_all(name, '\\.[1-3]$'),
         expt = case_when( group == 'A2' ~ 'CIC_ko_1',
                           group == 'H9' ~ 'CIC_ko_2',
                           group == 'B82' ~ 'ATXN1L_ko_1',
                           group == 'B16' ~ 'ATXN1L_ko_2',
                           group == 'NHA' ~ 'NHA'))  %>% 
  group_by(gene, expt, group) %>% 
  summarise(mn = mean(log2(value), na.rm = T),
            sd = sd(log2(value), na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes(expt, mn))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, gene == tep.tg)) +
#   geom_linerange(data = subset(r, gene == tep.tg),
#                  aes(ymin = mn - sd, ymax = mn + sd))+
#   coord_flip()

all = tibble(
  geo = 'GSE140471',
  cell = 'NHA',
  data = list(r) ) %>% 
  bind_rows(all,.)

# GSE179882 - HDAC1 knock-down in astrocyte lines 
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179882/suppl/GSE179882_gexp_counts_sampleTitle.csv.gz')

r = read_csv(gzfile('GSE179882_gexp_counts_sampleTitle.csv.gz'), skip = 1) %>% 
  mutate(ensg = str_remove_all(Title, '\\..*'))

ids = gprofiler2::gconvert(r$ensg, organism = 'hsapiens', target = 'HGNC' )

r <- r %>% 
  left_join(., ids %>% select(input, symbol = target), by = c('ensg'= 'input')) %>% 
  pivot_longer(cols = c(-Title, -ensg, -symbol)) %>% 
  # filter(log2(value) > 2.5, !is.na(symbol)) %>%
  mutate(group = str_remove_all(name, ' Rep[1-3]$')) %>% 
  group_by(symbol, group) %>% 
  summarise(mn = mean(log2(value), na.rm = T),
            sd = sd(log2(value), na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes(group, mn))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, symbol == tep.tg)) +
#   geom_linerange(data = subset(r, symbol == tep.tg), 
#                  aes(ymin = mn - sd, ymax = mn + sd))+
#   coord_flip()

all = tibble(
  geo = 'GSE179882',
  cell = 'NHA',
  data = list(r) ) %>% 
  bind_rows(all,.)

# Endothelium -------------------------------------------------------------


# GSE58663 - HUVEC + VEGFA or Histamine
# system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58663/suppl/GSE58663_Gene_Expression_RPKM.txt.gz')

r = read_tsv(gzfile('GSE58663_Gene_Expression_RPKM.txt.gz')) 

r <- r %>% 
  pivot_longer(cols = starts_with('HUVEC_RNA-Seq_')) %>% 
  mutate(name = str_remove_all(name, 'HUVEC_RNA-Seq_|_RPKM')) %>% 
  mutate(group = str_remove_all(name, '_rep[1-2]$')) %>% 
  # filter(log2(value) > 2.5, !is.na(Symbol))  %>% 
  group_by(Symbol, group) %>% 
  summarise(mn = mean(log2(value), na.rm = T),
            sd = sd(log2(value), na.rm = T)) %>% 
  ungroup()

# ggplot(r, aes( group, mn))+ theme_bw()+
#   geom_violin(fill = 'grey90', draw_quantiles = c(.25,.5, .75), trim = T)+
#   geom_point(data = subset(r, Symbol == tep.tg))+
#   geom_linerange(data = subset(r, Symbol == tep.tg), 
#                 aes(ymin = mn - sd, ymax = mn + sd))+
#   coord_flip()

all = tibble(
  geo = 'GSE58663',
  cell = 'HUVEC',
  data = list(r) ) %>% 
  bind_rows(all,.)


# filter ------------------------------------------------------------------

x = all %>% slice(1:2) %>% unnest(data)


for(i in 21:33){
  tep.tg = ranjita.list$Target[i]
  p = ggplot(x, aes( group, mn, fill = geo))+ theme_bw()+
    geom_violin( draw_quantiles = c(.25,.5, .75), trim = T, alpha = .5)+
    geom_point(data = subset(x, HGNC_gene_symbol == tep.tg))+
    geom_linerange(data = subset(x, HGNC_gene_symbol == tep.tg),
                  aes(ymin = mn - sd, ymax = mn + sd))+
    coord_flip() + labs(title = tep.tg)
  print(p)
}
