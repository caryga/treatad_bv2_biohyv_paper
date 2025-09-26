
# packages ----------------------------------------------------------------

syn.client <- reticulate::import('synapseclient')
syn.utils <- reticulate::import('synapseutils')
syn <- syn.client$Synapse()

library(tidyverse)

diagnostic_plots <- function(lmm.fit){
  
  # ── 1. Violin + box plots of residuals by batch ──────────────────────────────
  diag_df <- broom.mixed::augment(lmm.fit)                       # adds .resid and .fitted
  
  p_batch <- ggplot(diag_df,
                    aes(x = batch, y = .resid, fill = batch)) +
    geom_violin(trim = FALSE, width = 1, colour = NA, alpha = 0.4) +
    geom_boxplot(width = 0.12, outlier.shape = NA, colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Residual distribution by batch",
         x = "Batch", y = "Model residuals") +
    theme_bw() +
    theme(legend.position = "none")
  
  # ── 2. Residual-vs-fitted plot ───────────────────────────────────────────────
  p_rvf <- ggplot(diag_df, aes(x = .fitted, y = .resid)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Residuals vs fitted values",
         x = "Fitted values", y = "Residuals") +
    theme_bw()
  
  # ── 3. Forest plot of batch-level random effects + overall estimate ──────────
  # random-effect BLUPs with 95% CIs
  re_df <- broom.mixed::tidy(lmm.fit, effects = "ran_vals", conf.int = TRUE) %>% 
    filter(term == "(Intercept)", level != "(Intercept)")  # keep batch intercepts
  
  # overall (fixed) intercept
  fix_int <- broom.mixed::tidy(lmm.fit, effects = "fixed", conf.int = TRUE) %>% 
    filter(term == "(Intercept)") %>% 
    mutate(level = "Overall",  # to match column names
           grp = NA)           # placeholder
  
  # combine for plotting; order batches by BLUP
  forest_df <- re_df %>% 
    arrange(estimate) %>% 
    mutate(level = factor(level, levels = level))     
  
  p_forest <- ggplot(forest_df,
                     aes(x = estimate, y = level)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    # add diamond for overall estimate
    geom_point(data = fix_int, shape = 18, size = 4, fill = "black") +
    geom_vline(xintercept = fix_int$estimate, linetype = "dotted") +
    labs(title = "Batch intercepts (BLUPs) with overall estimate",
         x = "Log2 fold-change (BLUP)", y = NULL) +
    theme_bw()
  
  p <- cowplot::plot_grid( p_batch / p_rvf, p_forest, nrow = 1, rel_widths = c(1,.5))
  
  p
}

syn$login()

biohyv_targets = read_tsv( syn$get('syn65987195') %>% .$path, show_col_types = F) %>% 
  mutate(
    hyp_area = case_when(grepl('immune',hypothesis) ~ 'Immune Response targets',
                         grepl('mito', hypothesis) ~ 'Mitochondrial Metabolism targets')) %>% 
  select(target = targets, hyp_area) %>% 
  group_by(target) %>% 
  summarise(hyp_area = paste0(unique(hyp_area), collapse = ' | ')) %>% 
  distinct()

# alamarBlue --------------------------------------------------------------

syn$login()
synid = 'syn65941598'

alamar_blue <- readxl::read_xlsx( syn$get(synid) %>% .$path )

blank <- alamar_blue %>% 
  filter(value == 'blank') %>% 
  select(batch,blank = rep1)

alamar_blue <- alamar_blue %>% 
  filter(value != 'blank') %>% 
  left_join(., blank) %>% 
  mutate( across(starts_with('rep'), ~ .x - blank ))

vars <- alamar_blue %>% select(starts_with('rep')) %>% names()

norm_ab <- alamar_blue %>%
  rowwise() %>% 
  mutate( 
    avg = mean(c_across(all_of(vars))),
    sd = sd(c_across(all_of(vars))),
    upper = avg + sd,
    lower = avg - sd
  ) %>% 
  ungroup() %>% 
  mutate(
    ctrl.avg = mean(avg[siRNA == 'Control']),
    ctrl.upr = mean(upper[siRNA == 'Control']),
    ctrl.lwr = mean(lower[siRNA == 'Control']),
    .by = c(BV2.cell.line, batch)
  ) %>% 
  group_by(BV2.cell.line, batch, siRNA) %>% 
  mutate(

    tg.avg = mean(avg),
    fc = tg.avg / ctrl.avg,
    lfc = log2(fc),

    tg.lower = mean(lower),
    lower_fc = tg.lower / ctrl.lwr,
    log_lower = log2(lower_fc),

    tg.upper = mean(upper),
    upper_fc = tg.upper / ctrl.upr,
    log_upper = log2(upper_fc)
         ) %>% 
  ungroup() %>% 
  mutate(
    batch         = factor(batch),
    BV2.cell.line = factor(BV2.cell.line),
    siRNA         = factor(siRNA)        
  ) %>% 
  mutate(
    siRNA = relevel(siRNA, ref = "Control"),
    BV2.cell.line = relevel(BV2.cell.line, ref = "scramble")
  ) 

# Fit Linear Mixed Model

fit <- lme4::lmer(lfc ~ siRNA * BV2.cell.line + (1 | batch), data = norm_ab)

emm  <- emmeans::emmeans(fit, ~ siRNA | BV2.cell.line)
ctr  <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = "Control", adjust = "none")

tidy_ctr <- broom::tidy(ctr, conf.int = TRUE) %>% 
  mutate( p.adj.BH   = p.adjust(p.value, method = "BH") )

# # model diagnostics
# p <- diagnostic_plots(fit)

ablue_lmm = tidy_ctr %>% 
  mutate(siRNA = str_split_fixed(contrast, ' - ', 2)[,1]) %>% 
  left_join(norm_ab, ., by = c('siRNA','BV2.cell.line')) %>% 
  mutate(pheno = 'alamar_blue') %>% relocate(pheno)

# mitotracker -------------------------------------------------------------

syn$login()

synid <- 'syn65941595'

mito_tracker <- readxl::read_xlsx( syn$get(synid) %>% .$path , sheet = 'fluor_intensity' )
cell_number <-  readxl::read_xlsx( syn$get(synid) %>% .$path , sheet = 'total_cells' )

norm_mt <- mito_tracker 

n <- which(grepl('rep',colnames(mito_tracker)))

# norm_mt[,n] <- norm_mt[,n] / cell_number[,n]

vars <- norm_mt %>% select(starts_with('rep')) %>% names()

norm_mt <- norm_mt %>%
  # filter(!grepl('2023|2025.04',batch)) %>% 
  rowwise() %>% 
  mutate( 
    avg = mean(c_across(all_of(vars)), na.rm =T),
    sd = sd(c_across(all_of(vars)), na.rm =T),
    upper = avg + sd,
    lower = avg - sd
  ) %>% 
  ungroup() %>% 
  mutate(
    ctrl.avg = mean(avg[siRNA == 'Control']),
    ctrl.upr = mean(upper[siRNA == 'Control']),
    ctrl.lwr = mean(lower[siRNA == 'Control']),
    .by = c(BV2.cell.line, batch)
  ) %>% 
  group_by(BV2.cell.line, siRNA, batch) %>% 
  mutate(
    tg.avg = mean(avg),
    fc = tg.avg / ctrl.avg,
    lfc = log2(fc),

    tg.lower = mean(lower),
    lower_fc = tg.lower / ctrl.lwr,
    log_lower = log2(lower_fc),

    tg.upper = mean(upper),
    upper_fc = tg.upper / ctrl.upr,
    log_upper = log2(upper_fc)
  ) %>% 
  ungroup() %>% 
  mutate(
    batch         = factor(batch),
    BV2.cell.line = factor(BV2.cell.line),
    siRNA         = factor(siRNA)
  ) %>% 
  mutate(
    siRNA = relevel(siRNA, ref = "Control"),
    BV2.cell.line = relevel(BV2.cell.line, ref = "scramble")
    ) 

# Fit Linear Mixed Model

fit <- lme4::lmer(lfc ~ siRNA * BV2.cell.line + (1 | batch), data = norm_mt
                  # %>% filter((batch %in% c('2025.04.28','2023.01.11'))==F)
                  )

emm  <- emmeans::emmeans(fit, ~ siRNA | BV2.cell.line)
ctr  <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = "Control", adjust = "none") 


tidy_ctr <- broom::tidy(ctr, conf.int = TRUE) %>% 
  mutate( p.adj.BH   = p.adjust(p.value, method = "BH") )

# # model diagnostics
# p <- diagnostic_plots(fit)

mitotracker_lmm = tidy_ctr %>% 
  mutate(siRNA = str_split_fixed(contrast, ' - ', 2)[,1]) %>% 
  left_join(norm_mt, ., by = c('siRNA','BV2.cell.line')) %>% 
  mutate(pheno = 'mitotracker') %>% relocate(pheno)

# phagocytosis ------------------------------------------------------------

syn$login()

synid <- 'syn65941596'

phagocytosis <- readxl::read_xlsx( syn$get(synid) %>% .$path , sheet = 'pct_positive_cells' )
cell_number <-  readxl::read_xlsx( syn$get(synid) %>% .$path , sheet = 'total_cells' )

vars <- phagocytosis %>% select(starts_with('rep')) %>% names()

# phagocytosis <- phagocytosis %>% filter(passage_n != 8)

norm_ph <- phagocytosis %>% 
  rename(batch = passage_n) %>% 
  rowwise() %>% 
  mutate( 
    avg = mean(c_across(all_of(vars)), na.rm = T),
    sd = sd(c_across(all_of(vars)), na.rm = T),
    upper = avg + sd,
    lower = avg - sd
  ) %>% 
  ungroup() %>% 
  mutate(
    ctrl.avg = mean(avg[siRNA == 'Control'], na.rm=T),
    ctrl.upr = mean(upper[siRNA == 'Control'], na.rm=T),
    ctrl.lwr = mean(lower[siRNA == 'Control'], na.rm=T),
    .by = c(BV2.cell.line,batch)
  ) %>% 
  group_by(BV2.cell.line,batch,siRNA) %>% 
  mutate(
 
    tg.avg = mean(avg),
    fc = tg.avg / ctrl.avg,
    lfc = log2(fc),
 
    tg.lower = mean(lower),
    lower_fc = tg.lower / ctrl.lwr,
    log_lower = log2(lower_fc),
  
    tg.upper = mean(upper),
    upper_fc = tg.upper / ctrl.upr,
    log_upper = log2(upper_fc)
  ) %>% 
  ungroup() %>% 
  mutate(
    batch         = factor(batch),
    BV2.cell.line = factor(BV2.cell.line),
    siRNA         = factor(siRNA)
  ) %>% 
  mutate(
    siRNA = relevel(siRNA, ref = "Control"),
    BV2.cell.line = relevel(BV2.cell.line, ref = "scramble")
  ) 

# Fit Linear Mixed Model

fit <- lme4::lmer(lfc ~ siRNA * BV2.cell.line + (1 | batch), data = norm_ph)

emm  <- emmeans::emmeans(fit, ~ siRNA | BV2.cell.line)
ctr  <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = "Control", adjust = "none") 

tidy_ctr <- broom::tidy(ctr, conf.int = TRUE) %>% 
  mutate( p.adj.BH   = p.adjust(p.value, method = "BH") )

# # model diagnostics
# p <- diagnostic_plots(fit)

phagocytosis_lmm = tidy_ctr %>% 
  mutate(siRNA = str_split_fixed(contrast, ' - ', 2)[,1]) %>% 
  left_join(norm_ph, ., by = c('siRNA','BV2.cell.line')) %>% 
  mutate(pheno = 'phagocytosis') %>% relocate(pheno)

# nfkb --------------------------------------------------------------------

syn$login()

synid <- 'syn65941597'

nfkb <- readxl::read_xlsx( syn$get(synid) %>% .$path , sheet = 1 ) %>% 
  mutate(tx = if_else(LPS.dose == '100 ng/ml', 'LPS', 'none'))

vars <- nfkb %>% select(starts_with('rep')) %>% names()

# nfkb$rep1[which(
#   nfkb$BV2.cell.line == 'scramble' &
#     nfkb$siRNA == 'DLD' &
#     nfkb$LPS.dose == "0 ng/ml")] <- NA_real_

norm_nfkb <- nfkb %>%
  filter(!grepl('2023',batch)) %>%
  rowwise() %>% 
  mutate( 
    # batch = '2024',
    avg = mean(c_across(all_of(vars)), na.rm = T),
    sd = sd(c_across(all_of(vars)), na.rm = T),
    upper = avg + sd,
    lower = avg - sd,
    # tx = str_c(BV2.cell.line,'.',tx)
  ) %>% 
  ungroup() %>% 
  mutate(
    ctrl.avg = mean(avg[siRNA == 'Control']),
    ctrl.upr = mean(upper[siRNA == 'Control']),
    ctrl.lwr = mean(lower[siRNA == 'Control']),
    .by = c(BV2.cell.line, tx, batch)
  ) %>% 
  group_by(BV2.cell.line, tx, siRNA, batch) %>% 
  mutate(

    tg.avg = mean(avg),
    fc = tg.avg / ctrl.avg,
    lfc = log2(fc),

    tg.lower = mean(lower),
    lower_fc = tg.lower / ctrl.lwr,
    log_lower = log2(lower_fc),

    tg.upper = mean(upper),
    upper_fc = tg.upper / ctrl.upr,
    log_upper = log2(upper_fc)
  ) %>% 
  ungroup() %>% 
  mutate(
    batch         = factor(batch),
    BV2.cell.line = factor(BV2.cell.line),
    tx            = factor(tx),
    siRNA         = factor(siRNA)  
  ) %>% 
  mutate(
    siRNA = relevel(siRNA, ref = "Control"),
    tx = relevel(tx, ref = "none"),
    BV2.cell.line = relevel(BV2.cell.line, ref = "scramble")
  ) 

# Fit Linear Mixed Model

# No LPS

fit <- lme4::lmer(lfc ~ siRNA * BV2.cell.line + (1 | batch), data = norm_nfkb %>% filter(tx == 'none'))
# fit <- lm(lfc ~ siRNA * BV2.cell.line, data = norm_nfkb %>% filter(tx == 'none', !is.na(lfc)))

emm  <- emmeans::emmeans(fit, ~ siRNA | BV2.cell.line)
ctr  <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = "Control", adjust = "none")

tidy_ctr <- broom::tidy(ctr, conf.int = TRUE) |>
  mutate(
    p.adj.BH   = p.adjust(p.value, method = "BH"),          # FDR (Benjamini-Hochberg)
    p.adj.bonf = p.adjust(p.value, method = "bonferroni")   # Bonferroni (optional)
  )

tidy_ctr_noLPS <- tidy_ctr %>% mutate(tx = 'none')

# LPS treated

fit <- lme4::lmer(lfc ~ siRNA * BV2.cell.line + (1 | batch), data = norm_nfkb %>% filter(tx == 'LPS'))
# fit <- lm(lfc ~ siRNA * BV2.cell.line, data = norm_nfkb %>% filter(tx == 'none', !is.na(lfc)))

emm  <- emmeans::emmeans(fit, ~ siRNA | BV2.cell.line)
ctr  <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = "Control", adjust = "none")

tidy_ctr <- broom::tidy(ctr, conf.int = TRUE) |>
  mutate(
    p.adj.BH   = p.adjust(p.value, method = "BH"),          # FDR (Benjamini-Hochberg)
    p.adj.bonf = p.adjust(p.value, method = "bonferroni")   # Bonferroni (optional)
  )

tidy_ctr_wLPS <- tidy_ctr %>% mutate(tx = 'LPS')

# # model diagnostics
# p <- diagnostic_plots(fit)

nfkb_lmm = bind_rows(
    tidy_ctr_noLPS %>% mutate(LPS.dose = '0 ng/ml'),
    tidy_ctr_wLPS %>% mutate(LPS.dose = '100 ng/ml')
  ) %>% 
  mutate(siRNA = str_split_fixed(contrast, ' - ', 2)[,1]) %>%
  left_join(norm_nfkb, ., by = c('siRNA','BV2.cell.line','LPS.dose','tx')) %>%
  mutate(pheno = str_c('nfkb_',tx)) %>% relocate(pheno)


# join --------------------------------------------------------------------

biohyv.lmm <- bind_rows(
  
  mitotracker_lmm %>% 
    select(pheno, line = BV2.cell.line, target = siRNA, lfc = estimate, 
           conf.low, conf.high, p.value, p.adj.BH) %>% 
    distinct(),
  
  ablue_lmm %>% 
    select(pheno, line = BV2.cell.line, target = siRNA, lfc = estimate, 
           conf.low, conf.high, p.value, p.adj.BH) %>% 
    distinct(),
  
  phagocytosis_lmm %>% 
    select(pheno, line = BV2.cell.line, target = siRNA, lfc = estimate, 
           conf.low, conf.high, p.value, p.adj.BH) %>% 
    distinct(),
  
  nfkb_lmm %>%
    select(pheno, line = BV2.cell.line, target = siRNA, lfc = estimate,
           conf.low, conf.high, p.value, p.adj.BH) %>%
    distinct()

  # nfkb_nonpar %>%
  #   select(pheno, line = BV2.cell.line, target = siRNA, lfc,
  #          conf.low = log_lower, conf.high = log_upper,
  #          p.value) %>%
  #   mutate(p.adj.BH = p.value) %>%
  #   distinct()

  ) %>% 
  mutate(
    hit_class = case_when(
      # p.value <= 0.05 & lfc > log2(1.05) ~ 'Positive Hit', #
      # p.value <= 0.05 & lfc < log2(1/1.05) ~ 'Negative Hit', #
      # p.adj.BH <= 0.05 & lfc > log2(1.25) ~ 'Positive Hit', #
      # p.adj.BH <= 0.05 & lfc < log2(1/1.25) ~ 'Negative Hit', #
      p.adj.BH <= 0.05 & lfc > 0 ~ 'Positive Hit', #
      p.adj.BH <= 0.05 & lfc < 0 ~ 'Negative Hit', #
      
      T ~ 'Non-Hit',
    )
  ) %>% 
  filter(target != 'Control')

# tmp <- biohyv.lmm %>% 
#   filter(pheno != 'caspase') %>%
#   left_join(., biohyv_targets) %>% 
#   mutate(hyp_area = str_split(hyp_area, ' \\| ')) %>% 
#   unnest(cols = c(hyp_area)) %>% 
#   mutate(
#     line = factor(line, c('scramble','psen2_kd')),
#     pheno = case_when(pheno == 'nfkb_LPS' ~ 'NFkB +LPS',
#                       pheno == 'nfkb_none' ~ 'NFkB -LPS',
#                       pheno == 'alamar_blue' ~ 'alamarBlue',
#                       pheno == 'mitotracker' ~ 'MitoTracker',
#                       T ~ pheno),
#     pheno = factor(pheno, c('alamarBlue','MitoTracker',
#                             'NFkB -LPS', 'NFkB +LPS', 'caspase', 'phagocytosis') )
#   ) 
# 
# p <- ggplot(tmp, aes( line,target )) + 
#   geom_point( aes(fill = lfc, size = -log10(p.value), shape = hit_class
#                   , color = hit_class))+ #, alpha = hit_class
#   facet_grid(cols = vars(pheno), rows = vars(hyp_area), scales = 'free', space = 'free')+
#   theme( axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1) )+
#   guides( alpha = 'none',  color = 'none', 
#           size = guide_legend( override.aes=list(pch = 22)) ) +
#   scale_y_discrete(limits = rev)+ 
#   scale_size_continuous(guide = 'none')+ 
#   scale_color_manual(values = c('grey20','grey50','grey20'))+
#   scale_alpha_manual(values = c(1,1,1))+
#   scale_shape_manual('Hit Class', values = c(24,22,25), 
#                      breaks = c('Positive Hit','Non-Hit','Negative Hit') )+
#   scale_fill_gradient2('Effect Size\nsiRNA v Cont', 
#                        low = 'purple', mid = 'white', high = 'green', midpoint = 0)+
#   labs(x = '', y = '' 
#        # , title = 'BV2 siRNA phenotypes'
#        # , subtitle = 'significance assessed by Wilcoxon rank sum test'
#   )
# print(p)


write_tsv(biohyv.lmm, here::here('manuscript','data','bv2_assay_results_processed_lmm.tsv'))
save.image(here::here('manuscript','data','bv2_assay_data.rdata'))
