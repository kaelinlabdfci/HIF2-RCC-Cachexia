source('/xchip/beroukhimlab/kei/project/PTHLH_Muhannad/src/00_config.R')

library(tidyverse)
library(here)

clinical_df <- readRDS(here('output/01_makeSIF', 'clinical_df.rds'))

sif <- clinical_df %>%
  select(bcr_patient_barcode, project) %>%
  mutate(tumor_type=sub('^TCGA-', '', project))

## CNA data
cna_pthlh <- readRDS(here('output/03_PTHLH_CNV_vs_tumor_types', 'cna_pthlh.rds')) %>%
  mutate(SampleID=gsub('\\.', '-', barcode)) %>%
  mutate(PatientID=str_extract(SampleID, 'TCGA-.{2}-.{4}')) %>%
  mutate(CNV=case_when(cn==-2 ~ 'Homozygous deletion',
                        cn==-1 ~ 'Heterozygous deletion',
                        cn==0 ~ 'Diploid',
                        cn==1 ~ 'Gain',
                        cn==2 ~ 'Amplification')) %>%
  mutate(CNV=factor(.$CNV, levels=c('Amplification', 'Gain', 'Diploid', 'Heterozygous deletion', 'Homozygous deletion'))) %>%
  mutate(amp_del=case_when(cn %in% c(-2, -1) ~ 'Del',
                           cn %in% c(1, 2) ~ 'Amp',
                           cn==0 ~ 'Diploid')) %>%
  mutate(amp_del=factor(.$amp_del, levels=c('Del', 'Diploid', 'Amp')))

## Gene expression data
tpm <- readRDS(here('../Tangent/20221215_TCGA_RNAseq/output', 'tpm.all.RData')) %>%
  rownames_to_column('gene') %>%
  separate(gene, into=c('ensemblid', 'gene_symbol', 'gene_type', 'egid', 'chr'), sep='\\|') %>%
  mutate(egid=as.numeric(egid))

tpm_pthlh <- tpm %>%
  filter(gene_symbol=='PTHLH') %>%
  pivot_longer(names_to='id', values_to='tpm', cols=-c('ensemblid', 'gene_symbol', 'gene_type', 'egid', 'chr')) %>%
  mutate(SampleID=str_extract(id, '^TCGA-.{2}-.{4}-[0-9]{2}')) %>%
  group_by(SampleID, gene_symbol, egid) %>%
  summarize(tpm=mean(tpm)) %>%
  mutate(log2tpm=log2(tpm + 0.01)) %>%
  mutate(sample_type_code=str_extract(SampleID, '[0-9]{2}$')) %>%
  mutate(sample_type=case_when(sample_type_code=='01' ~ 'TP',
                                sample_type_code=='02' ~ 'TR',
                                sample_type_code=='03' ~ 'TB',
                                sample_type_code=='05' ~ 'TAP',
                                sample_type_code=='06' ~ 'TM',
                                sample_type_code=='07' ~ 'TAM',
                                sample_type_code=='11' ~ 'NT'))

cna_tpm_pthlh <- cna_pthlh %>%
  inner_join(tpm_pthlh, by='SampleID') %>%
  left_join(sif, by=c('PatientID'='bcr_patient_barcode'))

## KIRC
cna_tpm_pthlh_kirc <- cna_tpm_pthlh %>%
  filter(tumor_type=='KIRC') %>%
  filter(Hugo_Symbol=='PTHLH') %>%
  filter(amp_del %in% c('Diploid', 'Amp')) %>%
  group_by(CNV) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  mutate(CNV_n=paste0(CNV, ' (n=', n, ')')) %>%
  mutate(CNV=factor(.$CNV, levels=c('Diploid', 'Gain')))

stat.test <- cna_tpm_pthlh_kirc %>%
  rstatix::wilcox_test(log2tpm ~ CNV, paired=FALSE) %>%
  rstatix::add_xy_position(x='CNV') %>%
  mutate(y.position=y.position + 1)

g <- ggplot(cna_tpm_pthlh_kirc, aes(x=CNV, y=log2tpm)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2) +
  # ggbeeswarm::geom_beeswarm() +
  ggpubr::stat_pvalue_manual(stat.test, label='P = {p}', size=7) +
  scale_x_discrete(breaks=cna_tpm_pthlh_kirc$CNV, labels=cna_tpm_pthlh_kirc$CNV_n) +
  labs(y='Log2(TPM + 0.01)', title='PTHLH mRNA') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('output/05_PTHLH_mRNA_vs_CNV', 'PTHLH_mRNA_vs_CNV_KIRC.png'), dpi=100, width=6, height=8)
ggsave(g, file=here('output/05_PTHLH_mRNA_vs_CNV', 'PTHLH_mRNA_vs_CNV_KIRC.pdf'), width=6, height=8)

write.table(cna_tpm_pthlh_kirc, file=here('output/05_PTHLH_mRNA_vs_CNV', 'table_panelC.csv'), sep=',', row.names=FALSE, quote=FALSE, na='')

stat.test <- cna_tpm_pthlh %>%
  filter(!is.na(tumor_type)) %>%
  group_by(tumor_type) %>%
  rstatix::wilcox_test(log2tpm ~ amp_del, paired=FALSE) %>%
  rstatix::add_xy_position(x='amp_del') %>%
  mutate(y.position=y.position + 1)

## All tumor types
g <- ggplot(cna_tpm_pthlh %>% filter(!is.na(tumor_type)), aes(x=amp_del, y=log2tpm)) +
  geom_boxplot(aes(col=amp_del), outlier.shape=NA, show.legend=FALSE) +
  # geom_point(aes(col=CNV), position=position_jitterdodge(jitter.width=0.2)) +
  ggbeeswarm::geom_beeswarm(aes(fill=CNV), shape=21, stroke=0.1) +
  scale_color_manual(values=c('Del'='blue', 'Diploid'='black', 'Amp'='red')) +
  scale_fill_manual(values=c('Amplification'='#D7191C', 'Gain'='#FDAE61', 'Diploid'='#FFFFBF', 'Heterozygous deletion'='#ABD9E9', 'Homozygous deletion'='#2C7BB6')) +
  ggpubr::stat_pvalue_manual(stat.test %>% filter(p >= 0.05), label='P={p}', size=4) +
  ggpubr::stat_pvalue_manual(stat.test %>% filter(p < 0.05), label='P={p}', size=4, col='red') +
  # facet_grid(~tumor_type, scales='free_x', space='free_x') +
  facet_wrap(~tumor_type, nrow=3) +
  labs(y='Log2(TPM + 0.01)', fill='PTHLH CNV', title='PTHLH mRNA') +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('output/05_PTHLH_mRNA_vs_CNV', 'PTHLH_mRNA_vs_CNV_all_tumor_type.png'), dpi=100, width=18, height=12)
ggsave(g, file=here('output/05_PTHLH_mRNA_vs_CNV', 'PTHLH_mRNA_vs_CNV_all_tumor_type.pdf'), width=18, height=12)
