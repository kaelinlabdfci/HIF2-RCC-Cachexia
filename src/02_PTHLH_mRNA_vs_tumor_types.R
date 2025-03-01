source('/xchip/beroukhimlab/kei/project/PTHLH_Muhannad/src/00_config.R')

library(tidyverse)
library(here)

clinical_df <- readRDS(here('output/01_makeSIF', 'clinical_df.rds'))

sif <- clinical_df %>%
  select(bcr_patient_barcode, project) %>%
  mutate(tumor_type=sub('TCGA-', '', project))

tpm <- readRDS(here('../Tangent/20221215_TCGA_RNAseq/output', 'tpm.all.RData'))
tpm_pthlh <- tpm %>%
  rownames_to_column('gene') %>%
  filter(grepl('\\|PTHLH\\|', gene)) %>%
  pivot_longer(names_to='id', values_to='tpm', cols=-gene) %>%
  mutate(log2tpm=log2(tpm + 0.01)) %>%
  mutate(barcode=str_extract(id, '^TCGA-.{2}-.{4}-[0-9]{2}')) %>%
  separate(barcode, into=c('SampleID', 'sample_type_code'), sep=-3) %>%
  mutate(sample_type_code=sub('-', '', sample_type_code)) %>%
  mutate(sample_type=case_when(sample_type_code=='01' ~ 'TP',
                                sample_type_code=='02' ~ 'TR',
                                sample_type_code=='03' ~ 'TB',
                                sample_type_code=='05' ~ 'TAP',
                                sample_type_code=='06' ~ 'TM',
                                sample_type_code=='07' ~ 'TAM',
                                sample_type_code=='11' ~ 'NT')) %>%
  filter(sample_type=='TP') %>%
  left_join(sif, by=c('SampleID'='bcr_patient_barcode')) %>%
  filter(!is.na(project)) %>%
  group_by(project) %>%
  mutate(median_log2tpm=median(log2tpm)) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  arrange(desc(median_log2tpm)) %>%
  mutate(project=factor(.$project, levels=.$project %>% unique())) %>%
  mutate(tumor_type_n=paste0(tumor_type, ' (n=', n, ')'))

write.table(tpm_pthlh, file=here('output/02_PTHLH_mRNA_vs_tumor_types', 'table_panelA.csv'), sep=',', row.names=FALSE, quote=FALSE, na='')

g <- ggplot(tpm_pthlh, aes(x=project, y=log2tpm)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2) +
  scale_x_discrete(breaks=tpm_pthlh$project, labels=tpm_pthlh$tumor_type_n) +
  labs(y='log2(TPM + 0.01)', title='PTHLH mRNA expression') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('output/02_PTHLH_mRNA_vs_tumor_types', 'PTHLH_mRNA_vs_tumor_types.png'), dpi=100, width=16, height=8)
ggsave(g, file=here('output/02_PTHLH_mRNA_vs_tumor_types', 'PTHLH_mRNA_vs_tumor_types.pdf'), width=16, height=8)
