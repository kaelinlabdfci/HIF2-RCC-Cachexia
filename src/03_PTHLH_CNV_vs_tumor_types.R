source('/xchip/beroukhimlab/kei/project/PTHLH_Muhannad/src/00_config.R')

library(tidyverse)
library(here)

clinical_df <- readRDS(here('output/01_makeSIF', 'clinical_df.rds'))

sif <- clinical_df %>%
  select(bcr_patient_barcode, project) %>%
  mutate(tumor_type=sub('TCGA-', '', project))

cna_files <- list.files(here('data/downloadedFromcBioPortal'), pattern='data_cna.txt', recursive=TRUE, full.names=TRUE)

for (i in seq_along(cna_files)) {
  cna_file_i <- cna_files[i]
  print(paste(i, dirname(cna_file_i)))

  cna_pthlh_i <- read.delim(cna_file_i) %>%
    filter(Hugo_Symbol=='PTHLH' & Entrez_Gene_Id==5744) %>%
    select(c('Hugo_Symbol', starts_with('TCGA'))) %>%
    pivot_longer(names_to='barcode', values_to='cn', cols=-Hugo_Symbol)

  if (i==1) {
    cna_pthlh <- cna_pthlh_i
  } else {
    cna_pthlh <- cna_pthlh %>% bind_rows(cna_pthlh_i)
  }
}
saveRDS(cna_pthlh, here('output/03_PTHLH_CNV_vs_tumor_types', 'cna_pthlh.rds'), compress=FALSE)
# cna_pthlh <- readRDS(here('output/03_PTHLH_CNV_vs_tumor_types', 'cna_pthlh.rds'))

cna_pthlh2 <- cna_pthlh %>%
  mutate(SampleID=sub('\\.[0-9]{2}$', '', barcode)) %>%
  mutate(SampleID=gsub('\\.', '-', SampleID)) %>%
  mutate(sample_type_code=str_extract(barcode, '[0-9]{2}$')) %>%
  mutate(sample_type=case_when(sample_type_code=='01' ~ 'TP',
                                sample_type_code=='03' ~ 'TB',
                                sample_type_code=='06' ~ 'TM')) %>%
  left_join(sif, by=c('SampleID'='bcr_patient_barcode')) %>%
  mutate(CNV=case_when(cn==-2 ~ 'Homozygous deletion',
                        cn==-1 ~ 'Heterozygous deletion',
                        cn==0 ~ 'Diploid',
                        cn==1 ~ 'Gain',
                        cn==2 ~ 'Amplification')) %>%
  mutate(CNV=factor(.$CNV, levels=c('Amplification', 'Gain', 'Diploid', 'Heterozygous deletion', 'Homozygous deletion'))) %>%
  mutate(amp.del=case_when(cn %in% c(-2, -1) ~ 'Del',
                           cn %in% c(1, 2) ~ 'Amp',
                           cn==0 ~ 'Diploid')) %>%
  mutate(amp.del=factor(.$amp.del, levels=c('Amp', 'Diploid', 'Del')))

cna_pthlh_summary <- cna_pthlh2 %>%
  filter(sample_type=='TP') %>%
  filter(!is.na(tumor_type)) %>%
  group_by(tumor_type, CNV) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from=CNV, values_from=n) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(names_to='CNV', values_to='n', cols=-tumor_type) %>%
  group_by(tumor_type) %>%
  mutate(total=sum(n)) %>%
  mutate(amplification_n=n[CNV=='Amplification']) %>%
  mutate(gain_n=n[CNV=='Gain']) %>%
  mutate(amp_total=amplification_n + gain_n) %>%
  mutate(amp_total_ratio=amp_total / total) %>%
  ungroup() %>%
  mutate(rate=n/total) %>%
  mutate(CNV=factor(.$CNV, levels=c('Amplification', 'Gain', 'Diploid', 'Heterozygous deletion', 'Homozygous deletion'))) %>%
  arrange(desc(amp_total_ratio), tumor_type, CNV) %>%
  mutate(tumor_type=factor(.$tumor_type, levels=.$tumor_type %>% unique())) %>%
  mutate(tumor_type_n=paste0(tumor_type, ' (n=', total, ')'))

write.table(cna_pthlh_summary, file=here('output/03_PTHLH_CNV_vs_tumor_types', 'table_panelB.csv'), sep=',', row.names=FALSE, quote=FALSE, na='')

g <- ggplot(cna_pthlh_summary, aes(x=tumor_type, y=rate)) +
  geom_bar(aes(fill=CNV), stat='identity') +
  scale_fill_manual(values=c('Amplification'='#D7191C', 'Gain'='#FDAE61', 'Diploid'='#FFFFBF', 'Heterozygous deletion'='#ABD9E9', 'Homozygous deletion'='#2C7BB6')) +
  scale_x_discrete(breaks=cna_pthlh_summary$tumor_type, labels=cna_pthlh_summary$tumor_type_n) +
  labs(y='Proportion', title='PTHLH CNV') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('output/03_PTHLH_CNV_vs_tumor_types', 'PTHLH_CNV_vs_tumor_types.png'), dpi=100, width=16, height=8)
ggsave(g, file=here('output/03_PTHLH_CNV_vs_tumor_types', 'PTHLH_CNV_vs_tumor_types.pdf'), width=16, height=8)
