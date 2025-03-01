source('/xchip/beroukhimlab/kei/project/PTHLH_Muhannad/src/00_config.R')

library(tidyverse)
library(here)

clinical_df <- readRDS(here('output/01_makeSIF', 'clinical_df.rds'))

sif <- clinical_df %>%
  select(bcr_patient_barcode, project)

## Genomic position data
gene_info_file <- here('../RCC_Nitin/20220906_TCGA_KIRC/data', 'gene_result_chr12.txt')
gene_info <- read.delim(gene_info_file)
gene_info_chr12p <- gene_info %>%
  filter(grepl('12p', map_location)) %>%
  select(GeneID, Symbol, chromosome, map_location, start_position_on_the_genomic_accession, end_position_on_the_genomic_accession) %>%
  rename(start='start_position_on_the_genomic_accession', end='end_position_on_the_genomic_accession')

## CNA data
kirc_cna_file <- here('data/downloadedFromcBioPortal/kirc_tcga_pan_can_atlas_2018', 'data_cna.txt')
kirc_cna <- read.delim(kirc_cna_file)

kirc_cna_chr12p <- kirc_cna %>%
  filter(Entrez_Gene_Id %in% gene_info_chr12p$GeneID) %>%
  pivot_longer(names_to='SampleID', values_to='cn', cols=-c('Hugo_Symbol', 'Entrez_Gene_Id')) %>%
  mutate(SampleID=gsub('\\.', '-', SampleID))

## Gene expression data
tpm <- readRDS(here('../Tangent/20221215_TCGA_RNAseq/output', 'tpm.all.RData')) %>%
  rownames_to_column('gene') %>%
  separate(gene, into=c('ensemblid', 'gene_symbol', 'gene_type', 'egid', 'chr'), sep='\\|') %>%
  mutate(egid=as.numeric(egid))

tpm_chr12p <- tpm %>%
  filter(egid %in% gene_info_chr12p$GeneID) %>%
  pivot_longer(names_to='id', values_to='tpm', cols=-c('ensemblid', 'gene_symbol', 'gene_type', 'egid', 'chr')) %>%
  mutate(log2tpm=log2(tpm + 0.01)) %>%
  mutate(SampleID=str_extract(id, '^TCGA-.{2}-.{4}-[0-9]{2}')) %>%
  mutate(sample_type_code=str_extract(SampleID, '[0-9]{2}$')) %>%
  mutate(sample_type=case_when(sample_type_code=='01' ~ 'TP',
                                sample_type_code=='02' ~ 'TR',
                                sample_type_code=='03' ~ 'TB',
                                sample_type_code=='05' ~ 'TAP',
                                sample_type_code=='06' ~ 'TM',
                                sample_type_code=='07' ~ 'TAM',
                                sample_type_code=='11' ~ 'NT'))

actb_tpm <- tpm %>%
  filter(egid==60) %>%
  select(-c('ensemblid', 'gene_symbol', 'gene_type', 'egid', 'chr')) %>%
  pivot_longer(names_to='id', values_to='actb_tpm', cols=everything())

cna_tpm <- kirc_cna_chr12p %>%
  inner_join(tpm_chr12p, by=c('SampleID', 'Entrez_Gene_Id'='egid', 'Hugo_Symbol'='gene_symbol')) %>%
  left_join(actb_tpm, by='id') %>%
  mutate(rel_exp=tpm / actb_tpm) %>%
  mutate(CNV=case_when(cn==-2 ~ 'Homozygous deletion',
                        cn==-1 ~ 'Heterozygous deletion',
                        cn==0 ~ 'Diploid',
                        cn==1 ~ 'Gain',
                        cn==2 ~ 'Amplification')) %>%
  mutate(CNV=factor(.$CNV, levels=c('Amplification', 'Gain', 'Diploid', 'Heterozygous deletion', 'Homozygous deletion'))) %>%
  mutate(amp_del=case_when(cn %in% c(-2, -1) ~ 'Del',
                           cn %in% c(1, 2) ~ 'Amp',
                           cn==0 ~ 'Diploid')) %>%
  mutate(amp_del=factor(.$amp_del, levels=c('Amp', 'Diploid', 'Del')))

genes <- cna_tpm$Hugo_Symbol %>% unique()
for (i in 1:length(genes)) {
  gene_i <- genes[i]
  egid_i <- cna_tpm %>%
    filter(Hugo_Symbol==gene_i) %>%
    pull(Entrez_Gene_Id) %>%
    unique()

  print(paste(i, gene_i, egid_i))
  
  df_i <- cna_tpm %>%
    filter(Hugo_Symbol==gene_i & Entrez_Gene_Id==egid_i)
  
  dip_exp <- df_i %>%
    filter(amp_del=='Diploid') %>%
    pull(tpm)
  
  amp_exp <- df_i %>%
    filter(amp_del=='Amp') %>%
    pull(tpm)

  fc <- mean(amp_exp) / mean(dip_exp)    

  if (length(dip_exp)>=2 & length(amp_exp)>=2) {
    welch_ttest <- t.test(x=dip_exp, y=amp_exp, var.equal=FALSE, paired=FALSE)
    student_ttest <- t.test(x=dip_exp, y=amp_exp, var.equal=TRUE, paired=FALSE)
    wilcoxon_test <- wilcox.test(x=dip_exp, y=amp_exp, paired=FALSE)
    welch_pval <- welch_ttest$p.value
    student_pval <- student_ttest$p.value
    wilcoxon_pval <- wilcoxon_test$p.value
  } else {
    welch_pval <- NA
    student_pval <- NA
    wilcoxon_pval <- NA
  }

  rel_exp <- df_i %>%
    filter(amp_del=='Amp') %>%
    pull(rel_exp) %>%
    mean()

  stat_i <- data.frame(
    gene=gene_i,
    egid=egid_i,
    gene_type=df_i$gene_type %>% unique(),
    homo_del=df_i %>% filter(CNV=='Homozygous deletion') %>% nrow(),
    hetero_del=df_i %>% filter(CNV=='Heterozygous deletion') %>% nrow(),
    diploid=df_i %>% filter(CNV=='Diploid') %>% nrow(),
    gain=df_i %>% filter(CNV=='Gain') %>% nrow(),
    amp=df_i %>% filter(CNV=='Amplification') %>% nrow(),
    dip_exp_mean=mean(dip_exp),
    amp_exp_mean=mean(amp_exp),
    fc=fc,
    log2fc=log2(fc),
    welch_pval=welch_pval,
    student_pval=student_pval,
    wilcoxon_pval=wilcoxon_pval,
    rel_exp=rel_exp
  )

  if (i==1) {
    stat <- stat_i
  } else {
    stat <- stat %>%
      bind_rows(stat_i)
  }
}

stat_info <- stat %>%
  filter(dip_exp_mean > 1 | amp_exp_mean > 1) %>%
  # mutate(welch_qval=p.adjust(welch_pval, method='BH')) %>%
  # mutate(student_qval=p.adjust(student_pval, method='BH')) %>%
  mutate(wilcoxon_qval=p.adjust(wilcoxon_pval, method='BH')) %>%
  left_join(gene_info_chr12p %>% select(GeneID, start, end, map_location), by=c('egid'='GeneID')) %>%
  mutate(Gene_type=case_when(gene_type %in% c('protein_coding', 'lncRNA', 'miRNA') ~ gene_type, TRUE ~ 'others')) %>%
  mutate(Gene_type=factor(.$Gene_type, levels=c('protein_coding', 'miRNA', 'lncRNA', 'others')))

write.table(stat_info, file=here('output/04_FC_on_chr12p_KIRC', 'table_panelD.csv'), sep=',', row.names=FALSE, quote=FALSE, na='')

hg38 <- rCGH::hg38 %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom)))

centromere_start <- hg38 %>%
  filter(chrom==12) %>%
  pull(centromerStart)

centromere_end <- hg38 %>%
  filter(chrom==12) %>%
  pull(centromerEnd)

pos_start <- 0
pos_end <- round(max(stat_info$start, na.rm=TRUE)/1e+6, digits=-1)

g <- ggplot(stat_info %>%
      filter(!is.na(fc)) %>%
      filter(fc>=1) %>%
      filter(Gene_type=='protein_coding') %>%
      filter(!is.na(start) | !is.na(end)),
    aes(x=start, y=log2fc, label=gene)) +
  geom_hline(yintercept=log2(1.5), linetype='dashed') +
  geom_point(aes(size=rel_exp, fill=-log10(wilcoxon_pval), shape=Gene_type)) +
  ggrepel::geom_text_repel(data=. %>% filter(fc > 1.5 & wilcoxon_pval < 0.05), aes(label=gene), col='blue', size=5, box.padding=1, max.overlaps=Inf) +
  geom_point(data=stat_info %>% filter(gene=="KRAS"), aes(size=rel_exp, fill=-log10(wilcoxon_pval), shape=Gene_type)) +
  ggrepel::geom_text_repel(data=stat_info %>% filter(gene=="KRAS"), aes(label=gene), col='blue', size=5, box.padding=1, max.overlaps=Inf) +
  scale_shape_manual(values=c(21, 22, 23, 24)) +
  scale_size_continuous(range=c(3, 9)) +
  scale_fill_gradient2(low='yellow', mid='white', high='red', midpoint=-log10(0.05)) +
  scale_x_continuous(breaks=seq(pos_start*1e+06, pos_end*1e+06, 10e+06), labels=paste0(seq(pos_start*1e+06, pos_end*1e+06, 10e+06)/1e+06, 'Mb')) +
  labs(x='Genomic position on chr12p', y='log2FC (Amp/Diploid)', size='RelExp (/ACTB)', fill='-log10pvalue', shape='Gene type') +
  guides(fill=guide_colorbar(order=1), size=guide_legend(order=2), shape=guide_legend(order=3)) +
  theme_classic(base_size=20) +
  theme(plot.margin = unit(c(2,1,1,1), 'lines'))
ggsave(g, file=here('output/04_FC_on_chr12p_KIRC', 'FC_vs_position_KIRC.png'), dpi=100, width=16, height=8)
ggsave(g, file=here('output/04_FC_on_chr12p_KIRC', 'FC_vs_position_KIRC.pdf'), width=16, height=8)
