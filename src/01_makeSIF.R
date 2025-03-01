source('/xchip/beroukhimlab/kei/project/PTHLH_Muhannad/src/00_config.R')

library(tidyverse)
library(here)

library(TCGAbiolinks)

## Retrieve clinical information including patient's weight and height
tcga_projects <- c('TCGA-ACC',
                    'TCGA-BLCA',
                    'TCGA-BRCA',
                    'TCGA-CESC',
                    'TCGA-CHOL',
                    'TCGA-COAD',
                    'TCGA-DLBC',
                    'TCGA-ESCA',
                    'TCGA-GBM',
                    'TCGA-HNSC',
                    'TCGA-KICH',
                    'TCGA-KIRC',
                    'TCGA-KIRP',
                    'TCGA-LAML',
                    'TCGA-LGG',
                    'TCGA-LIHC',
                    'TCGA-LUAD',
                    'TCGA-LUSC',
                    'TCGA-MESO',
                    'TCGA-OV',
                    'TCGA-PAAD',
                    'TCGA-PCPG',
                    'TCGA-PRAD',
                    'TCGA-READ',
                    'TCGA-SARC',
                    'TCGA-SKCM',
                    'TCGA-STAD',
                    'TCGA-TGCT',
                    'TCGA-THCA',
                    'TCGA-THYM',
                    'TCGA-UCEC',
                    'TCGA-UCS',
                    'TCGA-UVM')

retrieve_clinical_data <- function(tcga_project) {
  query <- GDCquery(
    project = tcga_project, 
    data.category = "Clinical", 
    data.format = "bcr xml"
  )
  GDCdownload(query)
  clinical <- GDCprepare_clinic(query, clinical.info = "patient")
  return(clinical)
}

# clinical_list <- lapply(tcga_projects, retrieve_clinical_data)
# clinical_df <- clinical_list %>%
#   bind_rows()

for (i in seq_along(tcga_projects)) {
for (i in 24:length(tcga_projects)) {
  tcga_project_i <- tcga_projects[i]

  clinical_i <- retrieve_clinical_data(tcga_project_i) %>%
    mutate(across(everything(), as.character))

  print(paste(i, tcga_project_i, nrow(clinical_i)))

  if (i==1) {
    clinical_df <- clinical_i
  } else {
    clinical_df <- clinical_df %>%
      bind_rows(clinical_i)
  }
}
clinical_df <- clinical_df %>% type.convert(as.is=TRUE)
saveRDS(clinical_df, here('output/01_makeSIF', 'clinical_df.rds'), compress=FALSE)

