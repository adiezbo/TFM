############# CODIGO PARA COAD FILTRANDO DESPUES 450K ###########
```{r COAD} 
#BiocManager::install("GenomicDataCommons")
library(GenomicDataCommons)
library(dplyr)

#Filtramos muestras de metilacion de tejido normal
nl_ge_files_COAD = files() %>%
    GenomicDataCommons::filter(~cases.samples.sample_type=='Solid Tissue Normal' &
               cases.project.project_id == 'TCGA-COAD' &
               data_type == "Methylation Beta Value") %>%
               #platform == "Illumina Human Methylation 450") 
               
    expand(c('cases','cases.samples')) %>%
    results_all() %>%
    as_tibble()

#Filtramos muestras de metilacion de tejido tumoral
tm_ge_files_COAD = files() %>%
    GenomicDataCommons::filter(~cases.samples.sample_type=='Primary Tumor' &
               cases.project.project_id == 'TCGA-COAD' &
               data_type == "Methylation Beta Value") %>%
               #platform == "Illumina Human Methylation 450") 
               
    expand(c('cases','cases.samples')) %>%
    results_all() %>%
    as_tibble()

#Nos quedamos solo con muestras de pacientes con tejido normal y tumoral
nl_cases_COAD = bind_rows(nl_ge_files_COAD$cases, .id='file_id')
tm_cases_COAD = bind_rows(tm_ge_files_COAD$cases, .id='file_id')
matchedcases_COAD = intersect(nl_cases_COAD$case_id, tm_cases_COAD$case_id)
COAD <- length(matchedcases_COAD) 

#Creamos data frame con muestras n y t apareadas
matched_nl_files_COAD = nl_cases_COAD[nl_cases_COAD$case_id %in% matchedcases_COAD, 'file_id']
matched_tm_files_COAD = tm_cases_COAD[tm_cases_COAD$case_id %in% matchedcases_COAD, 'file_id']

matched_tn_ge_file_info_COAD = rbind(subset(nl_ge_files_COAD,file_id %in% matched_nl_files_COAD),
                                subset(tm_ge_files_COAD,file_id %in% matched_tm_files_COAD))
head(matched_tn_ge_file_info_COAD)

#Nos quedamos solo con muestras de Illumina Human Methylation 450
matched_tn_ge_file_info_COAD <- matched_tn_ge_file_info_COAD[!matched_tn_ge_file_info_COAD$platform == "Illumina Human Methylation 27",]

#Guardamos data frame
muestras_COAD = data.frame(lapply(matched_tn_ge_file_info_COAD, as.character), stringsAsFactors=FALSE)

write.table(muestras_COAD, file = "muestras_COAD.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
```

########################### CODIGO COAD FILTRANDO AL PRINCIPIO" #################################

```{r COAD2} 
#BiocManager::install("GenomicDataCommons")
library(GenomicDataCommons)
library(dplyr)

#Filtramos muestras de metilacion de tejido normal
nl_ge_files_COAD_2 = files() %>%
    GenomicDataCommons::filter(~cases.samples.sample_type=='Solid Tissue Normal' &
               cases.project.project_id == 'TCGA-COAD' &
               #data_type == "Methylation Beta Value") 
               platform == "Illumina Human Methylation 450") %>%  
               
    expand(c('cases','cases.samples')) %>%
    results_all() %>%
    as_tibble()

#Filtramos muestras de metilacion de tejido tumoral
tm_ge_files_COAD_2 = files() %>%
    GenomicDataCommons::filter(~cases.samples.sample_type=='Primary Tumor' &
               cases.project.project_id == 'TCGA-COAD' &
               #data_type == "Methylation Beta Value") 
               platform == "Illumina Human Methylation 450") %>% 
               
    expand(c('cases','cases.samples')) %>%
    results_all() %>%
    as_tibble()

#Nos quedamos solo con muestras de pacientes con tejido normal y tumoral
nl_cases_COAD_2 = bind_rows(nl_ge_files_COAD_2$cases, .id='file_id')
tm_cases_COAD_2 = bind_rows(tm_ge_files_COAD_2$cases, .id='file_id')
matchedcases_COAD_2 = intersect(nl_cases_COAD_2$case_id, tm_cases_COAD_2$case_id)
COAD2 <- length(matchedcases_COAD_2) 

#Creamos data frame con muestras n y t apareadas
matched_nl_files_COAD_2 = nl_cases_COAD_2[nl_cases_COAD_2$case_id %in% matchedcases_COAD_2, 'file_id']
matched_tm_files_COAD_2 = tm_cases_COAD_2[tm_cases_COAD_2$case_id %in% matchedcases_COAD_2, 'file_id']

matched_tn_ge_file_info_COAD_2 = rbind(subset(nl_ge_files_COAD_2,file_id %in% matched_nl_files_COAD_2),
                                subset(tm_ge_files_COAD_2,file_id %in% matched_tm_files_COAD_2))
head(matched_tn_ge_file_info_COAD_2)

#Nos quedamos solo con muestras de Illumina Human Methylation 450
#matched_tn_ge_file_info_COAD <- matched_tn_ge_file_info_COAD[!matched_tn_ge_file_info_COAD$platform == "Illumina Human Methylation 27",]

#Guardamos data frame
muestras_COAD_2 = data.frame(lapply(matched_tn_ge_file_info_COAD_2, as.character), stringsAsFactors=FALSE)
write.table(muestras_COAD_2, file = "muestras_COAD_2.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
```
############## CODIGO LIHC 450K #################################
```{r LIHC} 
#BiocManager::install("GenomicDataCommons")
library(GenomicDataCommons)
library(dplyr)

#Filtramos muestras de metilacion de tejido normal
nl_ge_files_LIHC = files() %>%
    GenomicDataCommons::filter(~cases.samples.sample_type=='Solid Tissue Normal' &
               cases.project.project_id == 'TCGA-LIHC' &
               #data_type == "Methylation Beta Value") %>%
               platform == "Illumina Human Methylation 450") %>%  
    expand(c('cases','cases.samples')) %>%
    results_all() %>%
    as_tibble()

#Filtramos muestras de metilacion de tejido tumoral
tm_ge_files_LIHC = files() %>%
    GenomicDataCommons::filter(~cases.samples.sample_type=='Primary Tumor' &
               cases.project.project_id == 'TCGA-LIHC' &
               #data_type == "Methylation Beta Value") %>%
               platform == "Illumina Human Methylation 450") %>% 
    expand(c('cases','cases.samples')) %>%
    results_all() %>%
    as_tibble()

#Nos quedamos solo con muestras de pacientes con tejido normal y tumoral
nl_cases_LIHC = bind_rows(nl_ge_files_LIHC$cases, .id='file_id')
tm_cases_LIHC = bind_rows(tm_ge_files_LIHC$cases, .id='file_id')
matchedcases_LIHC = intersect(nl_cases_LIHC$case_id, tm_cases_LIHC$case_id)
LIHC <- length(matchedcases_LIHC) 

#Creamos data frame con muestras n y t apareadas
matched_nl_files_LIHC = nl_cases_LIHC[nl_cases_LIHC$case_id %in% matchedcases_LIHC, 'file_id']
matched_tm_files_LIHC = tm_cases_LIHC[tm_cases_LIHC$case_id %in% matchedcases_LIHC, 'file_id']

matched_tn_ge_file_info_LIHC = rbind(subset(nl_ge_files_LIHC,file_id %in% matched_nl_files_LIHC),
                                subset(tm_ge_files_LIHC,file_id %in% matched_tm_files_LIHC))
head(matched_tn_ge_file_info_LIHC)

#Nos quedamos solo con muestras de Illumina Human Methylation 450
#matched_tn_ge_file_info_LIHC <- matched_tn_ge_file_info_LIHC[!matched_tn_ge_file_info_LIHC$platform == "Illumina Human Methylation 27",]

#Guardamos data frame
muestras_LIHC = data.frame(lapply(matched_tn_ge_file_info_LIHC, as.character), stringsAsFactors=FALSE)
#write.csv(muestras_COAD,"tx.csv")
write.table(muestras_LIHC, file = "muestras_LIHC_2.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
```
############## CODIGO LUAD 450K #################################
```{r LUAD} 
#BiocManager::install("GenomicDataCommons")
library(GenomicDataCommons)
library(dplyr)

#Filtramos muestras de metilacion de tejido normal
nl_ge_files_LUAD = files() %>%
    GenomicDataCommons::filter(~cases.samples.sample_type=='Solid Tissue Normal' &
               cases.project.project_id == 'TCGA-LUAD' &
               #data_type == "Methylation Beta Value") %>%
                platform == "Illumina Human Methylation 450") %>%
    expand(c('cases','cases.samples')) %>%
    results_all() %>%
    as_tibble()

#Filtramos muestras de metilacion de tejido tumoral
tm_ge_files_LUAD = files() %>%
    GenomicDataCommons::filter(~cases.samples.sample_type=='Primary Tumor' &
               cases.project.project_id == 'TCGA-LUAD' &
               #data_type == "Methylation Beta Value") %>%
                platform == "Illumina Human Methylation 450") %>%       
    expand(c('cases','cases.samples')) %>%
    results_all() %>%
    as_tibble()

#Nos quedamos solo con muestras de pacientes con tejido normal y tumoral
nl_cases_LUAD = bind_rows(nl_ge_files_LUAD$cases, .id='file_id')
tm_cases_LUAD = bind_rows(tm_ge_files_LUAD$cases, .id='file_id')
matchedcases_LUAD = intersect(nl_cases_LUAD$case_id, tm_cases_LUAD$case_id)
LUAD <- length(matchedcases_LUAD) 

#Creamos data frame con muestras n y t apareadas
matched_nl_files_LUAD = nl_cases_LUAD[nl_cases_LUAD$case_id %in% matchedcases_LUAD, 'file_id']
matched_tm_files_LUAD = tm_cases_LUAD[tm_cases_LUAD$case_id %in% matchedcases_LUAD, 'file_id']

matched_tn_ge_file_info_LUAD = rbind(subset(nl_ge_files_LUAD,file_id %in% matched_nl_files_LUAD),
                                subset(tm_ge_files_LUAD,file_id %in% matched_tm_files_LUAD))
head(matched_tn_ge_file_info_LUAD)

#Nos quedamos solo con muestras de Illumina Human Methylation 450
#matched_tn_ge_file_info_LUAD <- matched_tn_ge_file_info_LUAD[!matched_tn_ge_file_info_LUAD$platform == "Illumina Human Methylation 27",]

#Guardamos data frame
muestras_LUAD = data.frame(lapply(matched_tn_ge_file_info_LUAD, as.character), stringsAsFactors=FALSE)
#write.csv(muestras_COAD,"tx.csv")
write.table(muestras_LUAD, file = "muestras_LUAD.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
```
