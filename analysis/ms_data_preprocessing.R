####################################################################################################
# R script used for the normalization and QC filtering of raw results from the lipidomics analyses  
#
# 08.01.2019 / Bo Burla 
# Singapore Lipidomics Incubator (SLING), National University of Singapore
####################################################################################################

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)


# Function to read, normalize and QC filter a dataset with peak areas
# -------------------------------------------------------------------
  
process_rawdata_panel <- function(data_file, panel_text, ISTDconc_sheet, PeakAreas_sheet, ISTDmap_sheet, 
                                  sample_vol, istd_vol, 
                                  min_pqc_area_treshold, min_sbr_treshold, max_cva_treshold, species_to_keep){ 
  
  # Inport data from Excel workbook
  concISTD <- read_excel(data_file, sheet = ISTDconc_sheet,
                         col_names = TRUE,trim_ws = TRUE,na = c("#N/A", "NULL"), col_types = "guess")
  datWide <- read_excel(data_file, sheet = PeakAreas_sheet,
                        col_names = TRUE,trim_ws = TRUE,na = c("#N/A", "NULL"), col_types = "guess")
  mapISTD <- read_excel(data_file, sheet = ISTDmap_sheet,
                        col_names = TRUE,trim_ws = TRUE,na = c("#N/A", "NULL"), col_types = "guess")
  
  # Transform data table format from wide to long (tidy) and convert type of Area from int to num 
  datLong <- gather(datWide, key="Compound", value="Area", 
                    factor_key=FALSE, -RawDataFileName, -SampleName, -AcqDateTime, 
                    -SampleType, -AnalyticalGroup)
  datLong <- datLong %>% select(-RawDataFileName,-AcqDateTime )
  datLong$Area <- as.numeric(datLong$Area) %>% replace(is.na(.), 0)

  # ISTD normalization and calculation of concentrations
  #  - Joins a column with ISTD for each compound (based on ISTDmap file), groups by ISTD and 
  #    normalizes with the value of the ISTD
  #  - Calculates concentrations (based on ISTD concentrations in the table concISTD and given 
  #    sample and ISTD volumes)

  dat_norm <- datLong  %>% 
    left_join(mapISTD, c("Compound"="Compound")) %>% 
    left_join(concISTD, by=c("ISTD"="ISTD")) %>%
    mutate(isISTD = (Compound == ISTD)) %>% 
    group_by(ISTD,SampleName) %>% 
          mutate(
            normArea = Area*CF/Area[isISTD]*CF[isISTD],
            uM = normArea*(istd_vol*ISTDconcNGML/ISTD_MW)/sample_vol) 

  # Calculate QC metrics (Signal/Blank ratio, analytical CV and median signal of PQCs) for each 
  # transition/compound per experimental group
  cmpd_QC <- dat_norm %>%
    filter(SampleType =="PQC" & AnalyticalGroup != "") %>%
    group_by(Compound, Quantifier, isISTD, AnalyticalGroup, SampleType) %>%
    do(
      summarise(.data = ., 
                SBR = median(Area,na.rm = TRUE)/max(dat_norm$Area[dat_norm$SampleType == "BLANK" & dat_norm$Compound == unique(Compound)],
                                                    na.rm = TRUE),
                CVA = sd(normArea ,na.rm = TRUE)/mean(normArea,na.rm = TRUE)*100,
                INT = median(Area,na.rm = TRUE)))

  # Get list of compound that fullfill the QC criteria passed as parameter of this funcion
  comp_QC_pass <- cmpd_QC %>% ungroup() %>%
                    filter((INT > min_pqc_area_treshold & 
                           SBR > min_sbr_treshold &
                           CVA < max_cva_treshold) | Compound %in% species_to_keep )
  
  cat("Median CVA in PRED Group of ", panel_text, ":" ,
      median(comp_QC_pass %>% filter(AnalyticalGroup=="PRED" & !isISTD & Quantifier & !str_detect(Compound, "\\(IS")) %>% .$CVA),"\n")
  cat("Median CVA in ACTH Group of ", panel_text, ":" ,
      median(comp_QC_pass %>% filter(AnalyticalGroup=="ACTH" & !isISTD & Quantifier & !str_detect(Compound, "\\(IS")) %>% .$CVA),"\n")
  
  # Filter the dataset with the list of compounds that passed QC
  dat_norm_filt <- dat_norm %>% 
    filter(Quantifier, !isISTD, !str_detect(Compound, "(IS)"), SampleType=="SPL") %>%
    semi_join(comp_QC_pass, by=c("AnalyticalGroup"="AnalyticalGroup", "Compound"="Compound"))

  return (dat_norm_filt)
}


# Process raw data from all 4 panels (PL/DG/CE, SL, TG and S1P) and combine
# -------------------------------------------------------------------------

process_all_panels <- function(data_file){
  
  sample_annot <- read_excel(data_file, sheet = "SampleAnnotation",
                             col_names = TRUE,trim_ws = TRUE,na = c("#N/A", "NULL"), col_types = "guess")
  
  d_PL <- process_rawdata_panel(data_file, panel_text = "PL-CE-DG", 
                                ISTDconc_sheet = "ISTDconc",
                                PeakAreas_sheet = "Panel1_PL-CE-DG_PeakAreas",
                                ISTDmap_sheet = "Panel1_PL-CE-DG_ISTDmap",
                          sample_vol = 10, istd_vol = 90,
                          min_pqc_area_treshold = 250,  min_sbr_treshold = 5, max_cva_treshold = 25, 
                          species_to_keep = "")
  
  # REMARK: Cer d18:1/18:0 and Hex1Cer d18:1/24:1 had CVs of 27.8% and 30.1% in the prednisolone group, but 
  #         were kept in the dataset as they belong to reported the CVD markers and their CVs were 7.3% and 
  #         10.9%, respectively, the tetracosactide group. 
  d_SL <- process_rawdata_panel(data_file, panel_text = "SL", 
                                ISTDconc_sheet = "ISTDconc",
                                PeakAreas_sheet = "Panel2_SL_PeakAreas",
                                ISTDmap_sheet = "Panel2_SL_ISTDmap",
                          sample_vol = 10, istd_vol = 90,
                          min_pqc_area_treshold = 250,  min_sbr_treshold = 5, max_cva_treshold = 25,
                          species_to_keep = c("Cer d18:1/18:0", "Hex1Cer d18:1/24:1"))
  
  d_S1P <- process_rawdata_panel(data_file, panel_text = "S1P", 
                                 ISTDconc_sheet = "ISTDconc",
                                 PeakAreas_sheet = "Panel3_S1P_PeakAreas",
                                 ISTDmap_sheet = "Panel3_S1P_ISTDmap",
                          sample_vol = 10, istd_vol = 90,
                          min_pqc_area_treshold = 250,  min_sbr_treshold = 5, max_cva_treshold = 25,
                          species_to_keep = "")
  
  d_TG <- process_rawdata_panel(data_file, panel_text = "TG",
                                ISTDconc_sheet = "ISTDconc",
                                PeakAreas_sheet = "Panel4_TG_PeakAreas",
                                ISTDmap_sheet = "Panel4_TG_ISTDmap",
                          sample_vol = 10, istd_vol = 90,
                          min_pqc_area_treshold = 250,  min_sbr_treshold = 5, max_cva_treshold = 25,
                          species_to_keep = "") 
  
  # Substracting interference from the M+2 isotope that were not sufficientlt separated and 
  # therefore included in the integrated peak
  #   PC O-36:3: rel. abundance of PC-O-36:3 M+2 ion without phosphocholine)
  #   S1P: see Burla et al, Methods Mol Biol., 2018)
  
  d_PL <- d_PL %>% mutate(uM = ifelse(Compound == "PC O-36:3", 
                                      uM[Compound == "PC O-36:3"] - 0.0982 * uM[Compound == "PC O-36:4"], 
                                      uM)) 
  d_S1P <- d_S1P %>% mutate(uM = ifelse(Compound == "S1P d18:0", 
                                        uM[Compound == "S1P d18:0"] - 0.03161 * uM[Compound == "S1P d18:1"], 
                                        uM))  
  
  
  d_final <- bind_rows(d_PL, d_SL, d_TG, d_S1P)
  d_final <- sample_annot %>% left_join(d_final, by=c("SampleName","SampleName"))
  
  cat(paste0("Number of species in final filtered dataset: "),length(unique(d_final$Compound)))
  return(d_final)
}

# Process raw peak area data from all MS analyses (normalize and QC filter)
# ==========================================================================================


process_ms_dataset <- function(rawdata_file){
  d_final <- process_all_panels(data_file = rawdata_file) %>% rename(Conc = uM)   
  d_final_wide_uM <- d_final %>% ungroup() %>% 
    select(SampleName, AnimalID, TreatmentGroup, ExperimentalGroup, TreatmentTimePoint, 
           Sex, Compound, Conc) %>% 
    spread(key = Compound,value = Conc, drop=TRUE) %>% 
    arrange(TreatmentGroup, AnimalID)
  
  return(d_final_wide_uM)
}



