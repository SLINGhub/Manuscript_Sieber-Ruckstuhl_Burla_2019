##################################################################################
# Misc functions
#
# 20.01.2019 / Bo Burla 
# Singapore Lipidomics Incubator (SLING), National University of Singapore
##################################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)

# Reformat lipid names
clean_lipidnames <- function(datLong){
  datLong_clean <- datLong %>% 
    mutate(Compound= str_replace(Compound, "\\/C(?=\\d)","\\/")) %>% 
    mutate(Compound= str_replace(Compound, "[:space:]C"," "))
  return(datLong_clean)  
}

# Retreive lipid class from lipid name
add_lipidclass_names <- function(datLong){
  datLong_temp <- clean_lipidnames(datLong)
  datLong_temp <- datLong_temp %>%  
                    mutate(lipidClassBase = str_trim(str_extract(Compound, "[A-z0-9]+[[:blank:]]*")),
                           lipidClass = str_trim(str_extract(Compound, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m][0-9]{2}[:]{1}[0-9]{1})")))
  
  datLong_temp <- datLong_temp %>% 
                    mutate(lipidClass = str_replace(lipidClass, "([[:blank:]]+)([^d|t|m]{1})", "-\\2"),
                           isPUFA = ifelse(str_detect(Compound, "\\:0|\\:1|\\:2|\\:3"),FALSE,TRUE)) %>% 
                    select(-Compound, -Conc, Compound, Conc)
  
  datLong_temp <- datLong_temp %>% mutate_at(.vars = c("lipidClass", "lipidClassBase"), as.factor)
  return(datLong_temp)
}

# Use within stat_summary function to plot median lines
plot.mean <- function(x) {
  m <- mean(x)
  c(y = m, ymin = m, ymax = m)
}

# ggplot theme used for PCA plots
mytheme <- theme_pubr() +
  theme(axis.line = element_line(colour = 'black',size=0.8),  
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10),
        legend.position="right")
