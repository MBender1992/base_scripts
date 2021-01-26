
##################################################
# This script contains several useful functions  #
# for the analysis of chronic irradiated SCC     #
# cell lines for my PhD                          #
##################################################




######################################################
#                                                    #
# 1. Functions to ingest data in a reproducible way  #
#                                                    #
######################################################

load_Fireplex_data_PhD <- function(filename = "200619_chronic_irr_normalized.csv", threshold)
{
  # check if required package is installed and load it
  if (!require(dplyr)) install.packages("dplyr")
  library(dplyr)
  
  #load data and change hyphens to underscores so r does not confuse a hyphen with minus
  dat <- read_csv(filename)
  names(dat) <- str_replace_all(names(dat), "-","_")
  # clean up the data, convert to tidy format and replace negative expression values with 0 as expression cannot be negative
  dat_f <- dat %>% 
    select(-Messung) %>%
    filter(cell_line != "HaCaT") %>%
    gather(miRNA, expression, -c(ID,cell_line, Irradiation)) %>% 
    mutate(expression = ifelse(expression < 0, 0, expression))  
  # index miRNAs with a median expression below the threshold
  ind <- dat_f %>%
    group_by(miRNA) %>%
    filter(median(expression, na.rm=TRUE) <= threshold) %>%
    .$miRNA %>% 
    unique() 
  # omit indexed miRNAs 
  dat_f %>% 
    filter(!miRNA %in% ind) %>%
    mutate(log_exp = log2(expression+1)) %>%
    mutate(miRNA = str_replace_all(miRNA, "hsa_", "")) %>%
    mutate(miRNA = str_replace_all(miRNA, "_","-")) %>%
    mutate(cell_line = str_replace_all(cell_line, "-","_"))
}




######################################################
#                                                    #
#        3. Plotting and significance Tests          #
#                                                    #
######################################################


# ..................................................................................................................
# this function is used to ensure a uniform theme and layout is used in each plot #

# Arguments                                                     #
#  - axis.text.size: Size of the axes                           #
#  - Legend: logical indicating if legend should be printed     #
theme_PhD <- function(axis.text.size=10, Legend = TRUE,...){
  
  theme_custom <- theme_pubr()+
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(face = "bold", size=axis.text.size),
          axis.text.y = element_text(face = "bold", size=axis.text.size),
          axis.title.y = element_text(face = "bold", size = axis.text.size+1),
          panel.grid.minor=element_blank(),
          strip.text.x = element_text(face = "bold", size = axis.text.size+1),
          strip.background=element_blank(),
          panel.spacing = unit(1, "lines"), ...) 
  
  
  if(Legend == TRUE) {
    theme_custom + theme(legend.key.size = unit(1,"line"))
  } else {
    theme_custom + theme(legend.position = "none")
  }
}

# ..................................................................................................................
