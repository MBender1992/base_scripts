#<<<<<<<<<<<HEAD 

setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/base_scripts")

##################################################
# This script contains several useful functions  #
# for the analysis of melanoma data gathered     #
# for the phd thesis of Ramona/Leonie            #
##################################################

######################################################
#                                                    #
# 1. Functions to ingest data in a reproducible way  #
#                                                    #
######################################################


#....................................................................................................................
# this function allows transposition of dataframes while retaining row and column names #
# (which are lost by default when transposing data as its converted into a matrix)      #

# Arguments                                                                                         #
#  - colnames: names of the dataframe columns are needed as the transpose argument omits colnames   #
#              so they need to be added manually later                                              #
#  - data: dataframe to be transposed                                                               #

transpose_dataframe <-  function(colnames, data){
  names <- c(colnames)
  dat_wide <- data[,-1] %>%
    t %>% 
    as_tibble(rownames = NA) %>%
    rownames_to_column %>% 
    setNames(names)
  return(dat_wide)
}
# ..................................................................................................................






# ....................................................................................................................
# this function uses the previously defined functions and a uniform pipeline #
# to combine data from the 4 tables with raw data and transform them to a    #
# uniform format                                                             #

# Arguments                                                                             #
#  - characterAsFactor: logical defining if characters should be converted to factors   #
#                       important for example for machine learning                      #

load_melanoma_data <- function(characterAsFactor = TRUE){
  
  dir <- getwd()
  
  if(dir != "Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/Projekte/Doktorarbeiten_Melanom_Liquid_Biopsies/Daten"){
    setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/Projekte/Doktorarbeiten_Melanom_Liquid_Biopsies/Daten")
  } 
  
  # load csv files
  dat_miR   <- read_csv("miRNA_Expression_Fireplex_2.csv") 
  dat_add   <- read_csv("Add_Info_grouped_age_sex.csv") 
  dat_res   <- read_csv("Responder_misc.csv") %>% 
    mutate(ID = as.character(ID)) %>% select(-c(therapy_start, Abnahmedatum)) %>%
    mutate(TRIM_PDL1_Expression = str_replace_all(TRIM_PDL1_Expression,"\\++","+")) %>% 
    mutate(TRIM_PDL1_Expression = ifelse(TRIM_PDL1_Expression == "o", NA,TRIM_PDL1_Expression)) %>% 
    mutate(breslow_thickness_mm = str_replace_all(breslow_thickness_mm, " mm", ""), 
           breslow_thickness_mm = parse_number(breslow_thickness_mm)) %>% 
    mutate(subtype = str_replace_all(subtype, "cutanes Melanom", "cutaneous"),
           subtype = str_replace_all(subtype, "Schleimhautmelanom", "mucosal"))
  dat_lab   <- read_csv("Patiententabelle.csv")
  
  # change ID column to uniform capital letters for later filtering
  names(dat_miR) <- c("miRNA", toupper(names(dat_miR)[-1]))
  names(dat_add) <- c("Sample", toupper(names(dat_add)[-1]))
  
  # define IDs to be dropped for further analyses
  controls <- c("K104_1", "K104_2", "K104_3A", "K104_3B")
  duplicates <- c("22B","38B","39B","47B")
  
  # wide miR data
  dat_miR_trans <- transpose_dataframe(colnames = c("ID",dat_miR$miRNA), data = dat_miR) %>%
    filter(!ID %in% controls & !ID %in% duplicates) %>%   #drop duplicate patient data 
    mutate(ID = parse_number(ID)) #convert ID to numeric
  
  
  # wide metadata
  dat_add_trans <- transpose_dataframe(colnames = c("ID",dat_add$Sample), data = dat_add) %>%
    mutate(ID = toupper(ID)) %>%
    filter(!ID %in% controls & !ID %in% duplicates)%>% #drop duplicate patient data
    mutate(ID = parse_number(ID), #convert ID to numeric
           age = parse_number(as.character(age))) %>% #convert age to numeric
    select(-Treatment)
  
  
  #filter out duplicate patients in responder dataset
  dat_res_filter <- dat_res %>%  
    filter(!ID %in% controls & !ID %in% duplicates)%>%
    mutate(ID = parse_number(ID))
  
  setwd(dir)
  
  # convert feature columns to factor for ML purpose
  if(characterAsFactor == TRUE)
  {
    inner_join(dat_miR_trans,dat_add_trans, by="ID") %>%
      inner_join(., dat_res_filter, by="ID") %>%
      left_join(., dat_lab, by =c("ID","age","sex","Responder")) %>%
      select(-c(Immuntherapie))  %>%
      mutate(sex = factor(sex),
             Responder = factor(ifelse(Responder == "ja", "pos","neg")),
             adjuvant_IFN = factor(adjuvant_IFN),
             Stadium = factor(ifelse(Stadium == "IV", "IV","<IV")))%>% # dichotomize outcome of Stadium variable
      # drop unused factor levels
      droplevels()
  }
  
  # keep character columns as character for data gathering
  else{
    inner_join(dat_miR_trans,dat_add_trans, by="ID") %>%
      inner_join(., dat_res_filter, by="ID") %>%
      left_join(., dat_lab, by =c("ID","age","sex","Responder")) %>%
      select(-c(Immuntherapie))
  }
  
  
}
# ..................................................................................................................








# ....................................................................................................................
# this function loads the data for overlap wang to be directly used in a heatmap

# Arguments                                                                             

dat_heatmap_wang <- function(add.data = NULL, join.by = "ID"){
  # change working directory
  setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/Projekte/Doktorarbeiten_Melanom_Liquid_Biopsies/Daten")
  
  # load miRNA data
  dat_melanoma <- load_melanoma_data(characterAsFactor = FALSE) %>% 
    select(-c(adjuvant_IFN, befallen_Organe))
  # load overlap data
  overlap_wang <- read_csv("overlap_our_panel_Wang.csv")
  
  #tidy data
  dat_tidy <- dat_melanoma %>% 
    gather(
      miRNA, 
      expression,
      -c(ID, Hirnmetastase, BRAF, Stadium,
         Group, Baseline,Responder, age, sex,
         CRP, LDH, S100, Eosinophile,TRIM_PDL1_Expression)
    ) %>% 
    mutate(Stadium = ifelse(Stadium == "IV", "IV", "<IV"),
           log_exp = log2(expression),
           BRAF = factor(BRAF, levels = c("wt", "mu", "nras")))
  
  # define index with miRNAs in overlap
  index_wang <- overlap_wang  %>% .$miRNA %>%
    str_replace_all(" ", "")
  
  # extract overlap data
  if(!is.null(add.data)){
  dat_wang <- dat_tidy %>% 
    filter(miRNA %in% index_wang) %>%
    inner_join(add.data, by = join.by)
  } else { 
    dat_wang <- dat_tidy %>% 
      filter(miRNA %in% index_wang)
    }
  
  #  prepare data for Heatmap
  dat_Heatmap <- dat_wang %>% 
    select(-expression) %>%
    spread(miRNA, log_exp) %>% 
    data.frame() %>%
    column_to_rownames("ID") 
  
  return(dat_Heatmap)
}
#.............................................................................................................







######################################################
#                                                    #
#        2. Functions for summary statistics         #
#                                                    #
######################################################


#....................................................................................................................
# this function calculates the proportions of a given characteristic #
# (for example number of patients with BRAF mutations)               #
# also calculates percentages of the whole population                #

# Arguments                                         #
#  - data: dataframe to be analysed                 #
#  - var: variable (column) that shall be analysed  #

patient_table <- function(data, var){
  nas <- table(is.na(data[,var]))
  dat <- table(data[,var])
  new_dat <- data.frame(var = c(dat[1:length(dat)],ifelse(length(nas)==2, nas[2],0))) 
  rownames(new_dat) <- c(names(dat[1:length(dat)]), "NA")
  new_dat <- new_dat %>% rownames_to_column()
  dat_percent <- new_dat %>% mutate(percent = round(var/sum(var)*100,1)) 
  colnames(dat_percent) <- c(var, "n", "percent")
  return(dat_percent)
}
# ..................................................................................................................





#.....................................................................................................................
# this function calculates the influence of catergorical variables on Responder status #

# Arguments                                                                                         #
#  - data:  data                                                                                    #
#  - dependent.var: dependent variable in this case Responder                                       #
#  - idependent.var: independent variable tested for influencen on dependent variable               #
#  - dep.var.yes: name of the condition where the dependent variable is present, e,g. "pos" or "ja" #
#  - dep.var.no: name of the condition where the dependent variable is absent, e,g. "neg" or "nein" #
#  - independent.var.condition.yes: condition of the ind.var to be tested against other condition   #

two_by_two_Responder <- function(data,dependent.var,independent.var, dep.var.yes, dep.var.no, independent.var.condition.yes){
  Responder_ja <- sum(data[,dependent.var] == dep.var.yes &  data[,independent.var] == independent.var.condition.yes & !is.na(data[,independent.var]))
  Responder_nein <- sum(data[,dependent.var] == dep.var.yes &  data[,independent.var] != independent.var.condition.yes & !is.na(data[,independent.var]))
  Nonresponder_ja <- sum(data[,dependent.var] == dep.var.no &  data[,independent.var] == independent.var.condition.yes & !is.na(data[,independent.var]))
  Nonresponder_nein <- sum(data[,dependent.var] == dep.var.no &  data[,independent.var] != independent.var.condition.yes & !is.na(data[,independent.var]))
  two_by_two <- tibble(condition = c("+","-"),
                       Responder = c(Responder_ja, Responder_nein),
                       NonResponder = c(Nonresponder_ja, Nonresponder_nein))
  return(two_by_two)
}
# ..................................................................................................................











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
theme_Melanoma <- function(axis.text.size=10, Legend = TRUE,...){
  
  theme_custom <- theme_pubr()+
    theme(axis.title.x=element_text(face = "bold", size = axis.text.size+1),
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






#.................................................................................................................................
# these functions are used to ensure a uniform color palette in each plot #

scale_fill_Melanoma <- function(...){
  scale_fill_manual(values = c("darkgrey", "cornflowerblue","forestgreen","#d73027","purple", ...)) 
}

scale_color_Melanoma <- function(...){
  scale_color_manual(values = c("darkgrey", "cornflowerblue","forestgreen","#d73027","purple",...)) 
}
# ..................................................................................................................






#.................................................................................................................................
# this function customizes breaks and limits to have a uniform number of breaks in each plot #

# Arguments:                                                                        #
#  - y: inherited y object from ggplot function (does not need to be specified)     #
#  - n.breaks: changes number of breaks displayed                                   #

limit_fun <- function(y) {
  ymin <- floor(min(y))
  ymax <- ceiling(max(y))
  c(ymin,ymax)
}

break_fun <- function(y,n.breaks=n_breaks) {
  min_break <- 2*floor(min(y/2))
  max_break <- 2*ceiling(max(y/2))
  breaks <- seq(min_break,max_break, by= sum(max_break-min_break)/n.breaks)
  if (sum(max_break - min_break) <= 3) {
    unique(ceiling(breaks/0.5) *0.5)
  } else if (sum(max_break - min_break) > 3 & sum(max_break - min_break) <= 6){
    unique(round(seq(min_break,max_break, by= sum(max_break-min_break)/n.breaks),0))
  } else if (sum(max_break - min_break) > 6 & sum(max_break - min_break) <= 14){
    unique(2*round(seq(min_break,max_break, by= sum(max_break-min_break)/n.breaks)/2,0))
  } else {
    unique(4*round(seq(min_break,max_break, by= sum(max_break-min_break)/n.breaks)/4,0))
  }
}
# ..................................................................................................................






#.................................................................................................................................
# this function calculates t.test or wilcoxon test statistics with fdr adjustment #
# and adds x and y positions to plot custom p values in plots                     #

# Arguments:
#  - data: data
#  - x: independent variable
#  - y: dependent variable
#  - group.var: grouping variable if applicable
#  - p.adj: p adjustment method 
#  - signif: significance level (other results are filtered)
#  - method: can either be "t.test" or "wilcoxon"  


signif_Melanoma <- function(data,x,y, group.var=NULL, p.adj = "fdr",var.equal = FALSE, filter.p.adj= FALSE,
                            signif= 1, method="t.test"){
  
  if(filter.p.adj == FALSE){
    warning("Results filtered by raw p-values")
  }
  
  formula <- paste(
    y,"~",
    x, sep=" "
  )
  
  data <- data %>% filter(!is.na(get(x)))
  
  pvals <- compare_means(
    as.formula(formula),
    data=data,
    group.by = group.var,
    var.equal = var.equal,
    method=method,
    p.adjust.method = p.adj
  ) 
  
  dat.signif <- pvals %>% 
    add_significance("p.adj") %>%
    
    {
      if (filter.p.adj == TRUE) filter(.,p.adj < signif)
      else filter(.,p < signif)
    }
  
  return(dat.signif)
}
# ..................................................................................................................





# ..................................................................................................................
# this function plots desired variables against each other and adds significance level #

# Arguments:
#  - data: data
#  - diff.expr: Object calculated with previous function (signif_Melanoma)
#  - x: independent variable
#  - y: dependent variable
#  - facet: grouping variable if applicable, used for facetting
#  - p.label: specifies how the p-value is displayed, can be manually changed to p.adj, p.format, etc. 
#             mathematic operations like "round" are also possible within brackets
#  - add: specifies additional geometrics
#  - scale: specify if all plots shall have the same scale or if scales are free
#  - width: controls width of boxes
#  - outlier.shape: if jitter is added this helps avoid plotting points twice
#  - add.params: define shape, size, etc of added objects (like jitter)
#  - axis.text.size: size of axes
#  - ylab: y axis label
#  - significance: logical, if true p-value is displayed
#  - tip.length: defines length of tips (on bracket showing significance)
#  - p.size: size of p-value
#  - n_breaks: number of y axis breaks
#  - Legend: logical, if TRUE legend will be displayed

boxplot_Melanoma <- function(data,
                             diff.expr = NULL,
                             x,
                             y,
                             facet = NULL,
                             p.label = "p.adj.signif",
                             add = "jitter",
                             scales = "free",
                             width = 0.4,
                             outlier.shape=19,
                             add.params = list(color = "black",size = 1, width = 0.05),
                             axis.text.size = 8,
                             ylab = "log-2 miRNA Expression (a.u.)",
                             significance = TRUE,
                             tip.length=0,
                             p.size=5,
                             n_breaks = 10,
                             Legend = TRUE,
                             nrow.facet = 2,
                             ...){

  if(is.null(diff.expr))
    warning("No statistics are calculated as argument diff.expr is missing")

  significance.new <- if(is.null(diff.expr)) {
    FALSE
  } else{
    significance
  }

  data <- data %>%
    {if(!is.null(facet) & !is.null(diff.expr)) filter(.,get(facet) %in% unlist(diff.expr[,1]))
      else .}

  p <- data  %>%
    filter(!is.na(get(x))) %>%
    ggboxplot(x = x,
              y = y,
              fill = x,
              add = add,
              scales = scales,
              width = width,
              outlier.shape= ifelse(add == "jitter", NA, outlier.shape),
              add.params = add.params, ...
    )

  p.facet <- facet(p, facet.by = facet, nrow = nrow.facet, scales = scales) +
    ylab(ylab) +
    theme_Melanoma(axis.text.size = axis.text.size,Legend =Legend) +
    scale_color_Melanoma() +
    scale_fill_Melanoma()
  
  
  
  p.facet + if(significance.new == TRUE){stat_pvalue_manual(diff.expr , label = p.label,
                                                            tip.length = tip.length,size=p.size,
                                                            y.position = max(data[,y])*1.05)}
}
# ..................................................................................................................





# ..................................................................................................................
# this function is analogous to the previous one plotting a dotplot#

dotplot_Melanoma <- function(data,
                             diff.expr = NULL,
                             x,
                             y,
                             facet = NULL,
                             p.label = "p.adj.signif",
                             add = "mean_ci",
                             scales = "free",
                             size = 0.7,
                             width = 0.3,
                             outlier.shape=19,
                             add.params = list(size = 0.7, width = 0.05,color ="black"),
                             axis.text.size = 8,
                             ylab = "log-2 miRNA Expression (a.u.)",
                             significance = TRUE,
                             tip.length=0,
                             p.size=5,
                             position = position_jitter(0.02,0.02),
                             n_breaks = 10,
                             crossbar.width = 0.3,
                             crossbar.size = 0.25,
                             Legend = TRUE,
                             nrow.facet =2,
                             ...){

  if(is.null(diff.expr))
    warning("No statistics are calculated as argument diff.expr is missing")

  significance.new <- if(is.null(diff.expr)) {
    FALSE
  }
  else{
    significance
  }

  data <- data %>%
    {if(!is.null(facet)) filter(.,get(facet) %in% unlist(diff.expr[,1]))
      else .}

  p <- data  %>%
    filter(!is.na(get(x))) %>%
    ggdotplot(x = x,
              y = y,
              color = x,
              add = add,
              size=size,
              width = width,
              position = position,
              error.plot = "errorbar",
              add.params = add.params,
              ...)

  p.facet <- facet(p, facet.by = facet, nrow = nrow.facet, scales = scales) +
    stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                 geom = "crossbar", width = crossbar.width, color = "black", size=crossbar.size) +
    ylab(ylab) +
    theme_Melanoma(axis.text.size = axis.text.size,Legend =Legend) +
    scale_color_Melanoma()+
    scale_fill_Melanoma()

  p.facet +  if(significance.new == TRUE){stat_pvalue_manual(diff.expr , label = p.label,
                                                             tip.length = tip.length,size=p.size,
                                                             y.position = max(data[,y]))}
}
# ..................................................................................................................




# ..................................................................................................................
# this function plots desired variables against each other and adds significance level #

# Arguments:
#  - data: data
#  - x: independent variable
#  - y: dependent variable
#  - facet: grouping variable if applicable, used for facetting
#  - signif: significance level (other results are filtered)
#  - method: method accepted by compare_means
#  - p.label: specifies how the p-value is displayed, can be manually changed to p.adj, p.format, etc. 
#             mathematic operations like "round" are also possible within brackets
#  - plot.type: specifies desired plot type. can be chosen between "dotplot" and "boxplot"
#  - data: data
#  - group.var: grouping variable if applicable
#  - p.adj: p adjustment method 
#  - signif: significance level (other results are filtered)


signif_plot_Melanoma <- function(data, x, y, signif=1, method="t.test",p.adj = "fdr", facet =NULL, var.equal = FALSE,
                                 nrow.facet = 2, filter.p.adj= FALSE, p.label="p = {round(p,4)}",
                                 plot.type = "boxplot",  ...){
  if (method == "t.test") {
    message("Using T-Test for statistical comparison")
  } else {
    message("Using Wilcoxon Rank Sum Test for statistical comparison")
  }
  data <- data %>% filter(!is.na(get(x)))
  significance <- signif_Melanoma(data=data, x=x, y=y, signif=signif,p.adj = p.adj, var.equal = var.equal,
                                  group.var = facet, method = method, filter.p.adj = filter.p.adj)
  if (!is.null(facet)) {
    names(significance)[1] <- facet
  }

  if (plot.type == "boxplot"){
    plot <- boxplot_Melanoma(data=data,diff.expr = significance, x=x, y=y,facet=facet, p.label = p.label,
                             nrow.facet = nrow.facet,...)
  } else if (plot.type == "dotplot"){
    plot <- dotplot_Melanoma(data=data,diff.expr = significance, x=x, y=y,facet=facet, p.label = p.label,
                             nrow.facet = nrow.facet,...)
  } else {
    stop("please specify plot type: either \"dotplot\" or \"boxplot\" ")
  }

  list(
    stat_test_results = significance,
    graph = plot
  )
}
# ......................................................................................................................









######################################################
#                                                    #
#               4. Machine Learning                  #
#                                                    #
######################################################

# ......................................................................................................................
# this function calculates cross validated AUC with 95% confidence intervals for a given caret model object #

ci.cv_ROC <- function(data){
  obs <-data$pred$obs
  pred <- data$pred$pos
  
  obs <-split(obs , f = data$pred$Resample)
  pred <-split(pred , f = data$pred$Resample)
  
  list(
    cv_AUC = cvAUC(pred,obs),
    ci.cvAUC = ci.cvAUC(pred,obs))
}
# ......................................................................................................................






#....................................................................................................................
# this function calculates the weighted importance of a model ensemble              #
# it calculates importance for single features in each model and then combines the  #
# results in a weighted approach by model importance in the stacked model           # 

# Arguments:
#  - models: vector with model names
#  - dat: data
weightedImp <- function(models, dat){
  
  #calculate importance
  imp <- sapply(models,data=dat, function(names,data){
    varImp(data$models[[names]])$importance
  }) %>% 
    #transform to data frame
    as.data.frame() %>% 
    # remove duplicate entries
    select(-contains("neg")) %>%
    setNames(models)
  
  
  #calculate weights based on the model importance of the individual models within the ensemble
  weights <- varImp(ensemble_3$ens_model)$importance %>%
    select(neg) %>% 
    mutate(neg = neg/100) %>%
    unlist() %>%
    setNames(models)
  
  #combine weights with the separate model importances
  weighted_imp <- t(t(imp)*weights) %>% 
    data.frame()%>%
    rownames_to_column("Feature") %>%
    mutate(overall = rowSums(.[,-1])) %>% 
    mutate(percent = overall/sum(overall)*100) %>%
    mutate(scaled = percent/max(percent)*100) %>%
    select(c(Feature, scaled)) %>%
    arrange((scaled)) %>% 
    mutate(Feature = str_replace_all(.$Feature,"X","")) %>% 
    mutate(Feature = str_replace_all(.$Feature,"hsa","")) %>% 
    mutate(Feature = str_replace_all(.$Feature,"\\.\\.","")) %>%
    mutate(Feature = str_replace_all(.$Feature,"\\.$","")) %>%
    mutate(Feature = str_replace_all(.$Feature,"\\.","-")) 
  
  return(weighted_imp)
}






## plot for exosome data
ggexosome <- function(data,
                      x,
                      y,
                      group = "response",
                      size = 1.7,
                      se = F,
                      method = "glm",
                      fullrange = T
){
  dat <- data
  dat %>%
    ggplot(
      aes_string(x=x,y=y,shape = group, color = group, lty = group)
    ) +
    geom_point(size=size) +
    geom_smooth(method = method, se = se, fullrange=fullrange) +
    theme_Melanoma() +
    theme(legend.key.size = unit(3,"line"),
          legend.title = element_blank())+
    scale_color_Melanoma()
    
}














#....................................................................................................................
# this function calculates accuracy, F1, spec, sens,... for the Folds    #
# within cross validation using the resamples with the optimal tuning    #
# parameters

cv_res <- function(modelNames){
  
  # show ROC of models
  model_results <- lapply(names(modelNames), function(x){
    modelNames[[x]]$results %>%
      filter(ROC == max(ROC)) %>%
      select(c(ROC, Sens, Spec)) %>% .[1,]
  }) 
  
  model_results <- do.call(rbind.data.frame, model_results) 
  
  
  ## extract those resample with the optimal tuning parameters from the model list for each model
  model_best_tune <- lapply(names(modelNames), function(x){
    tmp <- modelNames[[x]]$pred
    
    for (i in 1:length(modelNames[[x]]$bestTune)){
      tmp   <- tmp %>%
        filter(!!sym(names(modelNames[[x]]$bestTune)[i]) == unlist(modelNames[[x]]$bestTune[i]))
    }
    
    return(tmp)
  })
  
  
  # calculate accuracy, sensitivity and so on for the cross validated models 
  cv_results  <- sapply(1:length(modelNames), function(y){
    
    resamples <- modelNames[[1]]$pred$Resample %>% unique()
    
    sapply(resamples, function(x){
      tmp_res <- model_best_tune[[y]] %>%
        filter(Resample == x)
      
      conf_mat <- confusionMatrix(tmp_res$pred, tmp_res$obs)
      conf_mat$byClass
    }) %>%
      t() %>%
      data.frame() %>%
      summarize_all(mean, na.rm =T)
  })%>%
    data.frame() %>%
    setNames(names(modelNames)) %>%
    t()
  
  cbind(ROC = model_results$ROC,cv_results)
  
}












#....................................................................................................................
# this function averages results from cv_res function for different Seeds    #
modelSeed_meanSD <- function(y){
  model_metrics <- lapply(1:length(y), function(x){
    cv_res(y[[x]])
  })
  
  metricList <- lapply(names(y[[1]]),function(modelName){
    tmp <- do.call(rbind, model_metrics) %>%
      data.frame() %>%
      rownames_to_column("model")%>%
      filter(str_detect(model, modelName)) %>%
      select(-model) %>%
      mutate_all(as.numeric) 
    
    tmp_mean <- tmp %>%
      summarize_all(mean)
    tmp_sd <- tmp %>%
      summarize_all(sd)
    
    
    data.frame(mean = t(tmp_mean), sd = t(tmp_sd)) %>%
      rownames_to_column("metric") %>%
      mutate(error = qt(0.975, length(y) -1)*sd/sqrt(length(y)),
             lower = mean - error,
             upper = mean + error) %>%
      column_to_rownames("metric")
  })
  
  names(metricList) <- names(y[[1]])
  return(metricList)
}





#....................................................................................................................
# calculate model metrics 
calcTestMetrics <- function(data, models){
  metricList_test <- lapply(models,function(x){
    tmp <- do.call(rbind, data) %>%
      data.frame() %>%
      rownames_to_column("model")%>%
      filter(str_detect(model, x)) %>%
      select(-model) %>%
      mutate_all(as.numeric) 
    
    tmp_mean <- tmp %>%
      summarize_all(mean)
    tmp_sd <- tmp %>%
      summarize_all(sd)
    
    n <- length(names(model_list[[1]]))
    
    data.frame(mean = t(tmp_mean), sd = t(tmp_sd)) %>%
      rownames_to_column("metric") %>%
      mutate(error = qt(0.975,n-1)*sd/sqrt(n),
             lower = mean - error,
             upper = mean + error) %>%
      column_to_rownames("metric")
  })
  
  names(metricList_test) <- names(model_list[[1]])
  return(metricList_test)
}





#....................................................................................................................
# function to combine predicted values with real values from a list
# x: name of list element
# model: name of model
predictList <- function(x,model){
  data.frame(Sample = x,
             Responder = predict.train(model_list[[x]][[model]], newdata = dat_split[[x]]$testing,type="prob"),
             RealClass = dat_split[[x]]$testing$Responder) 
}



#....................................................................................................................
# function to calculate TPR and FPR for different cutoffs averaged over a distinct number of samples defined by "supgroup"
# x: predicted values
# class: class containing real values
mean_roc <- function(data, cutoffs = seq(from = 0, to = 1, by = 0.1)) {
  map_df(cutoffs, function(cp) {
    out <- cutpointr(data = data, x = PredictionValues, class = RealClass,
                     subgroup = Sample, method = oc_manual, cutpoint = cp,
                     pos_class = "neg", direction = ">=")
    data.frame(cutoff = cp, 
               sensitivity = mean(out$sensitivity),
               specificity = mean(out$specificity))
  })
}






#....................................................................................................................
# function to calculate AUC using predicted and observed values
# probs: predicted values (probabilities)
# true_Y: observed values
getROC_AUC <- function(probs, true_Y){
  probsSort <- sort(probs, decreasing = TRUE, index.return = TRUE)
  val <- unlist(probsSort$x)
  idx <- unlist(probsSort$ix)  
  
  roc_y <- true_Y[idx];
  stack_x <- cumsum(roc_y == "pos")/sum(roc_y == "pos")
  stack_y <- cumsum(roc_y == "neg")/sum(roc_y == "neg")    
  
  auc <- sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
  return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
}








#....................................................................................................................
# function to calculate the average AUC of multiple ROC models
# samples: number of different models
mean_AUC <- function(data,samples){
  tmp <- sapply(samples, function(x){
    dat <- data %>% 
      filter(Sample == x)
  getROC_AUC(dat$PredictionValues, dat$RealClass) %>%
    .$auc
  }) 
  data.frame(mean.auc = mean(tmp),
             error = confInt(tmp)) %>%
    mutate(lower = mean.auc- error, 
           upper = mean.auc+error)
  
}









#....................................................................................................................
# function to combine predictions from different models (generated by different seeds) within a list
# out object gives the predictions table, ROC and AUC 
pred_summary <-function(model, list.element = 1:length(n)){
  pred <- map_df(list.element,model = model, predictList) %>% 
    select(-Responder.pos) %>%
    setNames(c("Sample","PredictionValues", "RealClass"))
  mr <- mean_roc(pred)
  ma <- mean_AUC(pred, list.element)
  list(predictions = pred,
       ROC = mr,
       AUC = ma
  )
}



#....................................................................................................................
# load saved models
loadModels <- function(model){
  if(model == "M1"){
    readRDS(file = "model_base_seeds.Rds")
  } else if(model == "M2") {
    readRDS(file = "model_miR_seeds.Rds")
  } else if(model == "M3") {
    readRDS(file = "model_best_perf_seeds.Rds")
  } else {
    stop("Please specify model (one of M1 (base model), M2 (miR model), M3 (best performing model))")
  }
}


  


### 
ggexosome2 <- function(data,x,y,facet.by=NULL,col,  method = "t.test"){
  data %>%  filter(!is.na(!!sym(x))) %>%
    ggboxplot(x,y, facet.by = facet.by, 
              ylab = y, 
              add = "jitter",
              add.params = list(color = col),
              outlier.shape = NA, 
              ) +
    stat_compare_means(method = method) +
    ylim(min(data[y])/1.1, max(data[y])*1.1)
}

# arrange a grid of plots for exosome data
exo_arrange <- function(data,x,y,facet.by,col="black", ylab){
  p1 <- ggexosome2(data, x,y[[1]],facet.by,col) + ylab(ylab[[1]])
  p2 <- ggexosome2(data, x,y[[2]],facet.by,col) + ylab(ylab[[2]])
  p3 <- ggexosome2(data,x,y[[3]], facet.by,col) + ylab(ylab[[3]])
  ggarrange(p1,p2,p3,nrow=1)
}