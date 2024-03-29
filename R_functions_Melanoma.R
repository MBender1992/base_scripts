#<<<<<<<<<<<HEAD

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

load_melanoma_data <- function(){
  library(readxl)

  # load csv files
  dat_miR   <- read_csv("Data/miRNA_Expression_Fireplex_Melanoma_Study.csv")
  dat_meta  <- read_xlsx("Data/Metadata_Melanoma_Study.xlsx") %>%
    select(-c(therapy_start, Abnahmedatum)) %>%
    mutate(TRIM_PDL1_Expression = str_replace_all(TRIM_PDL1_Expression,"\\++","+")) %>%
    mutate(TRIM_PDL1_Expression = ifelse(TRIM_PDL1_Expression == "o", NA,TRIM_PDL1_Expression)) %>%
    mutate(Stadium = toupper(Stadium)) %>%
    mutate(Stadium = str_extract(Stadium, "^[IV]{1,3}")) %>%
    mutate(BRAF = str_replace_all(BRAF, "\\.", "")) %>%
    mutate(breslow_thickness_mm = parse_number(breslow_thickness_mm))

  # change ID column to uniform capital letters for later filtering
  names(dat_miR) <- c("miRNA", toupper(names(dat_miR)[-1]))

  # define IDs to be dropped for further analyses
  controls <- c("K104_1", "K104_2", "K104_3A", "K104_3B")
  duplicates <- c("22B","38B","39B","47B")

  # wide miR data (78 patients with miRNA data)
  dat_miR_trans <- transpose_dataframe(colnames = c("ID",dat_miR$miRNA), data = dat_miR) %>%
      filter(!ID %in% controls & !ID %in% duplicates) %>%   #drop duplicate patient data
      mutate(ID = parse_number(ID)) #convert ID to numeric

  # join both tables
  right_join(dat_miR_trans,dat_meta, by="ID") %>%
    filter(!ID %in% c(1,2)) %>% # no data available for patient 1 and 2 but still part of the source table
    mutate(miRExpAssess = ifelse(is.na(rowSums(.[,which(str_detect(names(.),"mir"))])), 0,1))  %>%# if no miRNA expression has been measure fill in 0
    arrange(ID) %>%
    mutate(Responder = factor(Responder, levels = c("nein", "ja"), labels = c("no", "yes"))) %>%
    mutate(prior_BRAF_therapy = ifelse(str_detect(Vorbehandlung,"Mek|Dabra|Tafinlar|Tefinlar|MEK|BRAF|Vemu|[zZ]ellboraf"), 1, 0)) %>%
    select(-Vorbehandlung)
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
  require(ggpubr)
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

  data <- data %>% filter(!is.na(get(x)) & !is.na(get(y)))

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
                             add.params = list(size = 0.7, width = 0.1,color ="black"),
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

# ......................................................................................................................
# ci_cv_AUC <- function(model){
#   dat <- filter(model$pred, lambda == model$finalModel$lambdaOpt)
#
#   obs <- dat$obs
#   event_name <- levels(obs)[2]
#   pred <- dat[[event_name]]
#
#   obs <- split(obs, f = dat$Resample)
#   pred <- splot(pred, f = dat$Resample)
#
#   ci.cvAUC(pred, obs)
# }





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



# function to display p-values in table1
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- dat_table1[[name]]
    ind <- !is.na(y)
    y <- y[ind]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ dat_table1$Responder[ind])$p.value
    } else {
      p <- chisq.test(table(y, droplevels(dat_table1$Responder[ind])))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}


#....................................................................................................................
## function to calculate AUC, sens and spec for training and test set and extract model coefficients from a model if applicable
## x = model matrix with predictors
## y = response variable as factor
## method = method used for cross validation
## number = number of folds if k-fold cv used
## reapeats = number of repeats if repeatedcv used
## train.method = model method used in caret::train function
## metric = metric used to assess model quality
## tuneGrid = grid with hyperparameters to be tuned



# calc.model.metrics.2 <- function(x.train, y.train, x.test, y.test, train.method = "glmnet",
#                                  cv.method = "repeatedcv", number = 10, repeats = 5, metric = "ROC", tuneGrid){
#
#   # define ctrl function
#   cctrl1 <- trainControl(method=cv.method, number=number,repeats = repeats, returnResamp="all",savePredictions = T,
#                          classProbs=TRUE, summaryFunction=twoClassSummary)
#
#   # run glmnet model
#   md <- train(x.train, y.train, method = train.method,preProcess = c("center","scale"),
#               trControl = cctrl1,metric = metric,tuneGrid = tuneGrid)
#
#   # obtain cv AUC of training folds
#   ci_cv <- ci.cv.AUC(md)
#
#   # train coefs
#   feat <- coef(md$finalModel, md$finalModel$lambdaOpt)
#
#   # obtain index from max metric
#   opt <- md$results[which(md$results$lambda == md$finalModel$lambdaOpt),]
#
#   # predict
#   pred <- predict(md, x.test, type="prob")
#
#   # object to return
#   res <- list(
#     predictions = data.frame(pred.ja = pred$yes, obs = y.test),
#     coefficients = rownames_to_column(data.frame(vals = feat[feat[,1] != 0, 1][-1]),"coefs"),
#     train.metrics = opt[which(opt$ROC == max(opt$ROC)),],
#     train.cv = data.frame(cvAUC = ci_cv$cvAUC,
#                           se = ci_cv$se,
#                           lower = ci_cv$ci[1],
#                           upper = ci_cv$ci[2]),
#     test.metrics = data.frame(AUC = auc(roc(y.test, pred[,1], direction = ">", levels = c("no", "yes"))),
#                               Sens = sensitivity(y.test, predict(md, x.test, type="raw"))  ,
#                               Spec = specificity(y.test, predict(md, x.test, type="raw")))
#   )
#   return(res)
# }





# function to select which columns to use for the model matrix as specified by the model argument
model.matrix.subset <- function(model, data){
  if(model == "complete"){
    mm <- model.matrix(Responder~., data = data)[,-1]
  } else if(model == "miRNA"){
    mm <- model.matrix(Responder~., data = select(data, c(contains("mir"),Responder)))[,-1]
  } else if(model == "baseline"){
    mm <- model.matrix(Responder~., data = select(data, contains(c("Eosinophile","LDH","S100","CRP")),Responder))[,-1]
  } else if(model == "signif"){
    mm <- model.matrix(Responder~., data = select(data, contains(readRDS("significant_features.rds")),Alter,prior_BRAF_therapy,Responder))[,-1]
  } else if(model == "relaxedLasso"){
    mm <- model.matrix(Responder~., data = select(data, c(feat.relaxed$coef, Alter, prior_BRAF_therapy,Responder)))[,-1]
  } else if(model == "relaxedLassomiRNA"){
    mm <- model.matrix(Responder~., data = select(data, c(feat.relaxed.miRNA$coef,Responder)))[,-1]
  } else {
    stop("Please specify 1 of the following 4 options:
    1. 'baseline' for a base model using conventional serum markers (LDH, CRP, S100, Eosinophile)
    2. 'miRNA' for a model using only miRNAs (reduced by lasso to informative features)
    3. 'signif' for a model with significantly different features between responders and non-responders
    4. 'complete' for a model with all predictors (reduced by lasso)
    5. 'relaxedLasso' for a model with the best predictors selected by the 'complete' model (afterwards reduced again with LASSO)
    6. 'relaxedLassomiRNA for a model with the best predictors selected by the 'miRNA' model (afterwards reduced again with LASSO")
  }
  return(mm)
}


# models a function based on a presepcified model and evaluates training and test test using ROC, Sensitivity and Specificity
# mlEval <- function(modelMatrix, dat, rep, k, train.method = "glmnet", tuneGrid = expand.grid(alpha = 1, lambda = seq(0.01,0.2,by = 0.01))){
#   # define model matrix with selected features
#   x <- modelMatrix
#
#   # activate parallel computing
#   cl <- makeCluster(detectCores(), type='PSOCK')
#   registerDoParallel(cl)
#
#   # generate 10 folds for outer loop
#
#   set.seed(12)
#   fold.train <- createMultiFolds(y, k = k, times = rep) # ensure that at least 10 samples are in each fold
#
#   # split data based on these folds (Fold1 means that Fold1 is used for testing)
#   train.test.folds <- lapply(c(1:rep), function(split){
#
#     # select only folds containing the specified repeat in each iteration
#     if(split == 10){
#       ind <- names(fold.train) %>% str_detect("Rep10")
#       dat <- fold.train[ind]
#     } else {
#       ind <- names(fold.train) %>% str_detect(paste("Rep0",split, sep =""))
#       dat <- fold.train[ind]
#     }
#
#     # split data into training and test set with each fold being the test set once
#     res <- lapply(c(1:k), function(fold){
#       list(x.test = x[-dat[[fold]],],
#            x.train = x[dat[[fold]],],
#            y.test = y[-dat[[fold]]],
#            y.train = y[dat[[fold]]]
#       )
#     })
#     return(res)
#   })
#
#   # define name of the list elements
#   reps <- paste0("Rep", 1:rep)
#   folds <- paste0("Fold", 1:k)
#   train.test.folds <- setNames(lapply(train.test.folds, setNames, folds), reps)
#
#   set.seed(849)
#   lapply(c(1:rep), function(split){
#     # select Data from 1 repeat
#     dat <- train.test.folds[[paste("Rep",split, sep ="")]]
#     # print message to follow progress
#     message(paste("Starting calculation of Rep", split,"... of", rep))
#     # apply model to all folds of that 1 repeat and test against the remaining fold not used for training
#     res <- pblapply(c(1:k), function(fold){
#       calc.model.metrics.2(x.train = dat[[fold]]$x.train, y.train = dat[[fold]]$y.train, x.test =dat[[fold]]$x.test,
#                            y.test = dat[[fold]]$y.test, train.method = train.method,
#                            tuneGrid = tuneGrid)
#     })
#   })
#   stopCluster(cl)
# }



# convert metrics from training data list to dataframe
# unlist_model <- function(dat,metric, element){
#   ls <- lapply(1:rep, function(split){
#     do.call(rbind.data.frame, sapply(dat[[split]], '[', element))
#   })
#   unlist(sapply(ls, "[", metric))
# }


# extract coefficients from model list
# extract_coefs <- function(data){
#   lapply(1:10, function(x){
#     tmp <- sapply(sapply(data[[x]], '[', 'coefficients'), '[', 'coefs') %>% unlist()
#     data.frame(coef = tmp)
#   })
# }




#
# rbind_model_ci <- function(model){
#
#   ls <- ls_cvAUC(model)
#   out <- cvAUC(ls$predictions, ls$labels)
#   res <-  ci.cvAUC(ls$predictions, ls$labels)
#
#   rbind(cv.AUC.inner = data.frame(cvAUC = mean(unlist_model(model, "cvAUC", "train.cv")),
#                                  lower = mean(unlist_model(model, "lower", "train.cv")),
#                                  upper = mean(unlist_model(model, "upper", "train.cv"))),
#         cv.AUC.outer = data.frame(cvAUC = res$cvAUC,
#                                 lower = min(res$ci),
#                                 upper = max(res$ci)
#                                 )
#   )
# }




## function to split data in a way that the function cvAUC can use it for calculation of ROC, se and confidence intervals
# ls_cvAUC <- function(data){
#   predictions <- lapply(1:10, function(d){
#     tmp <- sapply(data[[d]], "[", "predictions")
#     lapply(1:10, function(x){
#       1-tmp[[x]]$pred.ja
#     })
#   })
#
#
#   labels <- lapply(1:10, function(d){
#     labs <- sapply(data[[d]], "[", "predictions")
#     lapply(1:10, function(x){
#       labs[[x]]$obs
#     })
#   })
#
#   list(predictions = lapply(rapply(predictions, enquote, how="unlist"), eval),
#        labels = lapply(rapply(labels, enquote, how="unlist"), eval))
# }






##
stat_test_BRAF <- function(data, var,  p.adj.anova = "fdr", hccm = TRUE){

  formula <- as.formula(paste("logexp~",var))

  # test for significant differences
  aov_test <-data %>%
    group_by(miRNA) %>%
    anova_test(formula, white.adjust = hccm) %>%
    adjust_pvalue(method = p.adj.anova) %>%
    filter(p.adj < 0.05)

  # calculate GH post-hoc test
  GH_test <- data %>%
    group_by(miRNA) %>%
    games_howell_test(formula)%>%
    filter(p.adj < 0.05) %>%
    filter(miRNA %in% aov_test$miRNA)

  # calculate maximum value for each miRNA
  maxexp <- data %>%
    filter(miRNA %in% GH_test$miRNA) %>%
    group_by(miRNA) %>%
    summarize(y.position = max(expression)*1.06)

  # calculate y.positions
  ypos <- left_join(data.frame(miRNA =GH_test$miRNA),  data.frame(maxexp)) %>%
    mutate(y.position = ifelse(duplicated(y.position), y.position*1.11, y.position))

  #
  stat_test <- GH_test %>%
    add_xy_position(x = var) %>%
    mutate(y.position = ypos$y.position)

  return(stat_test)
}
