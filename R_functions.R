setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/base_scripts")
source("R_functions_PhD.R")
source("R_functions_Melanoma.R")

#4....................................................................................................................
### custom.Heatmap
custom.Heatmap <- function(data,cols=c("#003399", "white", "#990000"),legend.title = "row Z-score",...) {
  data.scale <- scale(t(data))
  data.scale <- t(data.scale)
  colors <- colorRampPalette(cols, space="rgb")(128)
  Heatmap(data.scale, 
          cluster_rows = T,
          cluster_columns = T,
          clustering_method_columns = "average",
          clustering_method_row = "average",
          clustering_distance_row = "pearson",
          clustering_distance_column = "pearson",
          col = colors,
          column_dend_height = unit(1, "cm"), 
          row_dend_width = unit(1, "cm"),
          row_names_side = "right",
          row_names_gp = gpar(fontsize = 6),
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = legend.title,
                                      title_gp = gpar(fontsize = 13),
                                      color_bar="continuous",
                                      title_position ="topcenter",
                                      legend_direction = "vertical",
                                      labels_gp = gpar(fontsize = 12),
                                      legend_width = unit(4, "cm"),
                                      at = c(-2,-1,0,1,2)),
          rect_gp = gpar(col = "white", lty = 1, lwd = 1),...)
}






#5....................................................................................................................
###nicer theme 
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library("ggthemes")
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing = unit(0.2, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
}



#6....................................................................................................................
###multiple replacement
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}




#7....................................................................................................................
###summary statistics (mean, sd) can be widened by median, sum, max, min, CI
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}







#8....................................................................................................................
###barplot Publication Cancers 2019
custom.barplot <- function(a,b,c,d,position=position_dodge(0.7),size=1,lty=2,yintercept=0.58,...) {
  p <- ggbarplot(a, x=b,y=c,fill=d,position=position,size=size,...) + theme_Publication() +
    theme(axis.text.y=element_text(size=12),
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=40, size=9,vjust=0.7),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(size=1),
          axis.ticks.length=unit(0.3,"cm"),
          axis.title.y=element_text(size= 12,face="plain"), 
          plot.title = element_text(size=16),
          legend.text=element_text(),
          legend.position = "none",
          legend.key.size = unit(0.9,"line")) +
    ylab(expression(paste("Log"[2], "FC"))) +
    guides(fill=guide_legend(title=d)) +ggtitle(gsub("_","-",c)) +
    scale_fill_manual(values=c("#9E0E0E", "#0E3E9E")) 
    geom_hline(yintercept=yintercept,lty=lty)
  
  return(p)
}





#9....................................................................................................................
# define n and print x labels using this n
nlabel_x <- function(data_column,a){
  ind_a <- which(data_column == a)
  len <- length(ind_a)
  paste(a," (n = ", len,")", sep="")
}




# 10....................................................................................................................
custom_hline <- function(data,a,b,...){
  geom_hline(data=filter(data, cell_line==a & cell_cycle ==b), aes(yintercept=mean[treatment =="KO"]), lty=2,...)
}






# 11 permutation p-value....................................................................................................................
#permutation functions
permutation_pval <- function(data,miR, cell, B=10000){
  con <- data %>% filter(irradiation == "control" & miRNA == miR & cell_line == cell) %>%
    select(expression) %>% unlist
  irr <- data %>% filter(irradiation == "KAUVIR" & miRNA == miR & cell_line == cell) %>%
    select(expression) %>% unlist
  
  N <- ifelse(length(con) < length(irr), length(con), length(irr))
  obs <- mean(con) - mean(irr)
  dat <- c(con,irr)
  
  set.seed(1)
  obsstar <- replicate(B,{
    shuffle <- sample( dat )
    irr_star <- shuffle[1:N]
    con_star <- shuffle[(N+1):(2*N)]
    mean(irr_star)-mean(con_star)
  })
  
  #pvalue
  (sum(abs(obsstar) >= abs(obs)) +1) / (length(obsstar)+1) 
}
permutation_pval_total <- function(cell, data){
  sapply(mirs, permutation_pval, cell= cell, data=data)
}



#12 ....................................................................................................................
####Threshold Function to determine best cutoff ROCR package
threshold1 <- function(predict, response) {
  perf <- ROCR::performance(ROCR::prediction(predict, response), "sens", "spec")
  df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@x.values[[1]], spec = perf@y.values[[1]])
  df[which.max(df$sens + df$spec), "cut"]
}

#13 custom theme....................................................................................................................
theme_custom <- function (base_size = 12, base_family = "Roboto Condensed") {
  half_line <- base_size/2
  theme(line = element_line(color = "black", size = 0.5, linetype = 1, lineend = "butt"),
        rect = element_rect(fill = "white", color = "black", size = 0.5, linetype = 1),
        text = element_text(family = base_family, face = "plain", color = "black",
                            size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5,
                            angle = 0, margin = margin(), debug = F),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
        axis.text = element_text(size = base_size * 1.1, color = "gray30"),
        axis.text.x = element_text(margin = margin(t = 0.8 * half_line/2), vjust = 1),
        axis.text.x.top = element_text(margin = margin(b = 0.8 * half_line/2), vjust = 0),
        axis.text.y = element_text(margin = margin(r = 0.8 * half_line/2), hjust = 1),
        axis.text.y.right = element_text(margin = margin(l = 0.8 * half_line/2), hjust = 0),
        axis.ticks = element_line(color = "gray30", size = 0.7),
        axis.ticks.length = unit(half_line / 1.5, "pt"),
        axis.title.x = element_text(margin = margin(t = half_line), vjust = 1,
                                    size = base_size * 1.3, face = "bold"),
        axis.title.x.top = element_text(margin = margin(b = half_line), vjust = 0),
        axis.title.y = element_text(angle = 90, margin = margin(r = half_line),
                                    vjust = 1, size = base_size * 1.3, face = "bold"),
        axis.title.y.right = element_text(angle = -90, vjust = 0,
                                          margin = margin(l = half_line)),
        legend.background = element_rect(color = NA),
        legend.spacing = unit(0.4, "cm"),
        legend.spacing.x = NULL,
        legend.spacing.y = NULL,
        legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
        legend.key = element_rect(fill = "gray95", color = "white"),
        legend.key.size = unit(1.2, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(size = rel(0.8)),
        legend.text.align = NULL,
        legend.title = element_text(hjust = 0),
        legend.title.align = NULL,
        legend.position = "right",
        legend.direction = NULL,
        legend.justification = "center",
        legend.box = NULL,
        legend.box.margin = margin(0, 0, 0, 0, "cm"),
        legend.box.background = element_blank(),
        legend.box.spacing = unit(0.4, "cm"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "gray30",
                                    fill = NA, size = 0.7),
        panel.grid.major = element_line(color = "gray90", size = 1),
        panel.grid.minor = element_line(color = "gray90", size = 0.5,
                                        linetype = "dashed"),
        panel.spacing = unit(base_size, "pt"),
        panel.spacing.x = NULL,
        panel.spacing.y = NULL,
        panel.ontop = F,
        strip.background = element_rect(fill = "white", color = "gray30"),
        strip.text = element_text(color = "black", size = base_size),
        strip.text.x = element_text(margin = margin(t = half_line,
                                                    b = half_line)),
        strip.text.y = element_text(angle = -90, margin = margin(l = half_line,
                                                                 r = half_line)),
        strip.placement = "inside",
        strip.placement.x = NULL,
        strip.placement.y = NULL,
        strip.switch.pad.grid = unit(0.1, "cm"),
        strip.switch.pad.wrap = unit(0.1, "cm"),
        plot.background = element_rect(color = NA),
        plot.title = element_text(size = base_size * 1.8, hjust = 0.5,
                                  vjust = 1, face = "bold",
                                  margin = margin(b = half_line * 1.2)),
        plot.subtitle = element_text(size = base_size * 1.3, hjust = 0.5, vjust = 1,
                                     margin = margin(b = half_line * 0.9)),
        plot.caption = element_text(size = rel(0.9), hjust = 1, vjust = 1,
                                    margin = margin(t = half_line * 0.9)),
        plot.tag = element_text(size = rel(1.2), hjust = 0.5, vjust = 0.5),
        plot.tag.position = "topleft",
        plot.margin = margin(base_size, base_size, base_size, base_size), complete = T)
}






# 14 Heatmap Colors....................................................................................................................

HeatmapColors <- function(cols = c("#003399", "white", "#990000"), min.z = -2, max.z=2, low.w=-0.7, up.w=0.7,n=129){
colors <- colorRampPalette(cols, space="rgb")(n);

# Contrast parameters
ncols <- 18                          # Number of colors (out of 128) used  
# for both ends.
# end1 = lower whisker to min value
# end2 = upper whisker to max value
min.z <- min.z;#min(as.vector(mat));        # The min intensity in the matrix
max.z <- max.z;#max(as.vector(mat));        # The max intensity in the matrix
#b <- boxplot(as.vector(mat),plot=F); # Get the distribution
low.w <- low.w;#b$stats[1,1];               # Lower whisker
up.w <- up.w;#b$stats[5,1];                # Upper whisker

# Now create the non-linear color bar values
# breaks is a numerical vector of length 128
# each value corresponding to a color
# breaks will be given to 'image' function
breaks <- c(
  seq(min.z,low.w,by=(low.w-min.z)/ncols),
  seq(low.w,up.w,by=(up.w-low.w)/(128-2*ncols)),
  seq(up.w,max.z,by=(max.z-up.w)/ncols)
);
length(breaks)

breaks <- unique(round(breaks,4));
# image(t(mat[nrow(mat):1,]),col=colors,breaks=breaks,axes=FALSE);

# colorbar (just to plot the color bar indexed by intensities)
par(omi=c(0,0.1,0.1,0.1),mar=c(4,0,0,0))
image(matrix(1:length(colors), ncol=1), col=colors, xaxt="n", yaxt="n");
axis(1, at = c(1,ncols,64,128-ncols,128)/128, labels =
       round(c(min.z,low.w,0,up.w,max.z),2),las=2)  

library("circlize")
colors2 <- colorRamp2(breaks=breaks, colors=colors, transparency = 0, space = "RGB")
}











# # 15 generate custom patient table....................................................................................................................
# patient_table <- function(data, var){
#   nas <- table(is.na(data[,var]))
#   dat <- table(data[,var])
#   new_dat <- data.frame(var = c(dat[1:length(dat)],ifelse(length(nas)==2, nas[2],0))) 
#   rownames(new_dat) <- c(names(dat[1:length(dat)]), "NA")
#   new_dat <- new_dat %>% rownames_to_column()
#   dat_percent <- new_dat %>% mutate(percent = round(var/sum(var)*100,1)) 
#   colnames(dat_percent) <- c(var, "n", "percent")
#   return(dat_percent)
# }
# 
# 
# 
# 
# 
# 
# 
# # 16 ....................................................................................................................
# ### Funktion, um chi square test für verschiedene Variablen in Bezug auf den Responder Status zu untersuchen
# two_by_two_Responder <- function(dat_var1,var1_cond1,var1_cond2,dat_var2, var2_cond){
#   Responder_ja <- sum(dat_var1 == var1_cond1 &  dat_var2 == var2_cond & !is.na(dat_var2))
#   Responder_nein <- sum(dat_var1 == var1_cond1 &  dat_var2 != var2_cond & !is.na(dat_var2))
#   Nonresponder_ja <- sum(dat_var1 == var1_cond2 &  dat_var2 == var2_cond& !is.na(dat_var2))
#   Nonresponder_nein <- sum(dat_var1 == var1_cond2 &  dat_var2 != var2_cond & !is.na(dat_var2))
#   two_by_two <- tibble(condition = c("+","-"),
#                        Responder = c(Responder_ja, Responder_nein),
#                        NonResponder = c(Nonresponder_ja, Nonresponder_nein))
#   return(two_by_two)
# }






# 17 cross validated ROC analysis....................................................................................................................
ROC_cv <- function(dat,formula,method){
  #define cross validation parameters
  set.seed(1)
  ctrl <- trainControl(
    method = "repeatedcv", ## cross validation
    number = 10,   ## 10-fold
    repeats = 5,  # repeat k-fold cv 5 times
    summaryFunction = twoClassSummary, ## NEW
    classProbs = TRUE, # IMPORTANT
    verboseIter = F,
    savePredictions = T # important for plotting
  )
  
  # logit regression with cross validation
  logregFit <- train(formula, 
                     dat, 
                     method= method, 
                     preProc=c("center", "scale"), 
                     trControl=ctrl)
  
  list(
    Fit = logregFit,
    print=print(logregFit),
    summary = summary(logregFit),
    conf_mat = confusionMatrix(data = logregFit$pred$pred, reference = logregFit$pred$obs)
  )
  
  
}




# 18 function for Ramona....................................................................................................................
boxplot_Ramona <- function(data,
                           var,
                           width = 0.5,
                           outlier.shape = NA,
                           width_jitter = 0.15,
                           size_jitter = 1,
                           legend.position ="none",
                           col1= "darkgrey",
                           col2= "cornflowerblue",
                           col3= "forestgreen", 
                           significance = TRUE){
  p <- data %>% 
    filter(!is.na(get(var))) %>%
    ggplot(aes(get(var), log_exp, fill = get(var))) +
    geom_boxplot(width= width, outlier.shape = outlier.shape) +
    geom_jitter(width = width_jitter, size= size_jitter) +
    theme_bw() + 
    theme(legend.position = legend.position) + 
    xlab(var) +
    ylab("")+
    scale_y_continuous(breaks = seq(2,10,2), limits = c(1.5,11.8)) +
    scale_fill_manual(values = c(col1, col2, col3)) 
  
  p + if(significance == TRUE){stat_compare_means(label.x=1.1, label.y = 10.5, size=4)}
  
}






##
#calculate ANOVA for each column 
ANOVA_list <- function(data, names, covariate,group.var){
  ANOVA <- sapply(c(1:length(names)), function(x){
    formula <- as.formula(paste(names[x],"~", covariate,"*", group.var))
    data %>% anova_test(formula)
  })
  list <- c()
  for(i in 1:length(names)){
    list[[i]] <- ANOVA[,i] %>% as.data.frame()
  }
  names(list) <- names
  return(list)
}



# geoMean
gm_mean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}


# facet limits
facet_limits <-function(data, y, group){
  
  dat <- data.table(data)
  
  dat[,y_min := eval(parse(text = y))*0.5, by = eval(parse(text = group))]
  dat[,y_max:= eval(parse(text = y))*1.25, by = eval(parse(text = group))]
  
  dat <- as_tibble(dat)
  return(dat)
}






#....................................................................................................................  
## function to calculate confidence intervals based on t-distribution if the population sd is unknown and we have to use the 
# sample sd "s" 
# n: elements in sample
# t: t quantile (analog to 1.96 for 95% ci in normal distribution)
# m: margin of error
confInt <- function(x, conf = 0.95){
  n <- length(x)
  t <- qt((1+conf)/2, df = n -1)
  s <- sd(x)
  se <- s/sqrt(n)
  m <- t*se
  return(m)
}

