##################################################################################
#       Statistical analysis of Sholl profiles based on mixed-effect models     
# 
# Adrian Gabriel Zucco
# Center for Translational Neuromedicine
# University of Copenhagen
#
# Adapted from:
# Wilson, M. D., Sethi, S., Lein, P. J. & Keil, K. P.                            
# Valid statistical approaches for analyzing sholl data: Mixed effects           
# versus simple linear models. Journal of Neuroscience Methods 279, 33-43 (2017).
#
##################################################################################

#installing needed packages if missing
list.of.packages <-
  c(
    "data.table",
    "lmerTest",
    "nlme",
    "reshape2",
    "ggplot2",
    "plyr"
  )
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) {
  install.packages(new.packages)
}

library(data.table)
library(lmerTest)
library(nlme)
library(reshape2)
library(ggplot2)
library(plyr)


options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))

########### Function to generate the error bars ##############
data_sd <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun = summary_func,
                    varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

############ Data handling ################################

# Open files
input_file = "example_data.tsv"
df_toanalyze = read.csv(input_file, sep = "\t")
colnames(df_toanalyze)

df_toanalyze$radius <- as.factor(df_toanalyze$radius)
df_toanalyze$mouse_id <- as.factor(df_toanalyze$mouse_id)
df_toanalyze$cell_id <- as.factor(df_toanalyze$cell_id)

# Open your tsv / csv file
# Calculate SD for error bars
sholl_plot = data_sd(df_toanalyze,
                          varname = "intersections",
                          groupnames = c("condition", "radius"))


# Sholl profile plot with error bars, adapt variables to your column names
ggplot(sholl_plot,
       aes(
         y = intersections,
         x = radius,
         color = condition ,
         group = condition
       )) +
  geom_line() +
  geom_point() +
  geom_errorbar(
    aes(ymin = intersections - sd, ymax = intersections + sd),
    width = .2,
    position = position_dodge(0.05)
  ) +
  scale_y_continuous(breaks = seq(0, 50, by = 5), limits = c(0, 50)) + # Change here the range of the y-axis
  ylab("Intersections") + xlab("Radius") +
  theme_classic() + scale_color_manual(values = c("red", "green")) # Change the colors of the lines here


# Overall test for differences between the conditions using mixed effects model
m <- nlme::lme(intersections ~  1 + condition,
         data = df_toanalyze, random = ~ 1 | cell_id/radius , correlation=corAR1())
anova(m)

# Sholl profile analysis with mixed effects models per radius
m <-
  lmer(
    intersections ~ 1 + radius  + condition + radius:condition +  (1 | mouse_id/cell_id),
    df_toanalyze,
    REML = FALSE
  )
print(summary(m, ddf = "Satterthwaite"), correlation = FALSE)


############### p-values adjustement ###########################
pvals = summary(m, ddf = "Satterthwaite")$coefficients[, 5]
cond_pvals = pvals[grep("condition",names(pvals))]
padj <- as.data.frame(p.adjust( cond_pvals , method = "BH"))
padj$significance <- symnum(padj[,1], corr = FALSE, na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
padj

