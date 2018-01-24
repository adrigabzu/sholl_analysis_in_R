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

############ Installing needed packages if missing ############################
list.of.packages <-
  c("data.table",
    "lmerTest",
    "nlme",
    "reshape2",
    "ggplot2",
    "plyr",
    "dplyr",
    "MESS")
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
library(dplyr)
library(MESS)

options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))

########### Function to generate the error bars for the sholl profile ###########
data_se <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(mean = mean(x[[col]], na.rm = TRUE),
      se = sd(x[[col]], na.rm = TRUE) / sqrt(length(x[[col]])))
  }
  data_sum <- ddply(data, groupnames, .fun = summary_func,
                    varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

############ Data parsing ################################

# Open files, CHECK THE DATA FORMAT
input_file = "example_data.tsv"
df_toanalyze = read.csv(input_file, sep = "\t")

# CHANGE THE VARIABLE NAMES ACCORDING TO THE COLUMN NAMES
# or name them:
# mouse_id. cell_id, radius, intersections, condition
colnames(df_toanalyze)

# Factorize categorical values
df_toanalyze$radius <- as.factor(df_toanalyze$radius)
df_toanalyze$mouse_id <- as.factor(df_toanalyze$mouse_id)
df_toanalyze$cell_id <- as.factor(df_toanalyze$cell_id)


################# Sholl profile plot and analyisis ######################
# Calculate SD for error bars
sholl_plot = data_se(
  data = df_toanalyze,
  varname = "intersections",
  groupnames = c("condition", "radius")
)

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
    aes(ymin = intersections - se, ymax = intersections + se),
    width = .2,
    position = position_dodge(0.05)
  ) +
  # Change here the range of the y-axis
  scale_y_continuous(breaks = seq(0, 50, by = 5), limits = c(0, 50)) +
  # Change the labels if necessary
  ylab("Dendritic Intersections") + xlab("Distance from the soma") +
  # Change the colors of the lines here
  theme_classic() + scale_color_manual(values = c("red", "green"))


# Overall test for differences between the conditions using mixed effects model
cond_test_corAR <- nlme::lme(
  intersections ~  1 + condition,
  data = df_toanalyze,
  random = ~ 1 | cell_id / radius ,
  correlation = corAR1()
)
summary(cond_test_corAR)
anova(cond_test_corAR)

# Sholl profile analysis with mixed effects models per radius
me_per_radius <-
  lmer(
    intersections ~ 1 + radius  + condition + radius:condition +  (1 |
                                                                     mouse_id / cell_id),
    df_toanalyze,
    REML = FALSE
  )

test_summ = summary(me_per_radius, ddf = "Satterthwaite")
print(test_summ, correlation = FALSE)

################# p-values adjustement ###########################
pvals = test_summ$coefficients[, 5]
cond_pvals = pvals[grep("condition", names(pvals))]
padj <- as.data.frame(p.adjust(cond_pvals , method = "BH"))
padj$significance <-
  symnum(
    padj[, 1],
    corr = FALSE,
    na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    symbols = c("***", "**", "*", " ")
  )

print("Adjusted p-values: ")
padj

################# AUC of sholl profile analysis ##################
df_auc <- df_toanalyze %>%
  group_by(mouse_id, cell_id, condition) %>%
  summarize(AUC = auc(intersections, radius, type = "spline"))

# Boxplot of AUC
ggplot(df_auc,
       aes(y = AUC,
           x = condition,
           fill  = condition)) +
  geom_boxplot() +
  ylab("AUC of Sholl profile") +
  # Change colors here
  theme_classic() + scale_fill_manual(values = c("red", "green"))


# Mixed-effect model of AUC
m <- nlme::lme(AUC ~  1 + condition,
               data = df_auc,
               random = ~ 1 | mouse_id / cell_id)
summary(m)
