## ---- warning=FALSE, error=FALSE-----------------------------------------
drug_class = 'PI' # Possible drug types are 'PI', 'NRTI', and 'NNRTI'. 


## ---- warning=FALSE, error=FALSE-----------------------------------------
base_url = 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
gene_url = paste(base_url, 'DATA', paste0(drug_class, '_DATA.txt'), sep='/')
tsm_url = paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep='/')

gene_df = read.delim(gene_url, na.string = c('NA', ''), stringsAsFactors = FALSE)
tsm_df = read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
names(tsm_df) = c('Position', 'Mutations')


## ---- warning=FALSE, error=FALSE-----------------------------------------
head(gene_df, n=6)


## ---- warning=FALSE, error=FALSE-----------------------------------------
head(tsm_df, n=6)


## ---- warning=FALSE, error=FALSE-----------------------------------------
# Returns rows for which every column matches the given regular expression.
grepl_rows <- function(pattern, df) {
  cell_matches = apply(df, c(1,2), function(x) grepl(pattern, x))
  apply(cell_matches, 1, all)
}

pos_start = which(names(gene_df) == 'P1')
pos_cols = seq.int(pos_start, ncol(gene_df))
valid_rows = grepl_rows('^(\\.|-|[A-Zid]+)$', gene_df[,pos_cols])
gene_df = gene_df[valid_rows,]


## ---- warning=FALSE, error=FALSE-----------------------------------------
# Flatten a matrix to a vector with names from concatenating row/column names.
flatten_matrix <- function(M, sep='.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep=sep)))
  x
}

# Construct preliminary design matrix.
muts = c(LETTERS, 'i', 'd')
X = outer(muts, as.matrix(gene_df[,pos_cols]), Vectorize(grepl))
X = aperm(X, c(2,3,1))
dimnames(X)[[3]] <- muts
X = t(apply(X, 1, flatten_matrix))
mode(X) <- 'numeric'

# Remove any mutation/position pairs that never appear in the data.
X = X[,colSums(X) != 0]

# Extract response matrix.
Y = gene_df[,4:(pos_start-1)]


## ---- warning=FALSE, error=FALSE-----------------------------------------
library(DT)
datatable(data.frame(X)[1:10, ], options = list(scrollX=T, pageLength = 10))


## ---- warning=FALSE, error=FALSE-----------------------------------------
head(Y, n=6)


## ---- warning=FALSE, error=FALSE, message=FALSE--------------------------
# Load the leaps library and tidyverse set of libraries
library(leaps)
library(tidyverse)

# Initialize empty lists to hold
#   - regsub: The result of using regsub on the linear regression of
#             resistance response vector against the design matrix 
#             for each drug
#   - regsub_mdl_cp: The model selected by regsubsets() with the lowest 
#                    Mallow's Cp for each drug
#   - regsub_mdl_bic: The model selected by regsubsets() with the BIC
#                     value for each drug
#   - opt_mdl_pos_bic: Hold the list of positions deemed significant
#                      in the linear regression of resistance against
#                      the positions indicated in regsub_mdl_cp
#   - opt_mdl_pos_cp: Hold the list of positions deemed significant
#                     in the linear regression of resistance against
#                     the positions indicated in regsub_mdl_bic
regsub = list()
regsub_mdl_cp = list()
regsub_mdl_bic = list()
opt_mdl_pos_bic = list()
opt_mdl_pos_cp = list()

# Loop over all drugs in the response matrix Y
for(drug in names(Y)){
  # Deal with the response vector corresponding to the current drug
  y = Y[[drug]]
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  x = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  x = x[,colSums(x) >= 3]
  # Remove duplicate predictors.
  x = x[,colSums(abs(cor(x)-1) < 1e-4) == 1]
  
  # Select the position-mutation pairs associated with drug-resistance using regsubsets() and the backward 
  # stepwise method
  regsub[[drug]] = summary(regsubsets(y~., data=data.frame(x), method='backward', nvmax=ncol(x), really.big=TRUE))
  
  # Select the entry in outmat with the lowest BIC value from the output of summary(regsubsets()) above
  regsub_mdl_bic[[drug]] = regsub[[drug]]$outmat[which.min(regsub[[drug]]$bic),]
  # Obtain the position-mutation pairs selected by regsubsets() for the entry in outmat with the lowest BIC value
  regsub_mdl_pos_bic = names(regsub_mdl_bic[[drug]][regsub_mdl_bic[[drug]] == '*'])
  # From the coefficient vector extracted above, remove the mutation from the position mutation pair, make all 
  # the resulting positions numeric values, and finally only keep the unique values of the resulting vector
  opt_mdl_pos_bic[[drug]] =  regsub_mdl_pos_bic %>%
    substring(2,3) %>%
    as.numeric() %>%
    unique()
  # Remove NA values from the vector of positions
  opt_mdl_pos_bic[[drug]] = opt_mdl_pos_bic[[drug]][!is.na(opt_mdl_pos_bic[[drug]])]
  
  # Select the entry in outmat with the lowest Mallows Cp value from the output of summary(regsubsets()) above
  regsub_mdl_cp[[drug]] = regsub[[drug]]$outmat[which.min(regsub[[drug]]$cp),]
  # Obtain the position-mutation pairs selected by regsubsets() for the entry in outmat with the lowest Mallows 
  # Cp value
  regsub_mdl_pos_cp = names(regsub_mdl_cp[[drug]][regsub_mdl_cp[[drug]] == '*'])
  # From the coefficient vector extracted above, remove the mutation from the position mutation pair, make all 
  # the resulting positions numeric values, and finally only keep the unique values of the resulting vector
  opt_mdl_pos_cp[[drug]] <- regsub_mdl_pos_cp %>%
    substring(2,3) %>%
    as.numeric() %>%
    unique()
  # Remove NA values from the vector of positions
  opt_mdl_pos_cp[[drug]] = opt_mdl_pos_cp[[drug]][!is.na(opt_mdl_pos_cp[[drug]])]
}


## ------------------------------------------------------------------------
opt_mdl_pos_bic


## ------------------------------------------------------------------------
opt_mdl_pos_cp


## ---- warning=FALSE, error=FALSE, message=FALSE--------------------------
# Load the cowplot library
library(cowplot)

# Initialize lists to hold
#   - drug_lm: The linear regressions between the resistance response vector
#              and design matrix for each drug
#   - drug_lm_plt_arr: The diagnostic plots for each linear regression created
drug_lm = list()
drug_lm_plt_arr = list()

# Loop over all drugs in the response matrix Y
for(drug in names(Y)){
  # Deal with the response vector corresponding to the current drug
  y = Y[[drug]]
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  x = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  x = x[,colSums(x) >= 3]
  # Remove duplicate predictors.
  x = x[,colSums(abs(cor(x)-1) < 1e-4) == 1]
  
  # Find the linear regression of the resistance response vector y against the design matrix 
  drug_lm[[drug]] <- lm(y~., data=data.frame(x))
  
  # Create the studentized deleted residuals qq plot for the linear regression
  p1 <- ggplot(data.frame('studentized_deleted_residuals'=rstudent(drug_lm[[drug]])), aes(sample=studentized_deleted_residuals)) +
    stat_qq() +
    stat_qq_line(color='red') +
    labs(title='QQ Plot of Studentized Deleted\nResiduals',
         x='Theoretical Quantiles',
         y='Sample Quantiles') +
    theme_light()
  
  # Create the studentized deleted residuals histogram for the linear regression
  p2 <- ggplot(data.frame('resistance'=y**0.06060606, 
                    'studentized_deleted_residuals'=rstudent(drug_lm[[drug]])), aes(x=studentized_deleted_residuals)) +
    geom_histogram(aes(y=..density..), binwidth=sd(rstudent(drug_lm[[drug]])), color='black') +
    stat_function(fun=dnorm, args=list(mean=mean(rstudent(drug_lm[[drug]])), 
                                       sd=sd(rstudent(drug_lm[[drug]]))), color='red', size=1) + 
    labs(title='Histogram of Studentized Deleted\nResiduals',
         subtitle='With Normal PDF Curve Overlaid in Red',
         x='Studentized Deleted Residuals',
         y='Frequency') +
    theme_light()
  
  # Create the studentized deleted residuals vs predicted values for the linear regression
  p3 <- ggplot(data.frame('predicted'=drug_lm[[drug]]$fitted.values,
                    'studentized_deleted_residuals'=rstudent(drug_lm[[drug]])), 
         aes(x=predicted, y=studentized_deleted_residuals)) +
    geom_point() +
    geom_hline(aes(yintercept=0), color='red') +
    labs(title='Studentized Deleted Residuals\nvs Predicted Resistance',
         x='Predicted Resistance',
         y='Studentized Deleted\nResiduals') +
    theme_light()
  
  # Create the studentized deleted residuals line plot for the linear regression
  p4 <- ggplot(data.frame('studentized_deleted_residuals'=rstudent(drug_lm[[drug]])),
         aes(x=1:length(studentized_deleted_residuals), y=studentized_deleted_residuals)) +
    geom_point(pch=21, cex=3) +
    geom_line(color='red') +
    geom_hline(aes(yintercept=0)) +
    labs(title='Line Plot of Studentized Deleted\nResiduals',
         y='Studentized Deleted\nResiduals',
         x='') +
    theme_light()
  
  # Arrange all 4 diagnostic plots created above in a 2x2 grid and store it in list drug_lm_plt_arr
  drug_lm_plt_arr[[drug]] <- plot_grid(p1, p2, p3, p4, nrow=2, ncol=2) +
    draw_figure_label(paste(c(drug, '\n'), collapse=''), position = "top", size=12, fontface='bold')
}


## ------------------------------------------------------------------------
drug_lm_plt_arr[1]


## ------------------------------------------------------------------------
drug_lm_plt_arr[2]


## ------------------------------------------------------------------------
drug_lm_plt_arr[3]


## ------------------------------------------------------------------------
drug_lm_plt_arr[4]


## ------------------------------------------------------------------------
drug_lm_plt_arr[5]


## ------------------------------------------------------------------------
drug_lm_plt_arr[6]


## ------------------------------------------------------------------------
drug_lm_plt_arr[7]


## ---- warning=FALSE, error=FALSE, message=FALSE--------------------------
# Load the MASS library
library(MASS)

# Initialize lists to hold
#   - bac.cox: The result of using boxcox on the linear regression of the
#              resistance response vector gainst the design matrix for each drug
#   - bac.lambda: The optimal lambda found for each linear regression of the
#                 resistance response vector against the design matrix for each drug
#   - bc_drug_lm: The linear regressions between the Box Cox transformed resistance 
#                 response vector and design matrix for each drug
#   - bc_drug_lm_plt_arr: The diagnostic plots for each linear regression created
bac.box = list()
bac.lambda=list()
bc_drug_lm = list()
bc_drug_lm_plt_arr = list()

# Loop through all drugs in response matrix Y
for(drug in names(Y)){
  # Deal with the response vector corresponding to the current drug
  y = Y[[drug]]
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  x = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  x = x[,colSums(x) >= 3]
  # Remove duplicate predictors.
  x = x[,colSums(abs(cor(x)-1) < 1e-4) == 1]
  
  # Use the function boxcox on the linear regression of the resistance response vector y 
  # against the design matrix for each drug
  bac.box[[drug]] = boxcox(lm(y~., data=data.frame(x)), plotit=FALSE)
  # Find the optimal lambda to transform the resistance response vector
  bac.lambda[[drug]] = bac.box[[drug]]$x[which(bac.box[[drug]]$y==max(bac.box[[drug]]$y))]
  
  # If the optimal lambda is 0, then perform a linear regression of the natural log of the resistance response 
  # vector against the design matrix for each drug. Otherwise, perform a linear regression of the 
  # resistance response vector to the power of lambda against the design matrix for each drug.
  if(bac.lambda[[drug]] == 0){
    bc_drug_lm[[drug]] <- lm(log(y)~., data=data.frame(x))
  } else {
    bc_drug_lm[[drug]] <- lm((y**bac.lambda[[drug]])~., data=data.frame(x))
  }

  # Create the studentized deleted residuals qq plot for the transformed linear regression
  p1 <- ggplot(data.frame('studentized_deleted_residuals'=rstudent(bc_drug_lm[[drug]])), aes(sample=studentized_deleted_residuals)) +
    stat_qq() +
    stat_qq_line(color='red') +
    labs(title='QQ Plot of Studentized Deleted\nResiduals',
         x='Theoretical Quantiles',
         y='Sample Quantiles') +
    theme_light()
  
  # Create the studentized deleted residuals histogram for the transformed linear regression
  p2 <- ggplot(data.frame('resistance'=y**0.06060606, 
                    'studentized_deleted_residuals'=rstudent(bc_drug_lm[[drug]])), aes(x=studentized_deleted_residuals)) +
    geom_histogram(aes(y=..density..), binwidth=sd(rstudent(bc_drug_lm[[drug]])), color='black') +
    stat_function(fun=dnorm, args=list(mean=mean(rstudent(bc_drug_lm[[drug]])), 
                                       sd=sd(rstudent(bc_drug_lm[[drug]]))), color='red', size=1) + 
    labs(title='Histogram of Studentized Deleted\nResiduals',
         subtitle='With Normal PDF Curve Overlaid in Red',
         x='Studentized Deleted Residuals',
         y='Frequency') +
    theme_light()
  
  # Create the studentized deleted residuals vs predicted values for the transformed linear regression
  p3 <- ggplot(data.frame('predicted'=bc_drug_lm[[drug]]$fitted.values,
                    'studentized_deleted_residuals'=rstudent(bc_drug_lm[[drug]])), 
         aes(x=predicted, y=studentized_deleted_residuals)) +
    geom_point() +
    geom_hline(aes(yintercept=0), color='red') +
    labs(title='Studentized Deleted Residuals\nvs Predicted Resistance',
         x='Predicted Resistance',
         y='Studentized Deleted\nResiduals') +
    theme_light()
  
  # Create the studentized deleted residuals line plot for the transformed linear regression
  p4 <- ggplot(data.frame('studentized_deleted_residuals'=rstudent(bc_drug_lm[[drug]])),
         aes(x=1:length(studentized_deleted_residuals), y=studentized_deleted_residuals)) +
    geom_point(pch=21, cex=3) +
    geom_line(color='red') +
    geom_hline(aes(yintercept=0)) +
    labs(title='Line Plot of Studentized Deleted\nResiduals',
         y='Studentized Deleted\nResiduals',
         x='') +
    theme_light()
  
  # Arrange all 4 diagnostic plots created above in a 2x2 grid and store it in list bc_drug_lm_plt_arr
  bc_drug_lm_plt_arr[[drug]] <- plot_grid(p1, p2, p3, p4, nrow=2, ncol=2) +
    draw_figure_label(paste(c(drug, '\n'), collapse=''), position = "top", size=12, fontface='bold')
}


## ------------------------------------------------------------------------
bc_drug_lm_plt_arr[1]


## ------------------------------------------------------------------------
bc_drug_lm_plt_arr[2]


## ------------------------------------------------------------------------
bc_drug_lm_plt_arr[3]


## ------------------------------------------------------------------------
bc_drug_lm_plt_arr[4]


## ------------------------------------------------------------------------
bc_drug_lm_plt_arr[5]


## ------------------------------------------------------------------------
bc_drug_lm_plt_arr[6]


## ------------------------------------------------------------------------
bc_drug_lm_plt_arr[7]


## ---- warning=FALSE, error=FALSE, message=FALSE--------------------------
# Initialize empty lists to hold
#   - bc_regsub: The result of using regsub on the Box Cox transformed
#                linear regression of resistance response vector against 
#                the design matrix for each drug
#   - bc_regsub_mdl_cp: The model selected by regsubsets() with the lowest 
#                       Mallow's Cp for each drug
#   - bc_regsub_mdl_bic: The model selected by regsubsets() with the BIC
#                        value for each drug
#   - bc_opt_mdl_pos_bic: Hold the list of positions deemed significant
#                         in the linear regression of resistance against
#                         the positions indicated in bc_regsub_mdl_cp
#   - bc_opt_mdl_pos_cp: Hold the list of positions deemed significant
#                        in the linear regression of resistance against
#                        the positions indicated in bc_regsub_mdl_bic
bc_regsub = list()
bc_regsub_mdl_cp = list()
bc_regsub_mdl_bic = list()
bc_opt_mdl_pos_bic = list()
bc_opt_mdl_pos_cp = list()

# Loop through all drugs in response matrix Y
for(drug in names(Y)){
  
  # Deal with the response vector corresponding to the current drug
  y = Y[[drug]]
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  x = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  x = x[,colSums(x) >= 3]
  # Remove duplicate predictors.
  x = x[,colSums(abs(cor(x)-1) < 1e-4) == 1]
  
  # If the optimal lambda is 0, then use regsubsets() on the linear regression of the natural log of the
  # resistance response vector against the design matrix for each drug. Otherwise, use regsubsets() on 
  # the linear regression of the resistance response vector to the power of lambda against the design 
  # matrix for each drug.
  if(bac.lambda[[drug]] == 0){
    bc_regsub[[drug]] = summary(regsubsets(log(y)~., data=data.frame(x), method='backward', nvmax=ncol(x), really.big=TRUE))
  } else {
    bc_regsub[[drug]] = summary(regsubsets((y**bac.lambda[[drug]])~., data=data.frame(x), method='backward', nvmax=ncol(x), really.big=TRUE))
  }
  
  # Select the entry in outmat with the lowest BIC value from the output of summary(regsubsets()) above
  bc_regsub_mdl_bic[[drug]] = bc_regsub[[drug]]$outmat[which.min(bc_regsub[[drug]]$bic),]
  # Obtain the position-mutation pairs selected by regsubsets() for the entry in outmat with the lowest BIC value
  bc_regsub_mdl_pos_bic = names(bc_regsub_mdl_bic[[drug]][bc_regsub_mdl_bic[[drug]] == '*'])
  # From the coefficient vector extracted above, remove the mutation from the position mutation pair, make all the 
  # resulting positions numeric values, and finally only keep the unique values of the resulting vector 
  bc_opt_mdl_pos_bic[[drug]] <- bc_regsub_mdl_pos_bic %>%
    substring(2,3) %>%
    as.numeric() %>%
    unique()
  # Remove NA values from the vector of positions
  bc_opt_mdl_pos_bic[[drug]] = bc_opt_mdl_pos_bic[[drug]][!is.na(bc_opt_mdl_pos_bic[[drug]])]
  
  # Select the entry in outmat with the lowest Mallows Cp value from the output of summary(regsubsets()) above
  bc_regsub_mdl_cp[[drug]] = bc_regsub[[drug]]$outmat[which.min(bc_regsub[[drug]]$cp),]
  # Obtain the position-mutation pairs selected by regsubsets() for the entry in outmat with the lowest Mallows 
  # Cp value
  bc_regsub_mdl_pos_cp = names(bc_regsub_mdl_cp[[drug]][bc_regsub_mdl_cp[[drug]] == '*'])
  # From the coefficient vector extracted above, remove the mutation from the position mutation pair, make all 
  # the resulting positions numeric values, and finally only keep the unique values of the resulting vector 
  bc_opt_mdl_pos_cp[[drug]] <- bc_regsub_mdl_pos_cp %>%
    substring(2,3) %>%
    as.numeric() %>%
    unique()
  # Remove NA values from the vector of positions
  bc_opt_mdl_pos_cp[[drug]] = bc_opt_mdl_pos_cp[[drug]][!is.na(bc_opt_mdl_pos_cp[[drug]])]
}


## ------------------------------------------------------------------------
bc_opt_mdl_pos_bic


## ------------------------------------------------------------------------
bc_opt_mdl_pos_cp


## ---- warning=FALSE, message=FALSE---------------------------------------
# Load the car library
library(car)
# Initialize lists to hold:
#   - init_vif: The result of calling vif() on the Box Cox transformed
#          linear regressions for each drug
#   - init_max_vif: Hold the max vif value from init_vif
#   - init_mean_vif: Hold the mean vif value from init_vif
#   - bc_vif: The result of calling vif() on the Box Cox transformed
#          linear regressions for each drug
#   - bc_max_vif: Hold the max vif value from bc_vif
#   - bc_mean_vif: Hold the mean vif value from bc_vif
init_vif = list()
init_max_vif = list()
init_mean_vif = list()

bc_vif = list()
bc_max_vif = list()
bc_mean_vif = list()

# Loop through all drugs in the response matrix Y
for(drug in names(Y)){
  # Deal with the response vector corresponding to the current drug
  y = Y[[drug]]
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  x = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  x = x[,colSums(x) >= 3]
  # Remove duplicate predictors.
  x = x[,colSums(abs(cor(x)-1) < 1e-4) == 1]
  
  # Use vif() on the linear regression of the resistance response vector against the design matrix for each drug
  init_vif[[drug]] = vif(lm(y~., data.frame(x)))
  # Determine the maximum VIF for the Box Cox transformed linear regression of the current drug
  init_max_vif[[drug]] = round(max(init_vif[[drug]]),2)
  # Determine the mean of the VIFs for the Box Cox transformed linear regression of the current drug
  init_mean_vif[[drug]] = round(mean(init_vif[[drug]]),2)
  
  # If the optimal lambda is 0, then use vif() on the linear regression of the natural log of the resistance 
  # response vector against the design matrix for each drug. Otherwise, use vif() on the linear regression 
  # of the resistance response vector to the power of lambda against the design matrix for each drug.
  if(bac.lambda==0){
    bc_vif[[drug]] = vif(lm(log(y)~., data=data.frame(x)))
  } else {
    bc_vif[[drug]] = vif(lm((y**bac.lambda[[drug]])~., data=data.frame(x)))
  }
  # Determine the maximum VIF for the Box Cox transformed linear regression of the current drug
  bc_max_vif[[drug]] = round(max(bc_vif[[drug]]),2)
  # Determine the mean of the VIFs for the Box Cox transformed linear regression of the current drug
  bc_mean_vif[[drug]] = round(mean(bc_vif[[drug]]),2)
}

# Create a dataframe for the max and mean of the init_vif values for each drug
init_vif_res = tibble('Drug'=names(Y),
                    'Sample Max VIF'=unlist(init_max_vif),
                    'Sample Mean VIF'=unlist(init_mean_vif))
# Show the dataframe bc_vif_res as a datatable
init_col_dt = datatable(init_vif_res)

# Create a dataframe for the max and mean of the bc_vif values for each drug
bc_vif_res = tibble('Drug'=names(Y),
                    'Sample Max VIF'=unlist(bc_max_vif),
                    'Sample Mean VIF'=unlist(bc_mean_vif))
# Show the dataframe bc_vif_res as a datatable
bc_col_dt = datatable(bc_vif_res)


## ------------------------------------------------------------------------
init_col_dt


## ------------------------------------------------------------------------
bc_col_dt


## ------------------------------------------------------------------------
# Initialize lists to hold
#   - init_total_cp_pred_pos: All positions selected by the regsubsets() subset with the
#                             lowest Mallows Cp value before Box Cox transform
#   - init_total_BIC_pred_pos: All positions selected by the regsubsets() subset with the
#                              lowest BIC value before Box Cox transform
#   - bc_total_cp_pred_pos: All positions selected by the regsubsets() subset with the
#                           lowest Mallows Cp value after Box Cox transform
#   - bc_total_BIC_pred_pos: All positions selected by the regsubsets() subset with the
#                            lowest BIC value after Box Cox transform
#   - init_regsubs_res: Holds lists of the number of correct and incorrect position selections 
#                       before Box Cox transformation
#   - bc_regsubs_res: Holds lists of the number of correct and incorrect position selections 
#                     after Box Cox transformation
init_total_cp_pred_pos = c()
init_total_bic_pred_pos = c()
bc_total_cp_pred_pos = c()
bc_total_bic_pred_pos = c()
init_regsubs_res = list()
bc_regsubs_res = list()

# Loop through all drugs in response matrix Y
for(drug in names(Y)){
  # Check which positions selected by the regsubsets() subset with the lowest Mallows
  # Cp value before Box Cox transform are in the list of ground true positions
  opt_mdl_pos_cp_chk = opt_mdl_pos_cp[[drug]] %in% tsm_df$Position
  cp_correct = sum(opt_mdl_pos_cp_chk)
  cp_incorrect = length(opt_mdl_pos_cp_chk) - cp_correct
  
  # Add the positions deemed significant for this drug from the lowest Mallow's Cp value 
  # subset after Box Cox transform to the list of all positions deemed significant from
  # the lowest Mallow's Cp value subsets 
  init_total_cp_pred_pos = c(init_total_cp_pred_pos, opt_mdl_pos_cp[[drug]])
  
  # Check which positions selected by the regsubsets() subset with the lowest BIC
  # value before Box Cox transform are in the list of ground true positions
  opt_mdl_pos_bic_chk = opt_mdl_pos_bic[[drug]] %in% tsm_df$Position
  bic_correct = sum(opt_mdl_pos_bic_chk)
  bic_incorrect = length(opt_mdl_pos_bic_chk) - bic_correct
  
  # Add the positions deemed significant for this drug from the lowest BIC value 
  # value subset before Box Cox transform to the list of all positions deemed
  # significant from the lowest BIC value subsets 
  init_total_bic_pred_pos = c(init_total_bic_pred_pos, opt_mdl_pos_bic[[drug]])
  
  # Create a vector holding the number of correct and incorrect selections for each criterion
  # before Box Cox transform and store it in the init_regsub_res list
  init_regsubs_res[[drug]] = c('Mallows_Cp|Correct'=cp_correct, 'Mallows_Cp|Incorrect'=cp_incorrect,
                              'BIC|Correct'=bic_correct, 'BIC|Incorrect'=bic_incorrect)
  
  # Check which positions selected by the regsubsets() subset with the lowest Mallows
  # Cp value after Box Cox transform are in the list of ground true positions
  bc_opt_mdl_pos_cp_chk = bc_opt_mdl_pos_cp[[drug]] %in% tsm_df$Position
  bc_cp_correct = sum(bc_opt_mdl_pos_cp_chk)
  bc_cp_incorrect = length(bc_opt_mdl_pos_cp_chk) - bc_cp_correct
  
  # Add the positions deemed significant for this drug from the lowest Mallow's Cp value 
  # value subset after Box Cox transform to the list of all positions deemed significant 
  # from the lowest Mallow's Cp value subsets 
  bc_total_cp_pred_pos = c(bc_total_cp_pred_pos, bc_opt_mdl_pos_cp[[drug]])
  
  # Check which positions selected by the regsubsets() subset with the lowest BIC
  # value before Box Cox transform are in the list of ground true positions
  bc_opt_mdl_pos_bic_chk = bc_opt_mdl_pos_bic[[drug]] %in% tsm_df$Position
  bc_bic_correct = sum(bc_opt_mdl_pos_bic_chk)
  bc_bic_incorrect = length(bc_opt_mdl_pos_bic_chk) - bc_bic_correct
  
  # Add the positions deemed significant for this drug from the lowest BIC value 
  # value subset after Box Cox transform to the list of all positions deemed
  # significant from the lowest BIC value subsets 
  bc_total_bic_pred_pos = c(bc_total_bic_pred_pos, bc_opt_mdl_pos_bic[[drug]])
  
  # Create a vector holding the number of correct and incorrect selections for each criterion
  # after Box Cox transform and store it in the init_regsub_res list
  bc_regsubs_res[[drug]] = c('Mallows_Cp|Correct'=bc_cp_correct, 'Mallows_Cp|Incorrect'=bc_cp_incorrect,
                              'BIC|Correct'=bc_bic_correct, 'BIC|Incorrect'=bc_bic_incorrect)
}

# Check which positions selected by the regsubsets() subsets selection with the lowest
# Mallow's Cp values before Box Cox transform are in the list of ground true positions.
init_total_cp_pred_pos_chk = unique(init_total_cp_pred_pos) %in% tsm_df$Position
init_total_cp_correct = sum(init_total_cp_pred_pos_chk)
init_total_cp_incorrect = length(init_total_cp_pred_pos_chk) - init_total_cp_correct

# Check which positions selected by the regsubsets() subsets selection with the lowest
# BIC values before Box Cox transform are in the list of ground true positions.
init_total_bic_pred_pos_chk = unique(init_total_bic_pred_pos) %in% tsm_df$Position
init_total_bic_correct = sum(init_total_bic_pred_pos_chk)
init_total_bic_incorrect = length(init_total_bic_pred_pos_chk) - init_total_bic_correct

# Check which positions selected by the regsubsets() subsets selection with the lowest
# Mallow's Cp values after Box Cox transform are in the list of ground true positions.
bc_total_cp_pred_pos_chk = unique(bc_total_cp_pred_pos) %in% tsm_df$Position
bc_total_cp_correct = sum(bc_total_cp_pred_pos_chk)
bc_total_cp_incorrect = length(bc_total_cp_pred_pos_chk) - bc_total_cp_correct

# Check which positions selected by the regsubsets() subsets selection with the lowest
# BIC values after Box Cox transform are in the list of ground true positions.
bc_total_bic_pred_pos_chk = unique(bc_total_bic_pred_pos) %in% tsm_df$Position
bc_total_bic_correct = sum(bc_total_bic_pred_pos_chk)
bc_total_bic_incorrect = length(bc_total_bic_pred_pos_chk) - bc_total_bic_correct

# Create a tibble with the number of correct and incorrect selections before or after 
# Box Cox and Criterion (BIC or Mallow's Cp). Then add the ratio of correct selections to total
# selections.
tot_preds <- tibble('Model'=c('Before_Box_Cox','Before_Box_Cox','Before_Box_Cox','Before_Box_Cox',
                              'After_Box_Cox', 'After_Box_Cox', 'After_Box_Cox', 'After_Box_Cox'),
                    'Criterion'=c('BIC', 'BIC', 'Mallows_Cp', 'Mallows_Cp',
                                  'BIC', 'BIC', 'Mallows_Cp', 'Mallows_Cp'),
                    'Selection'=c('Correct', 'Incorrect','Correct', 'Incorrect',
                                   'Correct', 'Incorrect','Correct', 'Incorrect'),
                    'Total'=c(init_total_bic_correct, init_total_bic_incorrect,
                              init_total_cp_correct, init_total_cp_incorrect,
                              bc_total_bic_correct, bc_total_bic_incorrect,
                              bc_total_cp_correct, bc_total_cp_incorrect)) %>%
  mutate(Model=factor(Model, levels=c('Before_Box_Cox', 'After_Box_Cox')),
         Selection=factor(Selection, levels=c('Incorrect', 'Correct')))

# Reformat the tibble created above to be put in a datatable object later
tot_res_dt <- tot_preds %>%
  group_by(Model) %>%
  mutate(Crit_Pred = paste0(Criterion, '_', Selection)) %>%
  pivot_wider(id_cols=c(Model, Crit_Pred, Total), names_from=Crit_Pred, values_from=Total) %>%
  mutate(Mallows_Cp_Precision = round(Mallows_Cp_Correct/sum(Mallows_Cp_Correct, Mallows_Cp_Incorrect), 2),
         BIC_Precision = round(BIC_Correct/sum(BIC_Correct, BIC_Incorrect), 2)) %>%
  dplyr::select(Model, BIC_Correct, BIC_Incorrect, BIC_Precision, 
                Mallows_Cp_Correct, Mallows_Cp_Incorrect, Mallows_Cp_Precision)

# Create a dataframe with the number of correct and incorrect selections before Box Cox by drug 
# and Criterion (BIC or Mallow's Cp). Then add the ratio of correct selections to total
# selections.
init_regsub_res_df = data.frame(matrix(unlist(init_regsubs_res), ncol=2, byrow=TRUE)) %>%
  rename('Correct'=X1, 'Incorrect'=X2)
init_regsub_res_df$Criterion = c("Mallows_Cp", "BIC")
init_regsub_res_df$Drug = c('APV', 'APV', 'ATV', 'ATV', 'IDV', 'IDV', 'LPV', 'LPV', 'NFV', 'NFV', 'RTV', 'RTV', 'SQV', 'SQV')
init_regsub_res_df = init_regsub_res_df %>% 
  pivot_longer(c(-Drug, -Criterion), names_to='Selection', values_to='Total') %>%
  mutate(Selection=factor(Selection, levels=c('Incorrect', 'Correct')))

# Reformat the dataframe created above to be put in a datatable object later
init_regsub_res_dt <- init_regsub_res_df %>%
  group_by(Drug) %>%
  mutate(Crit_Pred = paste0(Criterion, '_', Selection)) %>%
  pivot_wider(id_cols=c(Drug, Crit_Pred, Total), names_from=Crit_Pred, values_from=Total) %>%
  mutate(Mallows_Cp_Precision = round(Mallows_Cp_Correct/sum(Mallows_Cp_Correct, Mallows_Cp_Incorrect), 2),
         BIC_Precision = round(BIC_Correct/sum(BIC_Correct, BIC_Incorrect), 2)) %>%
  dplyr::select(Drug, BIC_Correct, BIC_Incorrect, BIC_Precision, 
                Mallows_Cp_Correct, Mallows_Cp_Incorrect, Mallows_Cp_Precision)

# Create a dataframe with the number of correct and incorrect selections after Box Cox by drug 
# and Criterion (BIC or Mallow's Cp). Then add the ratio of correct selections to total
# selections.
bc_regsub_res_df = data.frame(matrix(unlist(bc_regsubs_res), ncol=2, byrow=TRUE)) %>%
  rename('Correct'=X1, 'Incorrect'=X2)
bc_regsub_res_df$Criterion = c("Mallows_Cp", "BIC")
bc_regsub_res_df$Drug = c('APV', 'APV', 'ATV', 'ATV', 'IDV', 'IDV', 'LPV', 'LPV', 'NFV', 'NFV', 'RTV', 'RTV', 'SQV', 'SQV')
bc_regsub_res_df = bc_regsub_res_df %>% 
  pivot_longer(c(-Drug, -Criterion), names_to='Selection', values_to='Total') %>%
  mutate(Selection=factor(Selection, levels=c('Incorrect', 'Correct')))

# Reformat the dataframe created above to be put in a datatable object later
bc_regsub_res_dt <- bc_regsub_res_df %>%
  group_by(Drug) %>%
  mutate(Crit_Pred = paste0(Criterion, '_', Selection)) %>%
  pivot_wider(id_cols=c(Drug, Crit_Pred, Total), names_from=Crit_Pred, values_from=Total) %>%
  mutate(Mallows_Cp_Precision = round(Mallows_Cp_Correct/sum(Mallows_Cp_Correct, Mallows_Cp_Incorrect), 2),
         BIC_Precision = round(BIC_Correct/sum(BIC_Correct, BIC_Incorrect), 2))  %>%
  dplyr::select(Drug, BIC_Correct, BIC_Incorrect, BIC_Precision, 
                Mallows_Cp_Correct, Mallows_Cp_Incorrect, Mallows_Cp_Precision)


## ---- warning=FALSE, error=FALSE, message=FALSE--------------------------
# Generate a stacked bar plot to visualize the correct and incorrect Selections by
# criterion per model (Before or After Box Cox transform)
ggplot(tot_preds, aes(x=Criterion, y=Total, fill=Selection)) +
  geom_bar(stat='identity', position='stack', width=.75) +
  scale_fill_manual('Selection', 
                    values = c('firebrick3', 'forestgreen'), 
                    labels = c('Incorrect', 'Correct')) +
  geom_hline(aes(yintercept=length(tsm_df$Position))) + 
  scale_y_continuous(breaks = sort(c(seq(0,60,by=10),length(tsm_df$Position)))) +
  geom_text(aes('BIC',length(tsm_df$Position),
                label = paste('Ground True Positions =',length(tsm_df$Position)), 
                vjust = -1), size=2.5) +
  labs(title='Total Correct and Incorrect Significant Position Selections',
       subtitle='Before and After Box Cox Transformation',
       x='Criterion',
       y='Number of Selected Significant Positions') +
  facet_wrap(~Model) +
  theme_light()

## ------------------------------------------------------------------------
# Create a table format for the total results
sketch = htmltools::withTags(table(
  class = 'cell-border',
  thead(
    tr(
      th(rowspan=2, 'Drug'),
      lapply(unique(tot_preds$Criterion), function(x) th(colspan=3, x))
    ),
    tr(
      lapply(rep(c('Correct', 'Incorrect', 'Correct/Total Selected'), 2), th)
    )
  )
))

# Show the total results as a datatable
datatable(tot_res_dt, container = sketch, 
          options=list(scrollX=TRUE), rownames=tot_res_dt$Model,
          caption='Comparing the Total Covariate Positions Selected by Both Criterion Before and After Box Cox Transformation')


## ---- fig.width=8, fig.height=11, warning=FALSE, error=FALSE-------------
# Generate a stacked bar plot to visualize the correct and incorrect Selections by
# criterion per drug from before the Box Cox Transform
g0 <- ggplot(init_regsub_res_df, aes(x=Criterion, y=Total, fill=Selection)) +
  geom_bar(stat='identity', position='stack', width=.75) +
  scale_fill_manual('Selection', 
                    values = c('firebrick3', 'forestgreen'), 
                    labels = c('Incorrect', 'Correct')) +
  ylim(0, 50) +
  labs(title='Initial\nSelection',
       x='',
       y='Number of Selected Significant Positions') +
  facet_wrap(~Drug, ncol=1) +
  theme_light()

# Generate a stacked bar plot to visualize the correct and incorrect Selections by
# criterion per drug from after the Box Cox Transform
g2 <- ggplot(bc_regsub_res_df, aes(x=Criterion, y=Total, fill=Selection)) +
  geom_bar(stat='identity', position='stack', width=.75) +
  scale_fill_manual('Selection', 
                    values = c('firebrick3', 'forestgreen'), 
                    labels = c('Incorrect', 'Correct')) +
  ylim(0, 40) +
  labs(title='Selection After Box Cox\nTransformation',
       x='',
       y='') +
  facet_wrap(~Drug, ncol=1) +
  theme_light()

# Show the two plots created above side-by-side
plot_grid(g0, g2, nrow=1) + 
  draw_figure_label('Criterion', position = 'bottom')


## ---- warning=FALSE, error=FALSE, message=FALSE--------------------------
# Create a table format for the results per drug
sketch = htmltools::withTags(table(
  class = 'cell-border',
  thead(
    tr(
      th(rowspan=2, 'Drug'),
      lapply(unique(init_regsub_res_df$Criterion), function(x) th(colspan=3, x))
    ),
    tr(
      lapply(rep(c('Correct', 'Incorrect', 'Correct/Total Selected'), 2), th)
    )
  )
))

# Create a datatable for the results per drug prior to Box Cox transform
init_dt <- datatable(init_regsub_res_dt[,2:ncol(init_regsub_res_dt)], container = sketch, 
                     options=list(scrollX=TRUE), rownames=init_regsub_res_dt$Drug,
                     caption='Comparing the Covariate Positions Chosen Prior to Box Cox Transform')

# Create a datatable for the results per drug after Box Cox transform
bc_dt <- datatable(bc_regsub_res_dt[,2:ncol(bc_regsub_res_dt)], container = sketch, 
                   options=list(scrollX=TRUE), rownames=bc_regsub_res_dt$Drug,
                   caption='Comparing the Covariate Positions Chosen After the Box Cox Transform')


## ------------------------------------------------------------------------
init_dt


## ------------------------------------------------------------------------
bc_dt


## ---- warning=FALSE, error=FALSE, message=FALSE--------------------------
# Initialize empty lists to hold
#   - bc_regsub_bf: The result of using regsub on the Box Cox transformed
#                   linear regression of resistance response vector against 
#                   the design matrix for each drug
#   - bf_sel_raw: The position-mutation pair p-values vectors for each
#                 linear regression in bc_regsub_bf
#   - bf_sel_pos: Vectors of the positions deemed significant by p-value
#                 after the Bonferroni correction was performed
bc_regsub_bf = list()
bf_sel_raw = list()
bf_sel_pos = list()

# Set significance level alpha to 0.05
alpha = 0.05

# Loop through all drugs in response matrix Y
for(drug in names(Y)){
  
  # Deal with the response vector corresponding to the current drug
  y = Y[[drug]]
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  x = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  x = x[,colSums(x) >= 3]
  # Remove duplicate predictors.
  x = x[,colSums(abs(cor(x)-1) < 1e-4) == 1]
  
  # If the optimal lambda is 0, then generate the linear regression of the natural log of the resistance response 
  # vector against the design matrix for each drug. Otherwise, generate the linear regression of the 
  # resistance response vector to the power of lambda against the design matrix for each drug.
  if(bac.lambda[[drug]] == 0){
    bc_regsub_bf[[drug]] = summary(lm(log(y)~., data=data.frame(x)))
  } else {
    bc_regsub_bf[[drug]] = summary(lm((y**bac.lambda[[drug]])~., data=data.frame(x)))
  }
  
  # Perform the Bonferroni correction on the p-values generated by the linear regression for the coefficients
  bf_p = p.adjust(bc_regsub_bf[[drug]]$coefficients[,4], 'bonferroni')
  # Remove p-values greater than alpha from the rsulting vector
  bf_sel_raw[[drug]] = bf_p[bf_p <= alpha]
  # From the coefficient vector extracted above, remove the mutation from the position mutation pair, make all
  # the resulting positions numeric values, and finally only keep the unique values of the resulting vector 
  bf_sel_pos[[drug]] <- names(bf_sel_raw[[drug]]) %>%
    substring(2,3) %>%
    as.numeric() %>%
    unique() %>%
    sort()
  # Remove NA values from the vector of positions
  bf_sel_pos[[drug]] = bf_sel_pos[[drug]][!is.na(bf_sel_pos[[drug]])]
}

bf_sel_pos


## ------------------------------------------------------------------------
# Combine all unique selected positions into one vector
all_bf_pos = unlist(bf_sel_pos) %>%
  unname() %>%
  unique()

# Check which positions deemed significant are in the list of ground true positions
bf_chk = all_bf_pos %in% tsm_df$Position
bf_correct = sum(bf_chk)
bf_incorrect = length(bf_chk) - bf_correct

bf_df = tibble('Method'= c('Bonferroni Correction', 'Bonferroni Correction'),
               'Selection'=c('Correct', 'Incorrect'),
               'Total'=c(bf_correct, bf_incorrect)) %>%
  mutate(Selection=factor(Selection, levels=c('Incorrect', 'Correct')))

# Generate a stacked bar plot to visualize the correct and incorrect selections
ggplot(bf_df, aes(x=Method, y=Total, fill=Selection)) +
  geom_bar(stat='identity', position='stack', width=.75) +
  scale_fill_manual('Selection', 
                    values = c('firebrick3', 'forestgreen'), 
                    labels = c('Incorrect', 'Correct')) +
  geom_hline(aes(yintercept=length(tsm_df$Position))) + 
  scale_y_continuous(breaks = sort(c(seq(0,60,by=10),length(tsm_df$Position)))) +
  geom_text(aes('Bonferroni Correction',length(tsm_df$Position),
                label = paste('Ground True Positions =',length(tsm_df$Position)), 
                vjust = -1), size=2.5) +
  labs(title='Total Correct and Incorrect Significant Position Selections',
       x='Method',
       y='Number of Selected Significant Positions') +
  theme_light()


## ---- warning=FALSE, error=FALSE, message=FALSE--------------------------
# Initialize empty lists to hold
#   - bc_regsub_bh: The result of using regsub on the Box Cox transformed
#                   linear regression of resistance response vector against 
#                   the design matrix for each drug
#   - bh_sel_raw: The position-mutation pair p-values vectors for each
#                 linear regression in bc_regsub_bh
#   - bh_sel_pos: Vectors of the positions deemed significant by p-value
#                 after the Benjamini-Hochberg procedure was performed
bc_regsub_bh = list()
bh_sel_raw = list()
bh_sel_pos = list()

# Set significance level alpha to 0.05
alpha = 0.05

# Loop through all drugs in response matrix Y
for(drug in names(Y)){
  
  # Deal with the response vector corresponding to the current drug
  y = Y[[drug]]
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  x = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  x = x[,colSums(x) >= 3]
  # Remove duplicate predictors.
  x = x[,colSums(abs(cor(x)-1) < 1e-4) == 1]
  
  # If the optimal lambda is 0, then generate the linear regression of the natural log of the resistance response 
  # vector against the design matrix for each drug. Otherwise, generate the linear regression of the 
  # resistance response vector to the power of lambda against the design matrix for each drug.
  if(bac.lambda[[drug]] == 0){
    bc_regsub_bh[[drug]] = summary(lm(log(y)~., data=data.frame(x)))
  } else {
    bc_regsub_bh[[drug]] = summary(lm((y**bac.lambda[[drug]])~., data=data.frame(x)))
  }
  
  # Perform the Benjamini-Hochberg procedure on the p-values generated by the linear regression for the coefficients
  bh_p = p.adjust(bc_regsub_bh[[drug]]$coefficients[,4], 'BH')
  # Remove p-values greater than alpha from the rsulting vector
  bh_sel_raw[[drug]] = bh_p[bh_p <= alpha]
  # From the coefficient vector extracted above, remove the mutation from the position mutation pair, make all 
  # the resulting positions numeric values, and finally only keep the unique values of the resulting vector 
  bh_sel_pos[[drug]] <- names(bh_sel_raw[[drug]]) %>%
    substring(2,3) %>%
    as.numeric() %>%
    unique() %>%
    sort()
  # Remove NA values from the vector of positions
  bh_sel_pos[[drug]] = bh_sel_pos[[drug]][!is.na(bh_sel_pos[[drug]])]
}

bh_sel_pos


## ------------------------------------------------------------------------
# Combine all unique selected positions into one vector
all_bh_pos = unlist(bh_sel_pos) %>%
  unname() %>%
  unique()

# Check which positions deemed significant are in the list of ground true positions
bh_chk = all_bh_pos %in% tsm_df$Position
bh_correct = sum(bh_chk)
bh_incorrect = length(bh_chk) - bh_correct

bh_df = tibble('Method'= c('Benjamini-Hochberg Procedure', 'Benjamini-Hochberg Procedure'),
               'Selection'=c('Correct', 'Incorrect'),
               'Total'=c(bh_correct, bh_incorrect)) %>%
  mutate(Selection=factor(Selection, levels=c('Incorrect', 'Correct')))

# Generate a stacked bar plot to visualize the correct and incorrect selections
ggplot(bh_df, aes(x=Method, y=Total, fill=Selection)) +
  geom_bar(stat='identity', position='stack', width=.75) +
  scale_fill_manual('Selection', 
                    values = c('firebrick3', 'forestgreen'), 
                    labels = c('Incorrect', 'Correct')) +
  geom_hline(aes(yintercept=length(tsm_df$Position))) + 
  scale_y_continuous(breaks = sort(c(seq(0,60,by=10),length(tsm_df$Position)))) +
  geom_text(aes('Benjamini-Hochberg Procedure',length(tsm_df$Position),
                label = paste('Ground True Positions =',length(tsm_df$Position)), 
                vjust = -1), size=2.5) +
  labs(title='Total Correct and Incorrect Significant Position Selections',
       x='Method',
       y='Number of Selected Significant Positions') +
  theme_light()

