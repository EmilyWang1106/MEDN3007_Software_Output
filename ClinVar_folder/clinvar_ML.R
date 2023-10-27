#This file has been created to build the ClinVar model

library(tidymodels)
library(readr)
library(dplyr)
library(parallel)
library(ranger)
master_joined <- read_csv("ClinVar_folder/master_joined.csv") # Load the master_joined table for 20 proteins. 

############

master_joined_clinvar <- master_joined %>%
  mutate(across(c(aa_ref, aa_alt, sec_structure), as.factor)) %>%
  select(-spotone)|>
  filter(aa_ref != aa_alt) %>% #Remove non-mutations
  filter(clinvar != 2) |> #Remove all rows with no information of pathogenity 
  drop_na(everything()) |> # Drop NAs 
  mutate(clinvar_cat = as.factor(if_else(clinvar==0, "benign", "pathogenic"))) # Creat a new coloumn for label
  



clinvar_split <- initial_split(master_joined_clinvar, prop = 0.8) # Split in to traning set and testing set 
clinvar_training = training(clinvar_split)
clinvar_testing = testing(clinvar_split)

clinvar_recipe <- clinvar_training %>% # Build the recipe will be used in the model
  recipe(clinvar_cat ~ .) %>%
  update_role(hotspot, new_role = "other_trainer") %>% # Ignore columns showing id
  update_role(protein_A, new_role = "id") %>%
  update_role(uniprot, new_role = "id") %>%
  update_role(clinvar, new_role = "id") |>
  step_corr(all_numeric_predictors()) %>% 
  step_center(all_numeric_predictors(), -all_outcomes()) %>%
  step_scale(all_numeric_predictors(), -all_outcomes()) %>%
  # step_naomit(everything()) %>%
  prep()

clinvar_model = rand_forest(trees = 100) %>% # initiate the model by random forest, 100 trees have been used. 
  set_engine("ranger", 
             write.forest = T,
             verbose = T, 
             keep.inbag = T, 
             num.threads = parallel::detectCores(),
             seed = 1) %>% 
  set_mode("classification")

clinvar_wf = workflow() %>%  # Build the work flow
  add_model(clinvar_model) %>%
  add_recipe(clinvar_recipe)

clinvar_fit = clinvar_wf %>% # build the model based on training set 
  fit(data = clinvar_training)

clinvar_predicted = predict(clinvar_fit, clinvar_testing, type="prob") %>% # run the model on validation set 
  bind_cols(clinvar_testing)

clinvar_predicted %>% # Make the roc curve 
  roc_curve(clinvar_cat, .pred_benign) %>%
  autoplot()

clinvar_predicted %>% # Make the gain curve 
  gain_curve(clinvar_cat, .pred_benign) %>%
  autoplot()

#######
#Do the cross validation 
set.seed(1) # Ensure data splitting in cross-validation are reproducible

data_folds <- vfold_cv(master_joined_clinvar, strata = "clinvar_cat", v = 5) # Do the 5-fold cross-validation 

clinvar_fit_rs <- clinvar_wf %>% # Fit the model to the cross-validations
  fit_resamples(data_folds, control = control_resamples(save_pred = TRUE))


clinvar_predictions_rs <- clinvar_fit_rs %>% # Collect the predictions for cross-validations
  collect_predictions() 

clinvar_predictions_rs %>% # Make the roc curve based on the cross-validations
  roc_curve(clinvar_cat, .pred_benign) %>%
  autoplot()

clinvar_predictions_rs %>%  # Make the gain curve based on the cross-validations
  gain_curve(clinvar_cat, .pred_benign) %>%
  autoplot()

collect_metrics(clinvar_fit_rs) # collect the averaged accuracy and roc-auc score for the cross-validations  

conf_mat(clinvar_predictions_rs, truth = clinvar_cat, estimate = .pred_class) # Make the confustion map for the cross-validations

