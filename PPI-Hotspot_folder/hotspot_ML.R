#This file has been created to build the ClinVar model

library(tidymodels)
library(readr)
library(dplyr)
library(parallel)
library(ranger)

############
hotspot_joined <-  read_csv("PPI-Hotspot_folder/hotspot_joined.csv") #Load the master table for 20 proteins to run hotspot model 

hotspot_joined_init <- hotspot_joined %>%
  mutate(across(c(aa_ref, aa_alt, sec_structure, label), as.factor)) %>%
  select(-spotone)|>
  filter(aa_ref != aa_alt) %>% #Remove non-mutations 
  filter(label != 'None') |>  #Remove columns with no information in the label column (not bengin in ClinVar and not PPI in PPI Hotspot database)
  drop_na(everything()) |> # Remove NAs
  mutate(label = droplevels(label))




hotspot_split <- initial_split(hotspot_joined_init, prop = 0.8) # Split into training and testing sets
hotspot_training = training(hotspot_split)
hotspot_testing = testing(hotspot_split)

hotspot_recipe <- hotspot_training %>% #Make the recipe for the machine learning model 
  recipe(label ~ .) %>%
  update_role(clinvar, new_role = "id") %>% # Drop columns providing id information 
  update_role(protein_A, new_role = "id") %>%
  update_role(uniprot, new_role = "id") %>%
  update_role(hotspot, new_role = "id") |>
  step_corr(all_numeric_predictors()) %>% 
  step_center(all_numeric_predictors(), -all_outcomes()) %>%
  step_scale(all_numeric_predictors(), -all_outcomes()) %>%
  # step_naomit(everything()) %>%
  prep()

hotspot_model = rand_forest(trees = 100) %>% #Initiate the random forest model, 100 trees were used 
  set_engine("ranger", 
             write.forest = T,
             verbose = T, 
             keep.inbag = T, 
             num.threads = parallel::detectCores(),
             seed = 1) %>% 
  set_mode("classification")

hotspot_wf = workflow() %>% #Build the workflow 
  add_model(hotspot_model) %>%
  add_recipe(hotspot_recipe)

hotspot_fit = hotspot_wf %>% #Build the model by fiting the workflow on the traing set 
  fit(data = hotspot_training)

hotspot_predicted = predict(hotspot_fit, hotspot_testing, type="prob") %>% # Apply the model to the validation set 
  bind_cols(hotspot_testing) 

hotspot_predicted %>% # Plot the roc curve 
  roc_curve(label, .pred_benign) %>%
  autoplot()

#######
#Make the cross-validations 

set.seed(1) #Ensure data spliting for cross-validation is reproducible 
data_folds <- vfold_cv(hotspot_joined_init, strata = "label", v = 5) #Make 5-fold cross-validation 

hotspot_fit_rs <- hotspot_wf %>% #Fit the workflow to corss-validations 
  fit_resamples(data_folds, control = control_resamples(save_pred = TRUE))


hotspot_predictions_rs <- hotspot_fit_rs %>% #Collect the predictions for the cross-validations 
  collect_predictions() 

hotspot_predictions_rs %>% #Plot the roc curve for cross-validations
  roc_curve(label, .pred_benign) %>%
  autoplot()

hotspot_predictions_rs %>%  #Plot the gain curve for cross-validations
  gain_curve(label, .pred_benign) %>%
  autoplot()

collect_metrics(hotspot_fit_rs)# Collect the averaged accuracy and roc-anu scores for cross validations

conf_mat(hotspot_predictions_rs, truth = label, estimate = .pred_class)#Make the confuion matrix. 

