#The file was created to apply the ClinVar and PPI-Hotspot model on ITK

library(data.table)

master.table.itk<-fread("ITK_test_folder/Master_table_ITK.csv") #Load the itk master table 

master.table.itk <- master.table.itk %>% # Process on the itk master table to ensure it is ready for applying teh models 
  mutate(across(c(aa_ref, aa_alt, sec_structure), as.factor)) %>%
  filter(aa_ref != aa_alt) %>% 
  drop_na(everything()) 

itk.clinvar.prediction <- predict(clinvar_fit, master.table.itk, type="prob") %>% #Apply the ClinVar model to ITK 
  bind_cols(master.table.itk)

itk.hotspot.prediction <- predict(hotspot_fit, master.table.itk, type="prob") %>% #Apply the PPI-Hotspot model to ITK 
  bind_cols(master.table.itk)

clinvar_path_ppi_test <- master_joined |>  # Make the table for ClinVar pathogenic variants 
  select(-spotone)|>
  filter(aa_ref != aa_alt) %>% 
  filter(clinvar != 2) |>
  drop_na(everything()) |>
  filter(clinvar != 0)

clinvar_path_hotspot_prediction <- predict(hotspot_fit, clinvar_path_ppi_test) |> #Apply the PPI-Hotspot model to clinvar pathogenic variants 
  bind_cols(clinvar_path_ppi_test)

sum(clinvar_path_hotspot_prediction$.pred_class== "PPI") #Return the number of pathogenic variants has been identified as PPI-disruppting residues
