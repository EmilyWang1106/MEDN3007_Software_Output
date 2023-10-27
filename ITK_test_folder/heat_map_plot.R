#This file has been created to plot the heatmap of ITK under ClinVar and PPI-Hotspot Models 


#install.packages("viridis")
library(viridis)

# Amino acids ordered by hydrophobicity
hydrophobic_order <- c("I", "V", "L", "F", "C", "M", "A", "G", "T", "S", "W", "Y", "P", "H", "E", "Q", "D", "N", "K", "R")

itk.clinvar.prediction <- itk.clinvar.prediction %>%
  mutate(aa_alt = factor(aa_alt, levels = hydrophobic_order))


itk_clinvar_plot <- ggplot(itk.clinvar.prediction, aes(x = aa_pos, y = aa_alt, fill = .pred_pathogenic))+
  geom_tile()+
  scale_fill_viridis_c() +
  labs(title = "", x = "Position", y ="Amino Acids", fill = "Pathogenicity")+
  theme_minimal()

print(itk_clinvar_plot)

itk.hotspot.prediction <- itk.hotspot.prediction %>%
  mutate(aa_alt = factor(aa_alt, levels = hydrophobic_order))


itk_hotspot_plot <- ggplot(itk.hotspot.prediction, aes(x = aa_pos, y = aa_alt, fill = .pred_PPI))+
  geom_tile()+
  scale_fill_viridis_c() +
  labs(title = "", x = "Position", y ="Amino Acids", fill = "PPI Disruption")+
  theme_minimal()

print(itk_hotspot_plot)