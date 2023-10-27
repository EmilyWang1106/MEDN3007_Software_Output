library(data.table)

AA_alt <-DDG_table$Amino.Acids
AA_ref <- itk_ref
Master_table<- data_frame(DDG_table$X.Position, AA_ref, AA_alt,DDG_table$DDG)
Master_table <- as.data.table(Master_table)

colnames(Master_table)[1] <-"Protein_Position"
colnames(Master_table)[4] <-"DDG"

pp_table <- fread("/Users/yiyangwang/Desktop/MEDN3007/ITK_VEP/ITK_polyphen.txt")
Master_table[pp_table, Polyphen:=pph2_prob, on = c(AA_ref="aa1",AA_alt="aa2", Protein_Position = "pos")]

Master_table[is.na(Master_table)] <- 0

Master_table[SA_table, ASA:= asa, on = c(AA_ref="seq", Protein_Position = "n")]
Master_table[SA_table, RSA:= rsa, on = c(AA_ref="seq", Protein_Position = "n")]

Master_table[CSM_table, CSM_Potential_Score:= Score, on = c(Protein_Position = "Position")]


