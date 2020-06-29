# Moisès Bernabeu
# Started Burjassot February 24th, 2020
# Single gene multifasta renaming for local database

library(seqinr)
library(phylotools)
library(stringr)
library(readxl)

file_aa <- "aa/aa_raw_seqs.fasta"
# write.csv(get.fasta.name(file_aa), file = paste(str_split(file_aa, "\\.")[[1]][1], "_names.csv", sep = ""))

db <- read_excel("db_names.xlsx"); db

tax_aa <- read.csv("aa/names.csv", row.names = 1); View(tax_aa)

old_names_aa <- tax_aa$x
taxid_aa <- tax_aa$TaxID
gene <- "rbcL"

accession_aa <- c()
for (i in 1:length(old_names_aa)) {
  x <- str_split(old_names_aa, " ")[[i]][1]
  y <- str_split(x, "\\.")[[1]][1]
  accession_aa[i] <- y
}; head(accession_aa)

new_names_aa <- paste(taxid_aa, accession_aa, gene, sep = "_"); head(new_names_aa)

namechange_aa <- data.frame(old_names_aa, new_names_aa); head(namechange_aa)

rename.fasta(infile = file_aa, namechange_aa,
             outfile = paste(str_split(file_aa, "\\.", n = 2)[[1]][1], "_renamed.fa", sep = ""))

write.csv(data.frame(taxa = new_names_aa, Group = db$Group,
                     Species = db$Specie,  Combined = db$Combined, TaxID = db$Taxid),
          file = paste(str_split(file_aa, "\\.", n = 2)[[1]][1], "_tree_names.csv", sep = ""),
          row.names = F)

#----Nucleotide data----
file_nt <- "nt_raw.fasta"
# write.csv(get.fasta.name(file_nt), file = paste(str_split(file_nt, "\\.")[[1]][1], "_names.csv", sep = ""))

tax_nt <- read.csv(paste(str_split(file_nt, "\\.")[[1]][1], "_names.csv", sep = ""), row.names = 1); View(tax_nt)

old_names_nt <- tax_nt$x
taxid_nt <- tax_nt$TaxID
gene <- "rbcL"

accession_nt <- c()
for (i in 1:length(old_names_nt)) {
  a <- str_split(old_names_nt, "protein_id=", n = 2)[[i]][2]; a
  b <- str_split(a, "\\]", n = 2)[[1]][1]; b
  c <- str_split(b, "\\.", n = 2)[[1]][1]
  accession_nt[i] <- c
}; head(accession_nt)

new_names_nt <- paste(taxid_nt, accession_nt, gene, sep = "_"); head(new_names_nt)

namechange_nt <- data.frame(old_names_nt, new_names_nt); head(namechange_nt)

rename.fasta(infile = file_nt, namechange_nt,
             outfile = paste(str_split(file_nt, "\\.", n = 2)[[1]][1], "_renamed.fasta", sep = ""))

write.table(data.frame(taxa = new_names_nt, Group = tax_nt$Order, Species = tax_nt$Specie,
                     Combined = paste(tax_nt$Order, tax_nt$Specie, sep = " | ", tax_nt$Type), TaxID = tax_nt$TaxID),
          file = paste(str_split(file_nt, "\\.", n = 2)[[1]][1], "_tree_names.txt", sep = ""),
          row.names = F, sep = "\t", quote = F)
