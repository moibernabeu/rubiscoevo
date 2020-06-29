# Moisès Bernabeu
# Started Cocentaina March 31st, 2020
# KEGG database sequences retrieval

library(KEGGREST)
library(seqinr)
library(Biostrings)

db1 <- c("loki:Lokiarch_30170", "pdl:Pyrde_1506", "apo:Arcpr_1634",
        "fpl:Ferp_1506", "mok:Metok_0719", "tpep:A0127_03925",
        "acy:Anacy_0029", "cmp:Cha6605_0645", "cthe:Chro_5313",
        "can:Cyan10605_0644")
db2 <-  c("cgc:Cyagr_0014", "cyn:Cyan7425_3422",
        "dsl:Dacsa_1767", "gen:GM3709_214", "cyc:PCC7424_1367",
        "hhg:XM38_011320", "len:LEP3755_50110", "mic:Mic7113_2336",
        "mar:MAE_47890", "nsh:GXM_06762")
db3 <-  c("pagh:NIES204_32340", "syc:syc0130_c", "tel:tll1506",
        "glj:GKIL_0669", "gvi:gvip295")

aa1 <- keggGet(db1, "aaseq")
aa2 <- keggGet(db2, "aaseq")
aa3 <- keggGet(db3, "aaseq")

aa <- c(aa1, aa2, aa3); aa
writeXStringSet(aa, format = "fasta", "Desktop/aa_kegg.fa")

nt1 <- keggGet(db1, "ntseq")
nt2 <- keggGet(db2, "ntseq")
nt3 <- keggGet(db3, "ntseq")

nt <- c(nt1, nt2, nt3); nt
writeXStringSet(aa, format = "fasta", "Desktop/nt_kegg.fa")

#----Miàu----
a <- c("tel:tll1506", "gvi:gvip295")
b <- keggGet(a, "ntseq")
writeXStringSet(b, format="fasta", "append.fa")
