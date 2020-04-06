# RuBisCoEvo
In this reppository you can find several files of the Bachelor Thesis project, the folders are organised by logical analysis procedure:

## Preliminar data
In this folder you will find all raw data from the databases and scripts to organise it. The `db.csv` is the database used to organise all taxonomical and sequences data just like these sequences, the `seq_rename.R` is a R script that allows to rename all sequences with a local database naming based on `>taxID_accession_no_coded_protein` fasta sequence name. In this folder organisation
```
|- preliminar
   |- seqs
      |- aa
      |- nt
   |- names
```
you will find the raw downloaded sequences in `seqs` folder, as well as the scripts to rename the sequences. In the root of this folder we held the `kegg_retrieval.R` script used to download KEGG sequences. `names` folder contains tables used to relate sequences names with taxonomical, useful data to develop the trees observation.

## Cured data
In this folder you will find all the generated data, organised in several folders for each step of the analysis and type of data used.

## Scripts
This contains the scripts used to develop all outputs you will see at cured data.

## Supplementary data
As provided in text, there are som tables that are supplemental to the maintext of the work. In this folder, these tables are stored in `.csv` file format.
