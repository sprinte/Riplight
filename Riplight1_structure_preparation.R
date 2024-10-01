### CLEAN UP SEQUENCING DATA ###

#### USER PART ####
library(tidyr)
library(dplyr)

path_in<-"F:/Documents/Studium/Promotion/EcoMigRip/RipLight1/Data/Structure & Abundance/Sequencing Raw Data"
path_out<-path_in

parameter<-"Bacteria"# either"Bacteria" or "Fungi"
project<-"RipLight1" # to name output files

#### END OF USER PART

#### DATA PREPARATION FROM PIPELINE ####
if (tolower(parameter)=="bacteria") {
  path_fragment<-"16S"
} else if (tolower(parameter)=="fungi") {
  path_fragment<-"ITS"
}

# directory
complete_path_in<-paste(path_in, path_fragment, "output", sep="/")
complete_path_out<-paste(path_out, path_fragment, "output", sep="/")

setwd(complete_path_in)

tax<-read.table("ASVs_taxonomy.tsv", header=T, sep="\t", stringsAsFactors = F)
taxcl<-tax

# eliminate NA rows
na_ASVs<-taxcl[rowSums(is.na(taxcl[-1])) == ncol(taxcl[-1]),][[1]] # ASVs that have all NAs
taxcl<-taxcl[rowSums(is.na(taxcl[-1])) != ncol(taxcl[-1]),] # counts without NA-ASVs
rownames(taxcl)<-1:nrow(taxcl)

# eliminate species column if all empty
if (all(is.na(taxcl$species))==TRUE) {
  colind<-which(colnames(taxcl) %in% "species")
  taxcl<-taxcl[, -colind]
}
names(taxcl)[1]<-"ASV_ID"

taxcl<-taxcl %>%
  replace(is.na(.), "unclassified")

#### CLEANING FOR BACTERIA ####
if (parameter=="Bacteria") {
for (i in 1:nrow(taxcl)) { # loop through rows
    taxcl[i, "genus"]<-gsub("\\([0-9]*\\)$", "", taxcl[i, "genus"]) # remove quality index (numbers in brackets)
    taxcl[i, "genus"]<-gsub(".*uncultured.*", "unclassified", taxcl[i, "genus"], ignore.case = T) # change "uncultured" to "unclassified"
    taxcl[i, "genus"]<-gsub(".*Candidatus.*", "unclassified", taxcl[i, "genus"], ignore.case = T)
    taxcl[i, "genus"]<-gsub(".*(_cl|_or|_fa|_ge)$", "unclassified", taxcl[i, "genus"])
    #taxcl[i, "genus"]<-gsub("_\\(Subgroup_[0-9]*\\).*", "", taxcl[i, "genus"], ignore.case = T) # either consider Subgroup or treat as "unclassified"
    taxcl[i, "genus"]<-gsub("_(marine|terrestrial|soil|sediment|termite|wastewater-sludge)_group.*", "", taxcl[i, "genus"], ignore.case = T)
  for (j in rev(2:ncol(taxcl[-1]))) { # for each row, go through columns from right (genus) to left (phylum)
      taxcl[i, j]<-gsub("\\([0-9]*\\)$", "", taxcl[i, j])
      taxcl[i, j]<-gsub(".*uncultured.*", "unclassified", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*Candidatus.*", "unclassified", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*(_cl|_or|_fa|_ge)$", "unclassified", taxcl[i, j])
      taxcl[i, j]<-gsub("_\\(Subgroup_[0-9]\\).*", "", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub("\\(.*group.*\\)", "", taxcl[i, j], ignore.case = T)
      #taxcl[i, j]<-gsub("^Subgroup_[0-9]*$", "unclassified", taxcl[i, j], ignore.case = T) # either consider Subgroup or treat as "unclassified"
      taxcl[i, j]<-gsub("_(marine|terrestrial|soil|sediment|termite|wastewater-sludge)_group.*", "", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*Unknown.*", "unclassified", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*Incertae_sedis.*", "unclassified", taxcl[i, j], ignore.case = T)
    if (taxcl[i, j]=="unclassified") {
      for (k in 1:(ncol(taxcl)-j)) {
        if(taxcl[i, j+k]=="unclassified") {
          next
        } else if (taxcl[i, j+k]!="unclassified") {
          taxcl[i, j]<-paste("no", names(taxcl[j]), sep="_")
        }
      }
    }
  }
}
  taxcl<-taxcl[taxcl$domain=="Bacteria"|taxcl$domain=="Archaea",]
  
  # remove Chloroplasts
  row_sub<-apply(taxcl, 1, function(row) all(!(row %in% c("Chloroplast", "chloroplast"))))
  taxcl<-taxcl[row_sub, ]
  
  #taxcl$ASV_ID<-gsub("p_", "", taxcl$ASV_ID)
}

#### CLEANING FOR FUNGI ####
if (parameter=="Fungi") {
  for (i in seq_along(1:nrow(taxcl))) { # loop through rows
    taxcl[i, "species"]<-gsub("^(s__|g__|f__|o__|c__|p__|k__)|^(s_|g_|f_|o_|c_|p_|k_)", "", taxcl[i, "species"]) # remove quality index (numbers in brackets)
    taxcl[i, "species"]<-gsub("\\([0-9]*\\)$", "", taxcl[i, "species"]) # remove quality index (numbers in brackets)
    taxcl[i, "species"]<-gsub(".*uncultured.*", "unclassified", taxcl[i, "species"], ignore.case = T) # change "uncultured" to "unclassified"
    taxcl[i, "species"]<-gsub(".*Candidatus.*", "unclassified", taxcl[i, "species"], ignore.case = T)
    taxcl[i, "species"]<-gsub(".*(_cl|_or|_fa|_ge|_sp)$", "unclassified", taxcl[i, "species"])
    taxcl[i, "species"]<-gsub("\\(.*group.*\\)", "", taxcl[i, "species"], ignore.case = T)
    taxcl[i, "species"]<-gsub("_\\(Subgroup_[0-9]*\\).*", "", taxcl[i, "species"], ignore.case = T) # either consider Subgroup or treat as "unclassified"
    taxcl[i, "species"]<-gsub("_(marine|terrestrial|soil|sediment|termite|wastewater-sludge)_group.*", "", taxcl[i, "species"], ignore.case = T)
    for (j in rev(2:ncol(taxcl[-1]))) { # for each row, go through columns from right (species) to left (phylum)
      taxcl[i, j]<-gsub("^(s__|g__|f__|o__|c__|p__|k__)|^(s_|g_|f_|o_|c_|p_|k_)", "", taxcl[i, j]) # remove quality index (numbers in brackets)
      taxcl[i, j]<-gsub("\\([0-9]*\\)$", "", taxcl[i, j])
      taxcl[i, j]<-gsub(".*uncultured.*", "unclassified", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*Candidatus.*", "unclassified", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*(_cl|_or|_fa|_ge|_sp)$", "unclassified", taxcl[i, j])
      taxcl[i, j]<-gsub("\\(.*group.*\\)", "", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub("_\\(Subgroup_[0-9]*\\)*.*", "", taxcl[i, j], ignore.case = T) 
      #taxcl[i, j]<-gsub("^Subgroup_[0-9]*$", "unclassified", taxcl[i, j], ignore.case = T) # either consider Subgroup or treat as "unclassified"
      taxcl[i, j]<-gsub("_(marine|terrestrial|soil|sediment|termite|wastewater-sludge)_group.*", "", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*Unknown.*", "unclassified", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*Incertae_sedis.*", "unclassified", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*unidentified.*", "unclassified", taxcl[i, j], ignore.case = T)
      taxcl[i, j]<-gsub(".*unclassified.*", "unclassified", taxcl[i, j], ignore.case = T)
      if (taxcl[i, j]=="unclassified") {
        for (k in 1:(ncol(taxcl)-j)) {
          if(taxcl[i, j+k]=="unclassified") {
            next
          } else if (taxcl[i, j+k]!="unclassified") {
            taxcl[i, j]<-paste("no", names(taxcl[j]), sep="_")
          }
        }
      }
    }
  }
  taxcl<-taxcl[taxcl$kingdom=="Fungi",]
  #taxcl$ASV_ID<-gsub("e_", "", taxcl$ASV_ID)
}

#### WRITE OUTPUT FILES ####
setwd(complete_path_out)
write.table(taxcl, file=paste("Taxonomy_cleaned_", parameter, "_", project, ".txt", sep=""), sep="\t", row.names=FALSE, quote = F)