
source("R/Functions.R")

Index <- IndexCal("Folder/subFolder/index.html")%>%data.table()

PSM <- fread("Folder/subFolder/PD_Export_FileName_PSMs.txt")

source("R/Process_HCD.R")

#write_tsv(Lib,"Test2.tsv")
#write_tsv(Lib[rank == 1 & XCorr >= 1.2,.(ProteinGroups,RelativeIntensity,FragmentMz,ModifiedPeptide,IntModifiedPeptide,StrippedPeptide,PrecursorCharge,PrecursorMz,iRT,FragmentNumber,FragmentType,FragmentCharge,FragmentLossType)],"Library_Name.tsv")

write_tsv(Lib[rank == 1 & XCorr >= 0.0,.(ProteinGroups,RelativeIntensity,FragmentMz,ModifiedPeptide,IntModifiedPeptide,StrippedPeptide,PrecursorCharge,PrecursorMz,iRT,FragmentNumber,FragmentType,FragmentCharge,FragmentLossType)],"Library_Name.tsv")


