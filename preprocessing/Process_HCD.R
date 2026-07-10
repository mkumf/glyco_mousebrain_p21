I2 <- copy(Index)

Index <- copy(I2)

colnames(Index) <- str_remove_all(colnames(Index),"Â")
colnames(Index) <- str_replace_all(colnames(Index)," ", " ")

Index[,Sequence := str_split_fixed(Sequence,"\\.",3)[,2]]%>%
  .[,Sequence := str_replace_all(Sequence, "m", "X")%>%toupper()]

Index[,`PSMs Peptide ID` := str_split_fixed(File,"_",4)[,4]]%>%
  .[,`PSMs Peptide ID` := str_remove(`PSMs Peptide ID`,".txt")]

Index[,Identifier := str_c(Sequence, `PSMs Peptide ID`, sep = "_")]

PSM[,Sequence := str_split_fixed(`Annotated Sequence`,"\\.",3)[,2]]%>%
  .[,Sequence := str_replace_all(Sequence, "m", "X")%>%toupper()]

PSM[,Identifier := str_c(Sequence, `PSMs Peptide ID`, sep = "_")]

Index <- merge(Index,
               PSM[,.(Identifier,`Spectrum File`)],by = "Identifier")

#Change path here

Index[,File := str_c("Folder/subFolder/",File)]

FragList <- data.table()

FragList <- adply(Index$File,1,SubFrag)%>%data.table()

FragList[,rank := rank(-i), by = File]

FragList <- FragList[rank <= 15 | matches != ""]%>%
  .[, matches := str_split_fixed(matches,",",2)[,1]]%>%
  .[, matches := str_replace(matches, "]"," (0)")]

FragList <- FragList[!(matches %like% "\\+H")]

FragList <- FragList[!(matches %like% "M" & matches %like% "NH3|H2O")]

FragList[,matches := str_replace(matches,"-"," ")]%>%
  .[,FragmentType := str_split_fixed(matches," ", 4)[,1]]%>%
  .[FragmentType %like% "M",FragmentType := "y"]%>%
  .[,FragmentNumber := str_split_fixed(matches," ", 4)[,2]]%>%
  .[,FragmentNumber := str_remove_all(FragmentNumber,"[\\(\\)\\+]")%>%as.integer()]%>%
  .[,FragmentCharge := str_split_fixed(matches," ", 4)[,3]]%>%
  .[,FragmentCharge := str_remove_all(FragmentCharge,"[\\(\\)\\+]")%>%as.integer()]%>%
  .[,FragmentLossType := str_split_fixed(matches," ", 4)[,4]]

Lib <- merge(Index,FragList,by = "File")

Lib[,StrippedPeptide := Sequence]%>%
  .[FragmentNumber == 0, FragmentNumber := str_length(StrippedPeptide)]%>%
  .[,RelativeIntensity := max(i), by = File]%>%
  .[,RelativeIntensity := i/RelativeIntensity]%>%
  .[,FragmentMz := m.z]%>%
  .[,PrecursorCharge := Charge]%>%
  .[,PrecursorMz := `Precursor m/z [Da]`]%>%
  .[,iRT := `RT [min]`]%>%
  .[FragmentLossType == "",FragmentLossType := "noloss"]%>%
  .[,ProteinGroups := `Protein Accessions`]%>%
  .[,Modification := str_split_fixed(`Spectrum File`, "G\\)_",2)[,2]]%>%
  .[,Modification := str_remove_all(Modification, ".mgf")]

Lib <- Lib[FragmentLossType != "H" & !(FragmentLossType %like% " O")]

Lib[!(Modification %like% "_") & Modification %like% "Tn", Modification := str_c(Modification, "0xT", sep = "_")]%>%
  .[!(Modification %like% "_"), Modification := str_c("0xTn", Modification, sep = "_")]

Lib[,Modification := str_remove_all(Modification,"[xTn]")]


Glycans <- data.table(Modification = c("0_1","0_2","0_3","0_4","0_5","1_0","1_1","1_2","1_3","1_4","2_0","2_1","2_2","2_3","3_0","3_1","3_2","4_0","4_1","5_0"))

Glycans[,Tn := str_split_fixed(Modification, "_" ,2)[,1]%>%as.integer()]%>%
  .[,`T` := str_split_fixed(Modification, "_" ,2)[,2]%>%as.integer()]%>%
  .[, Hex := `T`]%>%
  .[, HexNAc := `T` + Tn]

Glycans[,Synonyms := str_c("[Hex(",Hex,")HexNAc(",HexNAc,") (Any N-term)]")]%>%
  .[,Synonyms := str_remove_all(Synonyms,"Hex\\(0\\)")]%>%
  .[,Synonyms := str_remove_all(Synonyms,"HexNAc\\(0\\)")]%>%
  .[Synonyms == "[HexNAc(1) (Any N-term)]",Synonyms := "[HexNAc (Any N-term)]"]

Glycans[,C := 8 * `Tn` + 14 * `T`]%>%
  .[,H  := 13 * `Tn` + 23 * `T`]%>%
  .[,N  := 1 * `Tn` + 1 * `T`]%>%
  .[,O  := 5 * `Tn` + 10 * `T`]

# Glycans[,Composition := str_c("+C",C,"+H",H,"+N",N,"+O",O)]%>%
#   .[,Composition := str_replace_all(Composition, "N1\\+","N+")]

Glycans[,Mass := 365.132196 * `T` + 203.079373 * Tn]%>%
  .[,NomiMass := str_c("[+",round(Mass),"]")]

Lib <- merge(Lib,Glycans,by = "Modification")

Lib[,ModifiedPeptide := str_c("_",Synonyms,StrippedPeptide,"_")]%>%
  .[,IntModifiedPeptide := str_c("_",NomiMass,StrippedPeptide,"_")]

Lib[FragmentType == "b" & FragmentLossType == "H2O",H := H+2]%>%
  .[FragmentType == "b" & FragmentLossType == "H2O",O := O+1]%>%
  .[FragmentType == "b" & FragmentLossType == "NH3",N := N+1]%>%
  .[FragmentType == "b" & FragmentLossType == "NH3",H := H+3]

Lib[,Composition := str_c("1(","+C",C,"+H",H,"+N",N,"+O",O,")")]%>%
  .[,Composition := str_replace_all(Composition, "N1\\+","N+")]

Lib[FragmentType == "b",FragmentLossType := Composition]

Lib[,PrecursorMz := as.numeric(PrecursorMz) + Mass/as.numeric(PrecursorCharge)]

Lib <- Lib[is.na(FragmentCharge) == FALSE,.(ProteinGroups,RelativeIntensity,FragmentMz,ModifiedPeptide,IntModifiedPeptide,StrippedPeptide,PrecursorCharge,PrecursorMz,iRT,FragmentNumber,FragmentType,FragmentCharge,FragmentLossType,XCorr)]

Lib[,ModifiedPeptide := str_replace_all(ModifiedPeptide,"C","C[Carbamidomethyl]")]%>%
  .[,ModifiedPeptide := str_replace_all(ModifiedPeptide,"X","M[Oxidation]")]%>%
  .[,IntModifiedPeptide := str_replace_all(IntModifiedPeptide,"C","C[+57]")]%>%
  .[,IntModifiedPeptide := str_replace_all(IntModifiedPeptide,"X","M[+16]")]%>%
  .[,LabeledPeptide := ModifiedPeptide]%>%
  .[,StrippedPeptide := str_replace_all(StrippedPeptide,"X","M")]

Lib[,Identifier := str_c(LabeledPeptide,PrecursorCharge,sep = "_")]%>%.[,XCorr := as.numeric(XCorr)]

Rankings <- Lib[,.(Identifier,XCorr)]%>%unique()

Rankings[,rank := rank(-XCorr),by = Identifier]

Lib <- merge(Lib, Rankings, by = c("Identifier","XCorr"))
