setwd("E:/Project_GitHub/fetalGeneExpression")

refNCBI <- read.delim("NCBI_Ref.bed")
refHGNC <- read.delim("HGNC.bed") 

refNCBI$name <- sub("\\.\\d+","",refNCBI$name)
refNCBI <- refNCBI[!grepl("_",refNCBI$chrom, perl = T),]
mergeGenes <- merge(refHGNC, refNCBI, by.x = "Approved.symbol", by.y = "name2")
mergeGenes <- mergeGenes[,c(-2,-4,-5,-7,-11,-22:-24)]
#Write reliable genes in HGNC; 
write.table(mergeGenes, "HGNC_reliable.bed", sep = "\t", row.names = F, quote = F);

#write.table(TSES, "HGNC_reliable.bed", sep = "\t", row.names = F, quote = F);
TSES <- unique(mergeGenes[,c(1,8:11)])
for (i in 1:nrow(TSES)) {
    if (TSES$strand[i] == "-") {
        TSES$txStart[i] <- TSES$txEnd[i];
    } 
}

TSES <- unique(TSES[,-5]);
TSES$RangeB <- TSES$txStart - 2500;
TSES$RangeE <- TSES$txStart + 2500;
TSES$RangeB[TSES$RangeB < 0] <- 1;
TSES$id <- paste(TSES$Approved.symbol,TSES$chrom,TSES$txStart, sep = "_")
TSES <- TSES[,c(-1,-3,-4)]
write.table(TSES, "TSS.bed", row.names = F, col.names = F, sep = "\t", quote = F)

TSES_P_E$RangeB <- TSES_P_E$cdsEnd - 10000;
TSES_P_E$RangeE <- TSES_P_E$cdsEnd + 10000;
TSES_P_E$RangeB[TSES_P_E$RangeB < 0] <- 1;
TSES_P_E$id <- paste(TSES_P_E$name2,TSES_P_E$chrom,TSES_P_E$cdsEnd, sep = "_")
write.table(TSES_P_E[,c(3,24,25,26,5)], "TSE_protein.bed", row.names = F, col.names = F, sep = "\t", quote = F)

#no coding gene
TSES_N_S <- TSES_N_E <- TSES[grepl("^NR",TSES$name, perl = T),-1]

for (i in 1:nrow(TSES_N_S)) {
    if (TSES_N_S$strand[i] == "-") {
        TSES_N_S$cdsStart[i] <- TSES_N_S$cdsEnd[i]; 
        TSES_N_E$cdsEnd[i] <- TSES_N_E$cdsStart[i];
    }
}
