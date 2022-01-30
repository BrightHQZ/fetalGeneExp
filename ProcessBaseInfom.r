library(readxl)

setwd("F:/sampleList")
fileList <- dir(pattern = "*.xls")
baseInfo <- ""
for (i in 1:length(fileList)) {
    if (i == 1) {
        baseInfo <- read_xls(fileList[i])
    } else {
        baseInfo <- rbind(baseInfo, read_xls(fileList[i]))
    }
}

baseInfo <- unique(baseInfo)
colnames(baseInfo) <- baseInfo[1, ]
baseInfo <- baseInfo[-1, -1]

for (i in 1:length(colnames(baseInfo))) {
    baseInfo[is.na(baseInfo[, i]), i] <- ""
}
baseInfo <- as.data.frame(baseInfo)

baseInfo$modifiedWZ <- round(as.numeric(unlist(lapply(baseInfo[, 19], function(x) {
    gsub("\\+", "\\.", x)
}))))
baseInfo[, 16] <- as.numeric(baseInfo[, 16])
baseInfo[, 18] <- as.numeric(baseInfo[, 18])
baseInfo[, 22] <- as.numeric(baseInfo[, 22])
baseInfo[, 23] <- as.numeric(baseInfo[, 23])
baseInfo[, 24] <- as.numeric(baseInfo[, 24])

write.table(baseInfo, "combinedSample.txt", row.names = F, quote = F, sep = "\t")

caseDir <- read_xls("F:/sampleList/sampleDir/sampleDir.xls")

combinedRes <- merge(baseInfo, caseDir, by.x = "实验编号", by.y = "实验编号...3", all.x = T)

boys <- combinedRes[combinedRes$fetal2 != "NULL" & combinedRes$fetal2 != "", ]
boys <- boys[!is.na(boys[, 1]), ]
boys$YYS <- unlist(apply(boys,1,function(x) { return (x[4] == x[85])} ))
clearBoys <- boys[boys$YYS == TRUE, ]
girls <- combinedRes[combinedRes$fetal2 == "NULL", ]
girls <- girls[!is.na(boys[, 1]), ]
girls$YYS <- unlist(apply(girls,1,function(x) { return (x[4] == x[85])} ))
clearGirls <- girls[girls$YYS == TRUE, ]


gestation14_B <- clearBoys[clearBoys[,91] == "无风险" & clearBoys$modifiedWZ == 14, c(1,10,67,68) ]
gestation15_B <- clearBoys[clearBoys[,91] == "无风险" & clearBoys$modifiedWZ == 15, c(1,10,67,68) ]
gestation16_B <- clearBoys[clearBoys[,91] == "无风险" & clearBoys$modifiedWZ == 16, c(1,10,67,68) ]
gestation17_B <- clearBoys[clearBoys[,91] == "无风险" & clearBoys$modifiedWZ == 17, c(1,10,67,68) ]
gestation18_B <- clearBoys[clearBoys[,91] == "无风险" & clearBoys$modifiedWZ == 18, c(1,10,67,68) ]
gestation19_B <- clearBoys[clearBoys[,91] == "无风险" & clearBoys$modifiedWZ == 19, c(1,10,67,68) ]
gestation20_B <- clearBoys[clearBoys[,91] == "无风险" & clearBoys$modifiedWZ == 20, c(1,10,67,68) ]

gestation14_B <- gestation14_B[!is.na(gestation14_B[,1]),]
gestation15_B <- gestation15_B[!is.na(gestation15_B[,1]),]
gestation16_B <- gestation16_B[!is.na(gestation16_B[,1]),]
gestation17_B <- gestation17_B[!is.na(gestation17_B[,1]),]
gestation18_B <- gestation18_B[!is.na(gestation18_B[,1]),]
gestation19_B <- gestation19_B[!is.na(gestation19_B[,1]),]
gestation20_B <- gestation20_B[!is.na(gestation20_B[,1]),]


gestation14_G <- clearGirls[clearGirls[,91] == "无风险" & clearGirls$modifiedWZ == 14, c(1,10,67,68) ]
gestation15_G <- clearGirls[clearGirls[,91] == "无风险" & clearGirls$modifiedWZ == 15, c(1,10,67,68) ]
gestation16_G <- clearGirls[clearGirls[,91] == "无风险" & clearGirls$modifiedWZ == 16, c(1,10,67,68) ]
gestation17_G <- clearGirls[clearGirls[,91] == "无风险" & clearGirls$modifiedWZ == 17, c(1,10,67,68) ]
gestation18_G <- clearGirls[clearGirls[,91] == "无风险" & clearGirls$modifiedWZ == 18, c(1,10,67,68) ]
gestation19_G <- clearGirls[clearGirls[,91] == "无风险" & clearGirls$modifiedWZ == 19, c(1,10,67,68) ]
gestation20_G <- clearGirls[clearGirls[,91] == "无风险" & clearGirls$modifiedWZ == 20, c(1,10,67,68) ]

gestation14_G <- gestation14_G[!is.na(gestation14_G[,1]),]
gestation15_G <- gestation15_G[!is.na(gestation15_G[,1]),]
gestation16_G <- gestation16_G[!is.na(gestation16_G[,1]),]
gestation17_G <- gestation17_G[!is.na(gestation17_G[,1]),]
gestation18_G <- gestation18_G[!is.na(gestation18_G[,1]),]
gestation19_G <- gestation19_G[!is.na(gestation19_G[,1]),]
gestation20_G <- gestation20_G[!is.na(gestation20_G[,1]),]


fileDir <- read.table("sampleList.txt", header = F)

getLocation <- function (inData,inFile) {
    inData$fileLocation <- unlist(apply(inData, 1, function(x) { 
        t <- grep(paste(".*",x[3],".*IonXpress_",sub("0.","",as.character(as.numeric(x[4])/1000)),sep=""), fileDir$V1, perl = T);
        if (length(t) > 0) { return(fileDir[row.names(fileDir) == t[1],1]) } else { return("0") }
    } ))
    inData <- inData[inData$fileLocation != "0",]
    inData$fileLocation <- paste(inData$fileLocation," /data/gestation/",inFile,"/",inData[,1],".bam ",sep = "");
    inData$fileLocation <- paste("scp bioadmin@193.168.110.204:",sub("./","//",inData$fileLocation),sep = "")
    write.table(inData,paste("gestation/",inFile,".txt",sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
    return(inData);
}

gestation14_B <- getLocation(gestation14_B, "B14")
gestation15_B <- getLocation(gestation15_B, "B15")
gestation16_B <- getLocation(gestation16_B, "B16")
gestation17_B <- getLocation(gestation17_B, "B17")
gestation18_B <- getLocation(gestation18_B, "B18")
gestation19_B <- getLocation(gestation19_B, "B19")
gestation20_B <- getLocation(gestation20_B, "B20")


gestation14_G <- getLocation(gestation14_G, "G14")
gestation15_G <- getLocation(gestation15_G, "G15")
gestation16_G <- getLocation(gestation16_G, "G16")
gestation17_G <- getLocation(gestation17_G, "G17")
gestation18_G <- getLocation(gestation18_G, "G18")
gestation19_G <- getLocation(gestation19_G, "G19")
gestation20_G <- getLocation(gestation20_G, "G20")

mkdir -p  B14/childer
mkdir -p  B14/maternal

mkdir -p  B15/childer
mkdir -p  B15/maternal

mkdir -p  B16/childer
mkdir -p  B16/maternal

mkdir -p  B17/childer
mkdir -p  B17/maternal

mkdir -p  B18/childer
mkdir -p  B18/maternal

mkdir -p  B19/childer
mkdir -p  B19/maternal

mkdir -p  B20/childer
mkdir -p  B20/maternal



cat B14.txt | awk -F "\t" '{ print($5) }' | head -n 200 > b14.sh
cat B15.txt | awk -F "\t" '{ print($5) }' | head -n 200 > b15.sh
cat B16.txt | awk -F "\t" '{ print($5) }' | head -n 200 > b16.sh
cat B17.txt | awk -F "\t" '{ print($5) }' | head -n 200 > b17.sh
cat B18.txt | awk -F "\t" '{ print($5) }' | head -n 200 > b18.sh
cat B19.txt | awk -F "\t" '{ print($5) }' | head -n 200 > b19.sh
cat B20.txt | awk -F "\t" '{ print($5) }' | head -n 200 > b20.sh
cat b*.sh > b_all.sh

mkdir -p  G14/childer
mkdir -p  G14/maternal

mkdir -p  G15/childer
mkdir -p  G15/maternal

mkdir -p  G16/childer
mkdir -p  G16/maternal

mkdir -p  G17/childer
mkdir -p  G17/maternal

mkdir -p  G18/childer
mkdir -p  G18/maternal

mkdir -p  G19/childer
mkdir -p  G19/maternal

mkdir -p  G20/childer
mkdir -p  G20/maternal

cat G14.txt | awk -F "\t" '{ print($5) }' | head -n 200 > g14.sh
cat G15.txt | awk -F "\t" '{ print($5) }' | head -n 200 > g15.sh
cat G16.txt | awk -F "\t" '{ print($5) }' | head -n 200 > g16.sh
cat G17.txt | awk -F "\t" '{ print($5) }' | head -n 200 > g17.sh
cat G18.txt | awk -F "\t" '{ print($5) }' | head -n 200 > g18.sh
cat G19.txt | awk -F "\t" '{ print($5) }' | head -n 200 > g19.sh
cat G20.txt | awk -F "\t" '{ print($5) }' | head -n 200 > g20.sh
cat g*.sh > g_all.sh


rm -rf /data/gestation/B14/childer/*.*
rm -rf /data/gestation/B14/maternal/*.*

rm -rf /data/gestation/B15/childer/*.*
rm -rf /data/gestation/B15/maternal/*.*

rm -rf /data/gestation/B16/childer/*.*
rm -rf /data/gestation/B16/maternal/*.*

rm -rf /data/gestation/B17/childer/*.*
rm -rf /data/gestation/B17/maternal/*.*

rm -rf /data/gestation/B18/childer/*.*
rm -rf /data/gestation/B18/maternal/*.*

rm -rf /data/gestation/B19/childer/*.*
rm -rf /data/gestation/B19/maternal/*.*

rm -rf /data/gestation/B20/childer/*.*
rm -rf /data/gestation/B20/maternal/*.*



plotForGene <- function(inGene, inDirL, inColor) {
    dataRes <- data.frame();
    for (i in 1:length(inDirL)) {
        a <- read.delim(paste(inDirL[i],"/",inGene,".bed", sep = ""), header = F);
        if (i == 1) { dataRes <- a; }
        else { dataRes <- cbind(dataRes, a); }
    }
    dataRes <- dataRes[,seq(4, 30, 6)];
    colnames(dataRes) <- inDirL;
    return(dataRes);
}