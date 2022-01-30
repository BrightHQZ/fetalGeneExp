library(getopt)
library(ggpubr)

# 第一列：参数的长名称（多个字符） 
# 第二列：参数的短名称（一个字符）
# 第三列：这个flag对应的参数形式（0表示flag不接受参数；1表示可接可不接；2表示必须接参数）
# 第四列：参数的类型（logical；integer；double；complex；character；numeric
# 第五列：注释信息

# 首先将第一个参数的具体信息进行描述
spec <- matrix(
    # 每行五个，第五个可选，也就是说第五列可以不
    # byrow 按行填充矩阵的元
    # ncol  每行填充五个元素
    c(
        "rootDir", "r", 2, "character", "The root dir of storage data!",
        "from", "f", 2, "integer", "The first NO of sub dir with storage data!",
        "to", "t", 2, "integer", "The last NO of sub dir with storage data!",
        "sex", "s", 2, "character", "The fetus gender b/g!",
        "bed", "b", 2, "character", "The gene list which need caculated!",
        "cBI", "p", 2, "character", "Process type : all/bed/image!"
    ),
    byrow = TRUE, ncol = 5
)

# 使用getopt方法，注意这里的usage默认是关闭的，这里不能打开
opt <- getopt(spec = spec)

if (is.null(opt$rootDir) || is.null(opt$from) || is.null(opt$to) || is.null(opt$HGNC) || is.null(opt$cBI)) {
    # ... 这里你也可以自定义一些东放在里面
    cat(paste(getopt(spec = spec, usage = T), "\n"))
    quit()
}

if (!opt$cBI %in% c("all", "bed", "image")) {
    print("all/bed/image is need for program working!")
    cat(paste(getopt(spec = spec, usage = T), "\n"))
    quit()
}

checkFile <- function(infileL, inRoot, inDirL) {
    for (i in 1:length(inDirL)) {
        # print(paste(inRoot,"/",inDirL[i],"/",fileL,".bed", sep = ""));
        if (!file.exists(paste(inRoot, "/", inDirL[i], "/", infileL, ".bed", sep = ""))) {
            return("")
        }
    }
    return(infileL)
}

searchDataForGene <- function(inGene, inRoot, inDirL) {
    dataRes <- data.frame()
    for (i in 1:length(inDirL)) {
        a <- read.delim(paste(inRoot, "/", inDirL[i], "/", inGene, ".bed", sep = ""), header = F)
        if (i == 1) {
            dataRes <- a
        } else {
            dataRes <- cbind(dataRes, a)
        }
    }
    dataRes <- dataRes[, seq(4, 42, 6)]
    colnames(dataRes) <- inDirL
    dataRes$Median <- unlist(apply(dataRes, 1, sum)) / length(inDirL)
    res <- dataRes
    for (i in 1:(length(inDirL) + 1)) {
        res[, i] <- round(dataRes[, i] / mean(dataRes[, i]), 3)
    }
    res$TSS <- seq(-2500, 2500, 5)
    return(res)
    # return(dataRes);
}

checkDestroy <- function(inRes) {
    for (i in (length(colnames(inRes)) - 1)) {
        if (nrow(inRes[inRes[, i] == 0, ]) / nrow(inRes) > 0.3) {
            return(FALSE)
        }
    }
    return(TRUE)
}



# col <- c("black","red","green","blue","orange","drown")
if (toupper(opt$sex) == "B") {
    dirList <- paste("B", opt$from:opt$to, sep = "")
}
if (toupper(opt$sex) == "G") {
    dirList <- paste("G", opt$from:opt$to, sep = "")
}

if ( file.exists(opt$bed) ) {
    hgnc <- read.delim(opt$bed, header = F)
} else {
    print(paste("The gene list of scanning is not exists:", opt$bed, sep = " "));
    q();
}

if (!file.exists(paste(opt$rootDir, "/Data", sep = ""))) {
    dir.create(paste(opt$rootDir, "/Data", sep = ""))
}
if (!file.exists(paste(opt$rootDir, "/Image", sep = ""))) {
    dir.create(paste(opt$rootDir, "/Image", sep = ""))
}



destoryList <- c()
for (i in 1:nrow(hgnc)) {
    fileN <- checkFile(hgnc[i,4], opt$rootDir, dirList)
    if (fileN != "") {
        temp <- searchDataForGene(fileN, opt$rootDir, dirList)
        if (checkDestroy(temp)) {
            if (opt$cBI %in% c("all", "bed")) {
                write.table(temp, paste(opt$rootDir, "/Data/", fileN, ".bed", sep = ""), sep = "\t", row.names = F, quote = F)
                print(paste("The file of ", fileN, " is processed!", sep = ""))
            }
            # temp <- temp[temp$TSS >= -25 & temp$TSS <= 25, ]
            # temp$TSS <- temp$TSS * 100;
            if (opt$cBI %in% c("all", "image")) {
                p <- list()
                for (j in 1:length(dirList)) {
                    p[[j]] <- ggline(temp, x = "TSS", y = colnames(temp)[j], plot_type = "l", color = "steelblue") +
                        scale_x_discrete(breaks = seq(-2000, 2000, 1000)) + font("xy.text", size = 8) +
                        ylab(paste("Relative intensity (B", (13 + j), ")", sep = "")) + font("ylab", size = 10) +
                        geom_vline(aes(xintercept = 501), linetype = "dashed", color = "red", size = 0.6)
                }
                p[[length(dirList) + 1]] <- ggline(temp, x = "TSS", y = colnames(temp)[length(dirList) + 1], plot_type = "l", color = "#bd3005ce") +
                    scale_x_discrete(breaks = seq(-2000, 2000, 1000)) + font("xy.text", size = 7) +
                    geom_vline(aes(xintercept = 501), linetype = "dashed", color = "#5833ddd7", size = 0.6)

                ggarrange(p[[1]], p[[5]], p[[2]], p[[6]], p[[3]], p[[7]], p[[4]], p[[8]], ncol = 2, nrow = 4) %>%
                    ggexport(filename = paste(opt$rootDir, "/Image/", fileN, ".png", sep = ""), width = 1458, height = 1880, res = 250)
            }
        } else {
            destoryList <- c(destoryList, fileN)
            print(paste("Destoried the file of ", fileN, "!", sep = ""))
        }
    }
}

write.table(destoryList, paste(opt$rootDir, "/Data/destoryList.txt", sep = ""), sep = "\t", row.names = F, quote = F)
# a <- read.delim("Data/AMACR.bed")
# a <- a[41:61,1]
# b <- c();
# center <- c();
# if (length(a)/2 == 1) {
#    center <- length(a)%/%2 + 1
# } else { center[1] <- (length(a)%/%2) + 1}

# for (i in 1:length(a)) {
#    if (i == center) { b <- c(b,a[center]) }
#    else if (i < center) {
#        b <- c(b,((a[i] + a[i+1])/2));
#    } else if (i > center) {
#        b <- c(b,((a[i] + a[i-1])/2));
#    }
# }