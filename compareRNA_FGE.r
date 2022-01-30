# compare the gene expression in different gestation week for body
setwd("/Users/hegs/Desktop/Projects/fetalGenes")

hgncGenes <- read.delim("HGNC_reliable.bed") # nolint

nrowA <- rep("", 9)
for (i in 1:nrow(hgncGenes)) {
    if (hgncGenes[i, 4] == "-") {
        fileName <- paste(hgncGenes[i, c(2, 3, 8)], collapse = "_")
    } # nolint
    else {
        fileName <- paste(hgncGenes[i, c(2, 3, 7)], collapse = "_")
    } # nolint
    if (file.exists(paste("Data/", fileName, ".bed", sep = ""))) {
        temp <- read.delim(paste("Data/", fileName, ".bed", sep = ""))
        up <- which(temp$TSS == -150)
        down <- which(temp$TSS == 50)
        flagL <- seq(-2500, -1000, 5)
        flagR <- seq(1000, 2500, 5)
        if (i == 1) {
            res <- temp[, c(length(colnames(temp)), 1:(length(colnames(temp))))]
            colnames(res)[1] <- "GeneName"
            res <- res[-1:(-nrow(res)), ]
        }
        res <- rbind(res, nrowA)
        res[nrow(res), 1] <- hgncGenes[i, 2]
        res[nrow(res), 2] <- fileName
        for (j in 1:(length(colnames(temp)) - 1)) {
            flage_m <- mean(c(temp[temp$TSS %in% flagL, j], temp[temp$TSS %in% flagR, j]))
            center_m <- mean(temp[up:down, j])
            res[nrow(res), j + 2] <- round(center_m / flage_m, 5)
        }
    }
}
colnames(res) <- c("geneName", "fileName", paste("B", seq(14, 20), sep = ""), "mean")

for (i in 3:length(colnames(res))) {
    res[, i] <- as.numeric(res[, i])
}

res.h <- res[which(res$mean < 0.2), ]
# res.e <- res[which(res$mean != 0),]
library(pheatmap)
pheatmap(res.h[, 3:9])


res.h <- res[which(res$mean == 0), ]

res.e$cv <- res.e$pvalue <- res.e$cor <- 0
for (i in 1:nrow(res.e)) {
    t <- cor.test(t(res.e[i, 3:9]), 14:20)
    res.e$cor[i] <- round(t$estimate, 3)
    res.e$pvalue[i] <- t$p.value
    res.e$cv[i] <- round(sd(res.e[i, 3:9]) / mean(t(res.e[i, 3:9])) * 100, 2)
}
# PCA analysis



# Get high correlation genes.
res.matrix <- res.e[abs(res.e$cor) > 0.4 & res.e$pvalue < 0.05, c(1, 3:9)]
# res.matrix <- res.e[, c(1,3:9)]
row.names(res.matrix) <- res.matrix[, 1]
res.matrix <- res.matrix[, -1]

library(Mfuzz)
set.seed(123)
res.matrix <- data.matrix(res.matrix)
res.matrix <- new("ExpressionSet", exprs = res.matrix)
res.matrix <- standardise(res.matrix)

## 聚类个数
c <- 7
## 计算最佳的m值
m <- mestimate(res.matrix)
## 聚类
res.matrix.class <- mfuzz(res.matrix, c = c, m = m)
mfuzz.plot(res.matrix, res.matrix.class, mfrow = c(2, 3), new.window = FALSE)

# library(TCseq)

#
# cluster_num <- 6
# tcseq_cluster <- timeclust(res.matrix, algo = 'cm', k = cluster_num, standardize = TRUE)

# p <- timeclustplot(tcseq_cluster, value = 'z-score', cols = 3,
#    axis.line.size = 0.6, axis.title.size = 8, axis.text.size = 8,
#    title.size = 8, legend.title.size = 8, legend.text.size = 8)


## 过滤缺失超过25%的基因
gene.r <- filter.NA(eset, thres = 0.25)




eset <- new("ExpressionSet", exprs = data.matrix(res.f))
yeastF <- filter.NA(yeast)
gene.s <- standardise(est)