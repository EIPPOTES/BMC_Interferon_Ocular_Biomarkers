
# R分析GEO数据的示例代码
library(GEOquery)
library(limma)

# 下载GSE数据
gse <- getGEO("GSE29801", GSEMatrix = TRUE)

# 获取表达矩阵和样本信息
exprs <- exprs(gse[[1]])
pdata <- pData(gse[[1]])

# 差异表达分析
design <- model.matrix(~0 + factor(pdata$group))
colnames(design) <- c("control", "case")
fit <- lmFit(exprs, design)
contrast.matrix <- makeContrasts(case-control, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取结果
results <- topTable(fit2, number=1000, adjust.method="fdr")
