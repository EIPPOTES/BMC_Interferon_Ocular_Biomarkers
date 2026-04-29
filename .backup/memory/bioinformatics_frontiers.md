# 生物信息学前沿热点

*学习来源: Wikipedia | 学习时间: 2026-03-07*
*针对: 眼科与视网膜疾病研究*

---

## 第一部分: 生物信息学概述

### 1. 什么是生物信息学

**定义**: 生物信息学是一门跨学科领域，开发计算方法和软件工具来理解生物数据，特别是当数据集庞大而复杂时。

**核心领域**:
- **基因组学 (Genomics)**: DNA序列分析
- **转录组学 (Transcriptomics)**: RNA表达分析
- **蛋白质组学 (Proteomics)**: 蛋白质研究
- **代谢组学 (Metabolomics)**: 代谢物分析
- **表观基因组学 (Epigenomics)**: 表观遗传修饰

**在眼科研究中的应用**:
- 遗传性视网膜疾病的基因定位
- AMD、青光眼的多基因风险评分
- 眼部肿瘤的分子分型
- 药物靶点发现

---

## 第二部分: 单细胞测序 (Single-Cell Sequencing)

### 1. 技术概述

**定义**: 对单个细胞的核酸序列信息进行测序，提供细胞差异的高分辨率视图

**主要类型**:
- **scRNA-seq**: 单细胞转录组测序
- **scDNA-seq**: 单细胞基因组测序
- **scATAC-seq**: 单细胞染色质可及性测序
- **CITE-seq**: 单细胞表面蛋白测序

### 2. 技术优势

| 传统Bulk测序 | 单细胞测序 |
|-------------|-----------|
| 数百万细胞混合 | 单个细胞分辨率 |
| 平均信号 | 细胞异质性 |
| 丢失稀有细胞信息 | 发现稀有细胞类型 |
| 无法追踪细胞状态 | 细胞轨迹推断 |

### 3. 分析流程

```
单细胞分离 → 细胞裂解 → mRNA捕获 → cDNA合成 → 
文库制备 → 测序 → 数据分析
```

### 4. 核心分析方法

#### 4.1 降维与可视化
- **PCA**: 主成分分析
- **t-SNE**: t-分布随机邻域嵌入
- **UMAP**: 均匀流形近似和投影 (推荐)

#### 4.2 聚类分析
- 识别细胞类型
- 标记基因鉴定
- 细胞亚群分类

#### 4.3 轨迹推断
- **Monocle**: 伪时间分析
- **Slingshot**: 细胞分化轨迹
- **RNA velocity**: RNA速率分析

#### 4.4 细胞通讯分析
- **CellPhoneDB**: 配体-受体分析
- **CellChat**: 细胞间通讯推断

### 5. 眼科应用

**视网膜研究**:
- 视网膜细胞图谱构建
- 感光细胞发育轨迹
- 退行性病变中的细胞变化

**角膜研究**:
- 角膜缘干细胞鉴定
- 角膜上皮细胞异质性

**眼肿瘤**:
- 视网膜母细胞瘤细胞异质性
- 葡萄膜黑色素瘤微环境

### 6. 常用工具

| 工具 | 用途 | 语言 |
|------|------|------|
| **Seurat** | 单细胞分析全流程 | R |
| **Scanpy** | 单细胞分析全流程 | Python |
| **Cell Ranger** | 10x Genomics数据处理 | CLI |
| **Monocle3** | 轨迹推断 | R |
| **SingleR** | 细胞类型注释 | R |

### 7. R代码示例

```r
library(Seurat)

# 读取数据
data <- Read10X(data.dir = "filtered_gene_bc_matrices/")

# 创建Seurat对象
seu <- CreateSeuratObject(counts = data, project = "Retina", min.cells = 3, min.features = 200)

# 质控
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 标准化
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

# 降维
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunUMAP(seu, dims = 1:10)

# 聚类
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)

# 可视化
DimPlot(seu, reduction = "umap", label = TRUE)
FeaturePlot(seu, features = c("RHO", "OPN1SW", "GFAP"))

# 标记基因
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

---

## 第三部分: 空间转录组学 (Spatial Transcriptomics)

### 1. 技术概述

**定义**: 在完整组织中捕获转录活性的位置信息的方法

**技术类型**:

| 类型 | 代表技术 | 分辨率 | 特点 |
|------|---------|--------|------|
| 基于测序 | 10x Visium | 55μm | 全转录组 |
| 基于成像 | MERFISH, seqFISH | 单细胞/亚细胞 | 高分辨率 |
| 原位捕获 | Slide-seq, Stereo-seq | 10μm | 高分辨率 |

### 2. 技术优势

- **空间信息**: 保留组织空间结构
- **细胞环境**: 了解细胞邻域关系
- **组织结构**: 识别空间表达模式
- **疾病微环境**: 肿瘤-免疫相互作用

### 3. 分析方法

#### 3.1 空间聚类
- 识别空间域 (Spatial Domains)
- 空间可变基因鉴定

#### 3.2 细胞类型去卷积
- **SPOTlight**: 解卷积空间spot
- **Cell2location**: 细胞类型定位

#### 3.3 细胞-细胞相互作用
- 空间邻域分析
- 配体-受体共定位

### 4. 眼科应用

**视网膜**:
- 视网膜层特异性基因表达
- 黄斑 vs 周边视网膜差异
- 退行性病变的空间模式

**角膜**:
- 角膜上皮分层基因表达
- 角膜缘干细胞微环境

**眼肿瘤**:
- 肿瘤-免疫微环境
- 侵袭前沿基因表达

### 5. 常用工具

| 工具 | 用途 |
|------|------|
| **Seurat (Spatial)** | 空间数据分析 |
| **Squidpy** | 空间分析Python包 |
| **Giotto** | 空间转录组分析 |
| **STUtility** | 空间数据可视化 |

---

## 第四部分: CRISPR基因编辑

### 1. 技术概述

**定义**: 成簇规律间隔短回文重复序列 (Clustered Regularly Interspaced Short Palindromic Repeats)

**核心组件**:
- **Cas9**: 核酸酶，切割DNA
- **sgRNA**: 单向导RNA，引导Cas9到靶位点
- **PAM序列**: Cas9识别位点 (NGG)

### 2. CRISPR应用类型

| 应用 | 说明 | 用途 |
|------|------|------|
| **基因敲除** | 破坏基因功能 | 功能研究 |
| **基因插入** | 插入新序列 | 基因治疗 |
| **碱基编辑** | 单碱基替换 | 点突变修正 |
| **表观编辑** | 调控基因表达 | 可逆调控 |
| **CRISPR筛选** | 全基因组筛选 | 靶点发现 |

### 3. CRISPR筛选 (Perturb-seq)

**原理**: 结合CRISPR基因编辑和单细胞测序

**流程**:
```
sgRNA文库 → 细胞感染 → 单细胞捕获 → 
同时测sgRNA和转录组 → 基因功能分析
```

### 4. 眼科应用

**遗传性视网膜疾病**:
- ABCA4基因编辑治疗Stargardt病
- RHO基因编辑治疗视网膜色素变性
- CEP290编辑治疗Leber先天性黑矇

**疾病模型**:
- 构建视网膜疾病iPSC模型
- 基因功能研究

**药物靶点筛选**:
- 全基因组筛选保护视网膜的基因
- 筛选AMD易感基因

---

## 第五部分: 蛋白质组学 (Proteomics)

### 1. 技术概述

**定义**: 大规模研究蛋白质的科学

**主要技术**:
- **质谱 (MS)**: 蛋白质鉴定和定量
- **蛋白质芯片**: 高通量蛋白质检测
- **酵母双杂交**: 蛋白质相互作用

### 2. 质谱技术类型

| 技术 | 特点 | 应用 |
|------|------|------|
| **DDA** | 数据依赖采集 | 蛋白质鉴定 |
| **DIA** | 数据非依赖采集 | 高通量定量 |
| **TMT/iTRAQ** | 同位素标记 | 多组比较 |
| **LFQ** | 非标记定量 | 大样本队列 |

### 3. 眼科应用

**泪液蛋白质组**:
- 干眼症生物标志物
- 眼部感染诊断

**房水/玻璃体**:
- AMD发病机制
- 糖尿病视网膜病变

**视网膜组织**:
- 视网膜疾病蛋白质变化
- 光感受器蛋白质组

---

## 第六部分: 宏基因组学 (Metagenomics)

### 1. 技术概述

**定义**: 研究特定环境中所有生物的遗传物质

**应用**:
- 微生物组研究
- 环境基因组学
- 病原体检测

### 2. 眼科微生物组研究

**眼表微生物组**:
- 结膜菌群与干眼症
- 角膜感染病原体检测
- 隐形眼镜相关微生物

**分析流程**:
```
样本采集 → DNA提取 → 16S rRNA测序/宏基因组测序 → 
生物信息学分析 → 功能注释
```

### 3. 常用工具

| 工具 | 用途 |
|------|------|
| **QIIME2** | 微生物组分析 |
| **MetaPhlAn** | 物种分类 |
| **HUMAnN** | 功能分析 |
| **LEfSe** | 差异分析 |

---

## 第七部分: 多组学整合分析

### 1. 整合策略

**水平整合**:
- 同一样本的多组学数据
- 基因组 + 转录组 + 蛋白质组

**垂直整合**:
- 不同层次生物学信息
- DNA → RNA → 蛋白质 → 代谢物

### 2. 整合方法

| 方法 | 用途 |
|------|------|
| **MOFA+** | 多组学因子分析 |
| **SNF** | 相似性网络融合 |
| **mixOmics** | 多组学整合 |
| **XGBoost/RF** | 机器学习整合 |

### 3. 眼科多组学研究

**AMD研究**:
- GWAS + 转录组 + 蛋白质组
- 多组学风险评分

**青光眼**:
- 基因组 + 眼压 + 视野
- 整合预测模型

---

## 第八部分: 生物信息学数据库与资源

### 1. 基因组数据库

| 数据库 | 内容 | 网址 |
|--------|------|------|
| **NCBI/GenBank** | 序列数据 | https://www.ncbi.nlm.nih.gov |
| **Ensembl** | 基因组注释 | https://www.ensembl.org |
| **UCSC Genome Browser** | 基因组浏览 | https://genome.ucsc.edu |
| **gnomAD** | 人群变异 | https://gnomad.broadinstitute.org |

### 2. 眼科专用数据库

| 数据库 | 内容 |
|--------|------|
| **RetNet** | 遗传性视网膜疾病基因 |
| **EyeGEx** | 眼组织eQTL数据 |
| **Ocular Tissue Database** | 眼组织基因表达 |

### 3. 单细胞数据库

| 数据库 | 内容 |
|--------|------|
| **Cell Atlas** | 人类细胞图谱 |
| **Single Cell Portal** | Broad Institute |
| **ASCT+B Reporter** | 解剖结构细胞类型 |

---

## 第九部分: 眼科生信研究前沿

### 1. 当前热点 (2024-2025)

**单细胞视网膜图谱**:
- 人类视网膜细胞图谱完善
- 疾病状态下的细胞变化
- 发育过程中的细胞分化

**空间转录组应用**:
- AMD黄斑区基因表达
- 视网膜色素上皮异质性
- 脉络膜新生血管微环境

**基因治疗**:
- CRISPR治疗遗传性视网膜疾病
- AAV载体优化
- 基因编辑安全性评估

**AI+生信**:
- 深度学习预测基因功能
- 多组学数据整合预测
- 药物重定位

### 2. 推荐学习路径

**入门**:
1. R/Python编程基础
2. 生物统计学
3. Linux命令行

**进阶**:
1. RNA-seq数据分析
2. 单细胞测序分析
3. 基因组变异分析

**高级**:
1. 多组学整合
2. 机器学习应用
3. 网络生物学

---

## 学习资源

### 在线课程
- Coursera: Genomic Data Science
- edX: Introduction to Bioconductor
- Cold Spring Harbor: Single Cell Analysis

### 书籍
- 《Bioinformatics and Functional Genomics》
- 《Single-Cell Analysis》
- 《Computational Genomics with R》

### 网站
- Bioconductor: https://www.bioconductor.org
- Biostars: https://www.biostars.org
- Single Cell Blog: https://satijalab.org

---

*最后更新: 2026-03-07*
