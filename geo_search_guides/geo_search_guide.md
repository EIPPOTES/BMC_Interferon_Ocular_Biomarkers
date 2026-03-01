# 眼科GEO数据库全流程使用指南

**作者**: Based on user-provided professional guidance
**日期**: 2026-03-01
**目的**: 帮助眼科研究者高效使用NCBI GEO数据库

## GEO核心编号系统

- GSE (Series - 系列/项目): 最重要！整个实验项目的总称
- GSM (Sample - 样本): 单个样本的数据
- GPL (Platform - 平台): 测序或芯片的技术平台
- GDS (DataSet - 数据集): 官方整理的精选数据集

## 高级检索策略

- 布尔逻辑: AND, OR, NOT
- 字段限定: [Title], [Title/Abstract], [Organism], [Entry Type]
- 眼科检索公式示例:
-   '"Glaucoma"[Title/Abstract] AND "Trabecular meshwork"[All Fields] AND "Homo sapiens"[Organism] AND "gse"[Entry Type]'

## GSE页面关键信息

- Title & Summary: 快速判断研究内容
- Overall design: 实验分组情况
- Contributor(s): 数据上传作者，可找原始论文
- Platforms: 技术平台信息
- Samples: 所有GSM样本列表

## 数据下载指南

- Series Matrix File(s): 最推荐！已处理的表达矩阵
- SOFT formatted family file(s): 元数据，较少使用
- Supplementary file(s): 原始数据或计数矩阵

## GEO2R使用指南

- Define groups: 创建分组（如Disease和Control）
- Assign samples: 分配样本到对应组
- Analyze: 执行差异分析
- 查看结果: 关注P.Value和logFC列

## 眼科经典数据集推荐

- AMD: GSE29801, GSE115828
- Glaucoma: GSE27276, GSE138114
- Diabetic Retinopathy: GSE60436
- Keratoconus: GSE112155, GSE15020
- Dry Eye: GSE43671

## 实用技巧

- 使用具体疾病英文全称搜索，比"Ophthalmology"更精准
- 务必加上"Homo sapiens"[Organism]限定人类数据
- 查看原始论文了解实验设计和样本信息
- 下载Series Matrix File进行初步分析

## 搜索示例

### 年龄相关性黄斑变性 (AMD)
- **搜索式**: `"Age-related Macular Degeneration"[Title] OR "AMD"[Title] AND "Homo sapiens"[Organism] AND "gse"[Entry Type]`
- **预期结果**: GSE29801, GSE115828等

### 青光眼
- **搜索式**: `"Glaucoma"[Title/Abstract] AND "Trabecular meshwork"[All Fields] AND "Homo sapiens"[Organism]`
- **预期结果**: GSE27276等

### 糖尿病视网膜病变
- **搜索式**: `"Diabetic retinopathy"[Title] AND "fibrovascular membrane"[All Fields] AND "Homo sapiens"[Organism]`
- **预期结果**: GSE60436等

## 与本框架集成

**目的**: 将GEO数据与本分析框架结合使用

**步骤**:
1. 1. 使用本指南搜索合适的眼科GEO数据集
2. 2. 下载Series Matrix File或计数矩阵
3. 3. 使用本框架的`process_real_data.py`处理数据
4. 4. 运行`analysis_pipeline.py`进行分析
5. 5. 使用`optimized_analysis.py`进行优化分析
6. 6. 生成完整的分析报告

**优势**:
- 标准化的工作流程
- 优化的分析方法
- 透明的结果报告
- 可重现的分析过程
