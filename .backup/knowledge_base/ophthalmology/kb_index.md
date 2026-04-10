# 眼科专业知识库索引

## 知识库概述
- **创建日期**: 2026-03-23
- **更新频率**: 每月更新
- **验证状态**: 基于OpenAlex学术数据库
- **来源验证**: 所有文献均来自同行评审期刊

## 疾病分类结构

### 1. 青光眼 (Glaucoma)
- 原发性开角型青光眼 (POAG)
- 原发性闭角型青光眼 (PACG)
- 正常眼压性青光眼 (NTG)
- 继发性青光眼

### 2. 白内障 (Cataract)
- 年龄相关性白内障
- 先天性白内障
- 代谢性白内障
- 外伤性白内障

### 3. 视网膜疾病 (Retinal Diseases)
- 糖尿病视网膜病变 (DR)
- 年龄相关性黄斑变性 (AMD)
- 视网膜静脉阻塞 (RVO)
- 视网膜脱离 (RD)
- 视网膜母细胞瘤

### 4. 玻璃体疾病 (Vitreous Diseases)
- 玻璃体混浊
- 玻璃体出血
- 玻璃体后脱离

### 5. 视神经疾病 (Optic Nerve Diseases)
- 视神经炎
- 缺血性视神经病变
- 视神经萎缩

### 6. 眼表疾病 (Ocular Surface Diseases)
- 干眼综合征
- 角膜炎
- 结膜炎

## 数据来源

### 主要数据库
- OpenAlex (250M+ 学术作品)
- PubMed (生物医学文献)
- AAO Preferred Practice Patterns

### 权威期刊
- Ophthalmology
- JAMA Ophthalmology
- American Journal of Ophthalmology
- British Journal of Ophthalmology
- Retina
- Progress in Retinal and Eye Research

### 指南来源
- American Academy of Ophthalmology (AAO)
- European Society of Ophthalmology (SOE)
- International Council of Ophthalmology (ICO)

## 知识库文件结构

```
knowledge_base/ophthalmology/
├── kb_index.md                 # 本索引文件
├── SOURCES_VERIFICATION.md     # 来源验证记录
├── quick_reference.md          # 快速参考卡
└── diseases/
    ├── 01_glaucoma.md          # 青光眼
    ├── 02_retinal_diseases.md  # 视网膜疾病
    ├── 03_cataract.md          # 白内障
    └── 04_oct_imaging.md       # OCT成像技术
```

## 文档统计
- **总文档数**: 4个疾病专题 + 2个辅助文档
- **涵盖疾病**: 青光眼、糖尿病视网膜病变、AMD、白内障、RVO、OCT技术
- **引用文献**: 15+篇高引文献
- **数据验证**: ✅ 已通过OpenAlex验证

## 使用指南

### 快速查阅
使用 `quick_reference.md` 查看关键参数速查表

### 深入学习
查阅 `diseases/` 目录下的专题文档

### 验证来源
查看 `SOURCES_VERIFICATION.md` 了解文献来源和验证状态

## 更新日志

| 日期 | 更新内容 | 验证状态 |
|------|---------|---------|
| 2026-03-23 | 知识库初始化，创建4个疾病专题 | ✅ 已验证 |
