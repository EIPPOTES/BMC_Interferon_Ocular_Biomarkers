# OCT-MDD研究核心参考文献速查指南

**生成时间**: 2026-03-23  
**知识库版本**: v6.17  
**适用对象**: 崔医生 OCT-MDD 论文撰写

---

## 🎯 方法学金标准（必读）

### 1. 视网膜Meta分析方法学

**Petzold 2010 - MS视网膜OCT Meta分析** ⭐⭐⭐⭐⭐
```
Petzold A, de Boer JF, Schippling S, et al.
Optical coherence tomography in multiple sclerosis: 
a systematic review and meta-analysis.
Lancet Neurology. 2010;9(9):863-873.
DOI: 10.1016/s1474-4422(10)70168-x
Citations: 582
```
**核心价值**: 
- MS视网膜Meta分析金标准方法学
- OCT参数选择和统计方法可直接借鉴
- 为OCT-MDD未来Meta分析提供方法学模板

**Petzold 2017 - 视网膜分层Meta分析** ⭐⭐⭐⭐⭐
```
Petzold A, Balcer LJ, Calabresi PA, et al.
Retinal layer segmentation in multiple sclerosis: 
a systematic review and meta-analysis.
Lancet Neurology. 2017;16(10):797-812.
DOI: 10.1016/s1474-4422(17)30278-8
Citations: 547
Open Access: ✅ https://research.manchester.ac.uk/
```
**核心价值**:
- 视网膜分层（RNFL, GCL, IPL等）Meta分析方法
- 异质性分析和亚组分析策略
- 可用于OCT-MDD分层分析参考

---

## 🔍 OCT质量控制标准

### 2. OSCAR-IB共识 - OCT质量评估

**Tewarie 2012** ⭐⭐⭐⭐⭐
```
Tewarie P, Balk LJ, Costello F, et al.
The OSCAR-IB Consensus Criteria for Retinal OCT Quality Assessment.
PLoS ONE. 2012;7(4):e34823.
DOI: 10.1371/journal.pone.0034823
Citations: 587
Open Access: ✅ https://journals.plos.org/plosone/
```

**核心要点**:
| 标准 | 内容 | 应用建议 |
|-----|------|---------|
| **O** - Obvious problems | 图像伪影、信号强度<5 | 排除图像 |
| **S** - Score < 5 | 信号强度评分 | 纳入/排除阈值 |
| **C** - Centration | 扫描中心偏差 | 黄斑中心凹定位 |
| **A** - Algorithm failure | 自动分层错误 | 人工校正 |
| **R** - Retinal pathology | 视网膜病变影响 | 分层准确性 |
| **I** - Inadequate overlap | 重复扫描一致性 | 双眼数据 |
| **B** - Beam placement | 光束位置 | 扫描协议 |

**论文应用**:
- 方法学部分：描述OCT图像质量控制流程
- 附录：提供详细纳入/排除标准

---

## 🧠 神经退行性机制参考

### 3. 视网膜神经退行性机制

**Simó 2018 - DR神经退行性** ⭐⭐⭐⭐
```
Simó R, Stitt AW, Gardner TW.
Neurodegeneration in diabetic retinopathy: does it really matter?
Diabetologia. 2018;61(9):1902-1912.
DOI: 10.1007/s00125-018-4692-1
Citations: 549
Open Access: ✅ https://link.springer.com/
```

**核心发现**:
- 糖尿病视网膜病变不仅是微血管病，更是神经退行性疾病
- 早期神经节细胞凋亡先于微血管改变
- 视网膜神经纤维层（RNFL）变薄可作为神经退行性标志

**论文应用**:
- Discussion: 支持"视网膜作为中枢神经系统窗口"的理论
- 引用视网膜神经退行性改变早于结构性改变的观点
- 与抑郁症神经退行性假说建立联系

---

## 🤖 AI分析方法学参考

### 4. 大规模AI筛查验证

**Ruamviboonsuk 2019 - 深度学习DR筛查** ⭐⭐⭐⭐
```
Ruamviboonsuk P, Krause J, Chotcomwongse P, et al.
Deep learning versus human graders for classifying diabetic 
retinopathy severity in a nationwide screening program.
npj Digital Medicine. 2019;2:99.
DOI: 10.1038/s41746-019-0099-8
Citations: 218
Open Access: ✅ https://www.nature.com/
```

**方法学亮点**:
- 25,326张图像的大规模验证
- 全国筛查项目实际应用场景
- AI vs 人工分级的一致性评估

**论文应用**:
- 未来OCT-MDD多中心研究设计参考
- AI辅助诊断方法学依据

### 5. 可解释AI设计指南

**Chen 2022 - 医疗AI可解释性** ⭐⭐⭐⭐
```
Chen H, Gómez C, Huang CM, et al.
Explainable medical imaging AI needs human-centered design: 
guidelines and evidence from a systematic review.
npj Digital Medicine. 2022;5:156.
DOI: 10.1038/s41746-022-00699-2
Citations: 219
Open Access: ✅ https://www.nature.com/
```

**核心要点**:
- 医疗AI需要"以人为本"的可解释设计
- 透明度是可操作性（affordance）而非模型属性
- 提供可解释AI系统设计的循证指南

**论文应用**:
- 如果论文涉及AI分析，需参考可解释性要求
- 算法透明度描述方法

### 6. DeepDR-LLM系统

**Li 2024 - 深度学习+大语言模型** ⭐⭐⭐⭐
```
Li J, Guan Z, Wang J, et al.
Integrated image-based deep learning and language models 
for primary diabetes care.
Nature Medicine. 2024;30(10):2936-2946.
DOI: 10.1038/s41591-024-03139-8
Citations: 174
Open Access: ✅ https://www.nature.com/
```

**技术亮点**:
- 图像深度学习与大语言模型整合
- 低资源环境下的糖尿病护理解决方案
- 多模态数据融合方法

**论文应用**:
- 未来OCT-MDD研究的技术前瞻性参考
- 多模态数据分析方法（OCT+临床数据）

---

## 📊 糖尿病视网膜病变参考

### 7. UK DR/DME管理共识

**Amoaku 2020 - UK共识** ⭐⭐⭐⭐
```
Amoaku WMK, Ghanchi F, Bailey C, et al.
Diabetic retinopathy and diabetic macular oedema pathways 
and management: UK Consensus Working Group.
Eye. 2020;34(Suppl 1):1-51.
DOI: 10.1038/s41433-020-0961-6
Citations: 197
Open Access: ✅ https://www.nature.com/
```

**应用价值**:
- 疾病管理框架描述参考
- 对照疾病（DR）的临床背景介绍

### 8. Cochrane OCT检测DME评价

**Virgili 2015 - Cochrane系统评价** ⭐⭐⭐⭐
```
Virgili G, Menchini F, Casazza G, et al.
Optical coherence tomography (OCT) for detection of macular oedema 
in patients with diabetic retinopathy.
Cochrane Database Syst Rev. 2015;1:CD008081.
DOI: 10.1002/14651858.cd008081.pub3
Citations: 263
Open Access: ✅ https://www.cochranelibrary.com/
```

**应用价值**:
- Cochrane系统评价方法学参考
- OCT诊断准确性评估标准
- 未来OCT-MDD系统评价方法学依据

---

## 🌍 全球眼健康背景

### 9. 眼健康差异综述

**Elam 2022 - 眼健康差异** ⭐⭐⭐
```
Elam AR, Tseng VL, Rodriguez TM, et al.
Disparities in Vision Health and Eye Care.
Ophthalmology. 2022;129(10):e107-e114.
DOI: 10.1016/j.ophtha.2022.07.010
Citations: 170
Open Access: ✅ https://pmc.ncbi.nlm.nih.gov/
```

**应用价值**:
- Introduction: 全球眼健康背景
- 强调抑郁症患者眼健康筛查的重要性

---

## 📋 论文各章节应用建议

### Introduction
- [ ] 引用 **Elam 2022** - 眼健康差异背景
- [ ] 引用 **Simó 2018** - 视网膜神经退行性机制
- [ ] 引用 **Petzold 2010** - OCT在神经退行疾病中的应用

### Methods
- [ ] 引用 **Tewarie 2012** - OCT质量控制（OSCAR-IB标准）
- [ ] 引用 **Virgili 2015** - OCT诊断准确性评估方法

### Results
- [ ] 可参考 **Petzold 2017** - 视网膜分层分析方法描述

### Discussion
- [ ] 引用 **Simó 2018** - 神经退行性机制解释
- [ ] 引用 **Petzold 2010** - 与MS研究结果比较
- [ ] 引用 **Chen 2022** - 未来AI分析的可解释性考虑

### Future Directions
- [ ] 引用 **Li 2024** - AI+LLM整合技术前景
- [ ] 引用 **Ruamviboonsuk 2019** - 大规模筛查项目设计

---

## 📚 快速获取链接

| 文献 | PDF链接 | 获取状态 |
|-----|--------|---------|
| Petzold 2017 | https://research.manchester.ac.uk/... | ✅ 开放获取 |
| Tewarie 2012 | https://journals.plos.org/plosone/... | ✅ 开放获取 |
| Simó 2018 | https://link.springer.com/content/pdf/... | ✅ 开放获取 |
| Ruamviboonsuk 2019 | https://www.nature.com/articles/... | ✅ 开放获取 |
| Chen 2022 | https://www.nature.com/articles/... | ✅ 开放获取 |
| Li 2024 | https://www.nature.com/articles/... | ✅ 开放获取 |
| Amoaku 2020 | https://www.nature.com/articles/... | ✅ 开放获取 |
| Virgili 2015 | https://www.cochranelibrary.com/... | ✅ 开放获取 |
| Elam 2022 | https://pmc.ncbi.nlm.nih.gov/articles/... | ✅ 开放获取 |
| Petzold 2010 | - | 🔒 需订阅 |

---

## 🎯 阅读优先级建议

| 优先级 | 文献 | 原因 |
|-------|------|------|
| 🔴 最高 | Petzold 2010, 2017 | 方法学金标准 |
| 🔴 最高 | Tewarie 2012 | 质控标准必需 |
| 🟡 高 | Simó 2018 | 机制解释关键 |
| 🟡 高 | Chen 2022 | AI可解释性参考 |
| 🟢 中 | Ruamviboonsuk 2019 | 未来研究设计 |
| 🟢 中 | Li 2024 | 技术前瞻性 |
| ⚪ 参考 | Amoaku 2020, Virgili 2015 | 背景和方法学补充 |

---

## 📁 完整报告位置

- 详细检索报告: `knowledge_base/ophthalmology/retrieval_report_2026-03-23-2132.md`
- 综合知识库: `knowledge_base/ophthalmology/ophthalmology_knowledge_base_2025-03-23.md`
- AAO EyeWiki汇编: `knowledge_base/ophthalmology/aao_eyewiki_disease_compendium.md`

---

*生成时间: 2026-03-23 21:50*  
*下次更新: 2026-03-30*
