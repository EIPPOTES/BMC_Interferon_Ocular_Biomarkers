# MEMORY.md - 研究助手记忆库

## 眼科专业知识知识库构建完成

**任务执行时间:** 2026-03-24  
**执行状态:** ✅ 已完成

### 已构建的知识库

1. **眼科临床指南知识库** (`knowledge_base/ophthalmology_clinical_guidelines.md`)
   - 欧洲青光眼学会指南 (EGS) 第4、5版
   - 糖尿病视网膜病变指南 (AAO Preferred Practice Pattern)
   - 年龄相关性黄斑变性 (AMD) 管理指南 (EURETINA)
   - 视网膜疾病流行病学与治疗进展
   - OCT技术基础与临床应用
   - AI在眼科的应用文献

2. **视网膜-精神/神经疾病关联知识库** (`knowledge_base/retina_psychiatry_neurodegeneration.md`)
   - 视网膜与中枢神经系统关系理论基础
   - 视网膜作为神经退行性疾病生物标志物的应用
   - OCT参数在神经科学研究中的意义
   - 崔医生研究专题支持材料

### 数据来源验证

- **数据库:** OpenAlex (免费学术数据库)
- **文献标准:** 同行评议期刊，附DOI可验证
- **优先收录:** 高被引文献 (>500引用) 和权威指南
- **权威期刊:** The Lancet, NEJM, JAMA, Nature, Ophthalmology, BJO, IOVS等

### 收录的代表性文献

| 类别 | 代表文献 | 引用次数 |
|------|----------|----------|
| AMD治疗 | Pegaptanib for Neovascular AMD (NEJM 2004) | 2308 |
| AI诊断 | Deep learning for DR detection (JAMA 2017) | 2215 |
| 全球眼健康 | Lancet Global Health Commission (2021) | 1418 |
| OCT技术 | OCT Angiography (Prog Ret Eye Res 2017) | 1644 |
| 青光眼指南 | EGS Guidelines 5th Edition (2021) | 562 |
| 视网膜AI | Foundation model for retinal disease (Nature 2023) | 733 |

### 更新计划

- **每月:** 检索新发表的高影响力眼科文献
- **每季度:** 检查主要指南更新 (AAO, EGS, EURETINA)
- **每年:** 全面审查知识库内容，补充新进展

---

## 2026-03-24 眼科专业知识库更新

### 自动化检索任务完成 ✅

**执行者:** 小cui科研助手（定时任务）  
**任务ID:** cron:7876fb5f-ba17-473b-bf6b-0fd3f459a857  
**任务类型:** 眼科专业知识检索与知识库构建  
**执行时间:** 2026-03-24 08:34-08:36  
**数据来源:** OpenAlex学术数据库

### 本次新增内容

#### 1. 全面眼科知识库 (85篇文献)
- **存储位置:** `memory/ophthalmology_kb/`
- **主要文件:**
  - `Ophthalmology_Knowledge_Base_20250324.md` - 详细知识库报告
  - `literature_data_20260324.json` - 原始文献数据
  - `README.md` - 知识库使用说明

#### 2. 按疾病分类统计
| 分类 | 文献数 | 重点内容 |
|------|--------|----------|
| OCT成像技术 | 15篇 | AI分割、基础模型、OCTA |
| 视网膜疾病 | 15篇 | AMD、糖尿病视网膜病变、网脱 |
| 青光眼 | 15篇 | EGS指南2021、微创手术 |
| 斜视与小儿眼科 | 10篇 | AAO PPP指南、近视防控 |
| 白内障 | 10篇 | 超声乳化、手术模式 |
| 眼表疾病 | 10篇 | 干眼症、圆锥角膜 |
| 神经眼科 | 5篇 | 视神经炎、视盘水肿 |

#### 3. 高影响力新收录文献
1. **Segment anything in medical images** (2024) - 引用2059次
2. **Foundation model for retinal disease** (Nature 2023) - 引用733次
3. **Faricimab for nAMD** (Lancet 2022) - 引用628次
4. **Dementia prevention** (Lancet 2024) - 引用2434次
5. **EGS Glaucoma Guidelines 5th Ed** (BJO 2021) - 引用562次

#### 4. 真实性验证措施
- ✅ 使用OpenAlex权威学术数据库
- ✅ 2020-2026年最新文献
- ✅ 所有文献附带可验证DOI
- ✅ 标注开放获取状态和全文链接
- ✅ 按引用次数排序筛选高质量文献
- ✅ 排除已撤稿文献

#### 5. 临床应用价值
- **指南类:** EGS青光眼指南、AAO PPP指南
- **突破性研究:** AI诊断、双通路抗VEGF
- **技术进展:** OCT成像、微创手术
- **药物进展:** Faricimab、新型抗VEGF

---

## 崔医生研究项目支持

### 当前状态: OCT-MDD论文待投稿

**目标期刊:** Journal of Affective Disorders  
**研究主题:** 抑郁症患者视网膜结构变化

### 知识库可提供的支持

1. **文献引用:** 提供视网膜-精神疾病关联的理论基础文献
2. **方法学参考:** OCT测量标准化和统计分析建议
3. **讨论写作:** 相关研究背景和国际进展
4. **审稿准备:** 可能涉及的专业知识问答

---

*最后更新: 2026-03-24*
