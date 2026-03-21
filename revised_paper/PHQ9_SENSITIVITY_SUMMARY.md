# PHQ-9分层敏感性分析总结

**分析日期**: 2026-03-11
**目的**: 评估排除PHQ-9<5患者后主要发现的稳健性

---

## 📊 分析设计

### 分组
| 组别 | PHQ-9评分 | 眼数 | 占比 |
|------|-----------|------|------|
| 活跃抑郁 | ≥5 | 157眼 | 60.4% |
| 缓解期 | <5 | 103眼 | 39.6% |
| **总计** | - | **260眼** | **100%** |

### 分析方法
- 活跃抑郁 vs 对照
- 缓解期 vs 对照
- 比较效应量（β系数）和显著性

---

## 📈 主要结果

| 结局变量 | 活跃抑郁 β (P) | 缓解期 β (P) | 所有MDD β (P) | 结论 |
|----------|----------------|--------------|---------------|------|
| Mean Macular Thickness | -5.463 (0.001) | -6.689 (0.004) | -5.714 (0.000) | 两组均显著 |
| Outer Temporal Thickness | -6.390 (0.000) | -10.110 (0.000) | -7.383 (0.000) | 两组均显著 |
| Total Macular Volume | -0.158 (0.001) | -0.189 (0.004) | -0.163 (0.000) | 两组均显著 |
| CD Ratio | 0.002 (0.913) | 0.031 (0.227) | 0.007 (0.717) | 两组均不显著 |
| Rim Volume | 0.002 (0.939) | -0.020 (0.396) | -0.006 (0.725) | 两组均不显著 |

---

## 🔍 关键发现

### 1. 主要OCT参数在两组均显著
- **Mean Macular Thickness**: 活跃期β=-5.463，缓解期β=-6.689
- **Outer Temporal Thickness**: 活跃期β=-6.390，缓解期β=-10.110
- **Total Macular Volume**: 活跃期β=-0.158，缓解期β=-0.189

### 2. 缓解期效应量更大
- Outer Temporal Thickness: 缓解期效应量（β=-10.110）> 活跃期（β=-6.390）
- 提示视网膜改变可能不是单纯由急性症状驱动

### 3. 视盘参数在两组均不显著
- CD Ratio和Rim Volume在分层分析中均不显著
- 与主要分析结果一致

---

## 💡 解释与意义

### 主要结论
> 抑郁与视网膜结构改变的关联在**活跃抑郁**和**缓解期**患者中均存在，提示这些改变可能是**特质标志物（trait markers）**而非与急性症状严重程度相关的**状态改变（state-dependent changes）**。

### 临床意义
1. **诊断价值**: 即使在症状缓解期，视网膜改变仍可能存在
2. **机制研究**: 提示结构性改变可能反映长期的神经生物学改变
3. **治疗监测**: 视网膜厚度可能不适合作为短期疗效指标

---

## 📝 论文中添加的文本

已在Results 3.1节添加：

```markdown
**Sensitivity analysis by depression severity**: We conducted stratified analyses 
comparing (1) active depression patients (PHQ-9 ≥ 5, n=157 eyes) vs controls, and 
(2) remission patients (PHQ-9 < 5, n=103 eyes) vs controls. Both groups showed 
significant reductions in mean macular thickness (active: β=-5.463, P=0.001; 
remission: β=-6.689, P=0.004) and outer temporal thickness (active: β=-6.390, 
P<0.001; remission: β=-10.110, P<0.001) compared to controls. These findings 
suggest that retinal structural changes may represent trait markers of depression 
vulnerability rather than state-dependent changes associated with acute symptom 
severity.
```

---

## ✅ 审稿意见回应

| 审稿意见 | 回应 |
|----------|------|
| PHQ-9<5患者分类有争议 | ✅ 进行分层分析，证明两组均显著 |
| 需要排除PHQ-9<5的敏感性分析 | ✅ 完成分析，结果稳健 |
| 讨论trait vs state markers | ✅ 在Results中明确讨论 |

---

## 📁 相关文件

| 文件 | 位置 |
|------|------|
| `Sensitivity_Analysis_PHQ9_Stratified.xlsx` | `04_原始数据/` |
| 更新后的论文 | `manuscript_revised_with_Table5.md` |
