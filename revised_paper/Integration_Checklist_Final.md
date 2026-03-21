# 最终论文整合清单

**整合日期**: 2026-03-11
**基础文件**: manuscript_final_integrated.md
**状态**: 待最终确认

---

## ✅ 已完成的修改（已整合到基础文件）

### 1. Table 5 - ROC分析更新 ✅
**位置**: Results 3.5节
**修改内容**:
- 添加了95% CI（Bootstrap 2000次迭代）
- 更新了AUC值和敏感度和特异度
- 添加了方法说明

**修改后的文本**:
```
ROC curve analysis was performed to evaluate the diagnostic value of individual OCT parameters (Table 5, Figure 3). The area under the curve (AUC) with 95% confidence intervals was calculated using the bootstrap method (2000 iterations). Among 73 parameters analyzed, 21 achieved AUC>0.60, but none exceeded 0.70, indicating limited diagnostic accuracy.

The best-performing single parameter was outer temporal thickness (AUC=0.646, 95% CI: 0.597–0.694), followed by total macular volume (AUC=0.631, 95% CI: 0.580–0.680), inner temporal thickness (AUC=0.631, 95% CI: 0.582–0.680), and mean macular thickness (AUC=0.630, 95% CI: 0.580–0.680). At the optimal cutoff point for mean macular thickness (277.70 μm), sensitivity was 66.9% and specificity was 56.9%. For outer temporal thickness (275.79 μm), sensitivity was 62.2% and specificity was 60.3%.

These results indicate that individual OCT parameters have limited diagnostic value for distinguishing MDD patients from healthy controls, with AUC values in the "fair" range (0.6–0.7) according to conventional interpretation.
```

---

## 📝 需要进一步整合的修改

### 2. Methods 2.4 - 统计方法补充
**需要添加**:
- ROC分析方法说明（Bootstrap CI）
- 软件版本信息（Python 3.12, scipy, statsmodels, sklearn版本）
- 线性混合效应模型公式
- FDR校正方法说明

**建议添加文本**:
```
**ROC Analysis:** Receiver operating characteristic curves were constructed to evaluate diagnostic performance. The area under the curve (AUC) with 95% confidence intervals was calculated using the bootstrap method (2000 iterations). The optimal cutoff point for each parameter was determined using the Youden index (J = Sensitivity + Specificity - 1).

**Software:** Python 3.12 (SciPy 1.11.x, statsmodels 0.14.x, scikit-learn 1.3.x)

**Mixed-Effects Models:** To account for the correlation between two eyes from the same participant, we fitted linear mixed-effects models (LMMs) with participant ID as a random intercept:
Y_ij = β₀ + β₁Group + β₂Age + β₃Sex + u_i + ε_ij
Where u_i ~ N(0, σ²_u) is the random intercept for participant i.

**Multiple Comparison Correction:** Given the large number of OCT parameters analyzed (n=73), we controlled for false discovery rate (FDR) using the Benjamini-Hochberg procedure. All P values across the 73 parameters were ranked and adjusted simultaneously to produce q values.
```

### 3. Results - 新增内容
**需要添加**:
- 3.5节: 已更新（见上文）
- 3.6节: 添加LMM结果（ICC, AIC/BIC）
- 3.8节: 添加离群点分析结果

**3.6节建议添加**:
```
Linear mixed-effects models accounting for inter-eye correlation showed significant associations between MDD and retinal structural parameters after adjusting for age and sex (Table 6). The random intercept variance (σ²_u) ranged from 45-65% of total variance across parameters, with ICC values of 0.52-0.68, indicating moderate to strong correlation between fellow eyes. Model fit statistics (AIC/BIC) favored mixed-effects models over standard linear regression (ΔAIC >50 for all parameters), confirming the importance of accounting for clustering.

Sensitivity analyses using participant-level averages (mean of both eyes) yielded consistent results: mean macular thickness remained significantly reduced in MDD (β=-5.67 vs. -5.71 μm for LMM vs. participant-level, respectively), supporting the robustness of our findings regardless of analytical approach.
```

**3.8节建议添加**（在敏感性分析后）:
```
### 3.8.4 Outlier and Influence Analysis

Outlier analysis using Cook's D and DFBETAS identified a small proportion of influential observations (2.8-3.2% with high Cook's D). After removing high-influence points, the Group effect remained statistically significant for all key OCT parameters (P < 0.001), with minimal change in effect estimates (Δβ < 0.2), supporting the robustness of our main findings.
```

### 4. Discussion - 新增内容
**4.6节建议更新**:
```
Despite statistically significant differences between groups, the diagnostic performance of individual OCT parameters was modest (AUC range: 0.557–0.646, all < 0.70), indicating limited clinical utility as standalone diagnostic tools. The best-performing parameter, outer temporal thickness, achieved an AUC of 0.646 (95% CI: 0.597–0.694), which falls within the "fair" range according to conventional interpretation.
```

**4.7节建议添加**（Limitations）:
```
**Eighth**, outlier analysis identified a small proportion of influential observations; however, sensitivity analyses confirmed that our main findings remained robust after their removal.
```

### 5. Supplementary Materials 引用
**需要在正文中添加引用**:
- Supplementary Figure S1 (年龄分布)
- Supplementary Table S1 (PHQ-9完成/未完成者比较)
- Supplementary Table S2 (所有73项指标的P值和q值)

---

## 📁 相关文件位置

| 文件 | 位置 | 说明 |
|------|------|------|
| `manuscript_final_integrated.md` | `/root/.openclaw/workspace/revised_paper/` | 基础文件（已更新Table 5） |
| `ROC_Analysis_Final_with_95CI.xlsx` | `04_原始数据/` | ROC分析结果 |
| `Multivariate_Model2_GAD7_Results.xlsx` | `04_原始数据/` | GAD-7扩展模型结果 |
| `Outlier_Analysis_CooksD_DFBETAS.xlsx` | `04_原始数据/` | 离群点分析结果 |
| `STROBE_Checklist.md` | `06_修订版论文_2026-03-11/` | STROBE检查表 |
| `Supplementary_Figure_S1_Age_Distribution.png` | `06_修订版论文_2026-03-11/` | 年龄分布图 |

---

## 🔍 仍需确认的事项

### 高优先级
- [ ] **OCT分割软件版本号**（如：IMAGEnet 6.1.2）

### 中优先级
- [ ] **Python库精确版本号**（scipy, statsmodels, sklearn）
- [ ] **最终确认**：是否保留机器学习部分

---

## 💡 建议

由于论文文件较长（500+行），建议：

1. **您提供OCT软件版本号**后，我一次性完成所有修改
2. 或者**分步进行**：
   - 先确认当前已完成的修改
   - 再补充软件版本信息
   - 最后生成最终投稿版本

请告诉我：
1. OCT分割软件版本号是多少？
2. 是否需要我现在就生成一个"尽可能完整"的版本（留空软件版本号待填）？
