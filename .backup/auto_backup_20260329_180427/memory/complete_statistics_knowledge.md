# 完整统计学知识库

*学习来源: Wikipedia | 学习时间: 2026-03-07*
*涵盖: 参数检验、非参数检验、回归分析、高级模型*

---

## 第一部分: 基础统计分布 (已完成)

详见: `statistics_tables_learning.md`

- 标准正态分布 (Z)
- t分布
- 卡方分布
- F分布

---

## 第二部分: 非参数统计

### 1. 非参数统计概述

**定义**: 对数据分布做最少假设的统计分析方法

**适用场景**:
- 数据不满足正态分布假设
- 数据为等级/顺序尺度
- 存在异常值
- 样本量小且分布未知

**与参数检验对比**:

| 特征 | 参数检验 | 非参数检验 |
|------|---------|-----------|
| 分布假设 | 需要 (如正态分布) | 不需要 |
| 数据类型 | 连续变量 | 连续、等级、分类 |
| 检验效能 | 较高 | 稍低 |
| 稳健性 | 对异常值敏感 | 对异常值稳健 |

### 2. Mann-Whitney U检验 (Wilcoxon秩和检验)

**用途**: 两组独立样本的比较 (非参数t检验替代)

**假设**:
1. 两组观测值相互独立
2. 数据至少为顺序尺度 (可排序)
3. H0: 两组分布相同
4. H1: 两组分布不同

**眼科应用**:
- 比较两种治疗方式的视力改善 (分布非正态时)
- 比较两组患者的OCT厚度 (存在异常值时)

**R代码**:
```r
# 基础用法
wilcox.test(group1, group2)

# 带连续性校正
wilcox.test(group1, group2, correct=TRUE)

# 单侧检验
wilcox.test(group1, group2, alternative="greater")
```

### 3. Kruskal-Wallis检验

**用途**: 多组独立样本比较 (非参数ANOVA替代)

**特点**:
- 单因素方差分析的秩次版本
- 可处理2组或以上
- 样本量可不等

**后续分析**:
- Dunn's test (推荐)
- Bonferroni校正的Mann-Whitney检验
- Conover-Iman test

**眼科应用**:
- 比较3种抗VEGF药物的疗效
- 比较不同疾病分期的视网膜厚度

**R代码**:
```r
# Kruskal-Wallis检验
kruskal.test(outcome ~ group, data=df)

# 事后检验 (需要dunn.test包)
library(dunn.test)
dunn.test(df$outcome, df$group, method="bonferroni")
```

### 4. Wilcoxon符号秩检验

**用途**: 配对样本比较 (非参数配对t检验替代)

**与Mann-Whitney区别**:
- Wilcoxon符号秩: 配对/相关样本
- Mann-Whitney U: 独立样本

**假设**:
- 差值分布对称
- 考虑差值的大小和方向

**眼科应用**:
- 治疗前后视力比较
- 双眼对比 (左眼 vs 右眼)

**R代码**:
```r
# 配对检验
wilcox.test(before, after, paired=TRUE)

# 单样本检验 (与理论中位数比较)
wilcox.test(data, mu=median_value)
```

---

## 第三部分: 多重比较校正

### 1. Bonferroni校正

**原理**: 控制族错误率 (Family-Wise Error Rate)

**公式**:
```
α_adjusted = α / m
```
- α: 期望的整体显著性水平
- m: 比较次数

**示例**:
- 进行20次比较，期望整体α=0.05
- 每次比较的显著性水平: 0.05/20 = 0.0025

**p值调整方法**:
```
p_adjusted = min(p × m, 1)
```

**优点**:
- 简单易懂
- 控制FWER严格

**缺点**:
- 过于保守，检验效能降低
- 比较次数多时难以拒绝H0

**R代码**:
```r
# p值向量
p_values <- c(0.01, 0.03, 0.05, 0.001)

# Bonferroni校正
p.adjust(p_values, method="bonferroni")
```

### 2. Tukey HSD检验

**用途**: ANOVA事后多重比较

**全称**: Tukey's Honestly Significant Difference

**特点**:
- 基于学生化极差分布 (q分布)
- 比较所有可能的组间配对
- 控制整体错误率

**适用条件**:
- 已完成ANOVA且显著
- 各组样本量相等或相近

**眼科应用**:
- ANOVA发现3组间有差异后，确定哪两组不同

**R代码**:
```r
# ANOVA
aov_result <- aov(BCVA ~ treatment, data=df)

# Tukey HSD
TukeyHSD(aov_result)

# 可视化
plot(TukeyHSD(aov_result))
```

---

## 第四部分: 回归分析

### 1. 线性回归

**简单线性回归**:
```
Y = β₀ + β₁X + ε
```

**多元线性回归**:
```
Y = β₀ + β₁X₁ + β₂X₂ + ... + βₖXₖ + ε
```

**假设**:
1. 线性关系
2. 误差项独立
3. 误差项正态分布
4. 误差项方差齐性
5. 无多重共线性

**眼科应用**:
- 预测视力与OCT参数的关系
- 分析年龄、性别对眼压的影响

**R代码**:
```r
# 简单线性回归
lm(BCVA ~ CST, data=df)

# 多元线性回归
lm(BCVA ~ CST + Age + Gender, data=df)

# 模型诊断
model <- lm(BCVA ~ CST, data=df)
summary(model)
plot(model)  # 残差图
```

### 2. Logistic回归

**用途**: 二分类结局分析

**模型**:
```
logit(p) = ln(p/(1-p)) = β₀ + β₁X₁ + ... + βₖXₖ
```

**OR值解释**:
- OR = 1: 无关联
- OR > 1: 危险因素
- OR < 1: 保护因素

**眼科应用**:
- 预测疾病发生风险
- 分析治疗效果 (有效/无效)
- 识别危险因素

**R代码**:
```r
# 二分类Logistic回归
glm(disease ~ age + gender + smoking, 
    family=binomial(link="logit"), data=df)

# 获取OR和95%CI
model <- glm(disease ~ age, family=binomial, data=df)
exp(cbind(coef(model), confint(model)))
```

### 3. Cox比例风险模型

**用途**: 生存分析/时间-事件分析

**模型**:
```
h(t|X) = h₀(t) × exp(β₁X₁ + ... + βₖXₖ)
```

**关键概念**:
- **风险比 (HR)**: 类似OR，用于时间-事件数据
- **比例风险假设**: 各组风险比不随时间变化
- **基线风险函数 h₀(t)**: 所有协变量为0时的风险

**眼科应用**:
- 疾病进展时间分析
- 治疗失败时间分析
- 视力丧失时间分析

**R代码**:
```r
library(survival)

# 创建生存对象
surv_obj <- Surv(time, event)

# Cox回归
coxph(Surv(time, event) ~ treatment + age, data=df)

# 检验比例风险假设
cox.zph(model)

# Kaplan-Meier曲线
fit <- survfit(Surv(time, event) ~ treatment, data=df)
plot(fit)
```

---

## 第五部分: 高级模型

### 1. 混合效应模型 (Mixed Model / LMM)

**用途**: 处理相关/聚类数据

**模型结构**:
```
Y = Xβ + Zb + ε
```
- Xβ: 固定效应
- Zb: 随机效应
- ε: 随机误差

**适用场景**:
- 重复测量数据
- 双眼数据 (患者内相关)
- 多中心研究
- 层次/嵌套数据

**眼科应用**:
- 双眼分析 (患者为随机效应)
- 纵向随访数据 (时间点嵌套于患者)
- 多中心RCT (中心为随机效应)

**R代码**:
```r
library(lme4)

# 简单随机截距模型
lmer(BCVA ~ treatment + (1|patient_id), data=df)

# 随机斜率模型
lmer(BCVA ~ time + (time|patient_id), data=df)

# 模型摘要
model <- lmer(BCVA ~ treatment + (1|patient_id), data=df)
summary(model)

# 获取随机效应
ranef(model)
```

### 2. 广义估计方程 (GEE)

**用途**: 纵向/聚类数据的边际模型

**特点**:
- 关注总体平均效应 (population-averaged)
- 对相关性结构误设稳健
- 使用三明治标准误 (sandwich estimator)

**与混合模型区别**:
| 特征 | GEE | 混合模型 |
|------|-----|---------|
| 效应类型 | 总体平均 | 个体特定 |
| 随机效应 | 无 | 有 |
| 缺失数据 | 需要MCAR | 更灵活 |
| 解释 | 平均效应 | 个体效应 |

**眼科应用**:
- 大样本纵向研究
- 当关注总体效应而非个体变化时

**R代码**:
```r
library(geepack)

# GEE模型
gee_model <- geeglm(BCVA ~ treatment + time,
                    id=patient_id,
                    family=gaussian,
                    corstr="exchangeable",
                    data=df)

# 相关性结构选项:
# "independence" - 独立
# "exchangeable" - 可交换
# "ar1" - 一阶自相关
# "unstructured" - 无结构
```

---

## 第六部分: 统计功效与样本量

### 1. 统计功效 (Power)

**定义**: 当备择假设为真时，正确拒绝H0的概率

**公式**:
```
Power = 1 - β
```
- β: II型错误概率 (假阴性)

**影响因素**:
1. 效应量 (Effect size) - 越大越好
2. 样本量 - 越大越好
3. 显著性水平 (α) - 越大越好 (但通常固定)
4. 检验类型 - 单侧>双侧

**标准**:
- 通常要求 Power ≥ 0.80
- 高质量研究要求 Power ≥ 0.90

### 2. 样本量计算

**两均数比较**:
```
n = 2 × [(Zα/₂ + Zβ)² × σ²] / δ²
```

**两率比较**:
```
n = [Zα/₂√(2p̄(1-p̄)) + Zβ√(p₁(1-p₁)+p₂(1-p₂))]² / (p₁-p₂)²
```

**常用Z值**:
| 参数 | 值 |
|------|-----|
| Z₀.₀₅ (双侧) | 1.96 |
| Z₀.₀₁ (双侧) | 2.576 |
| Z₀.₂₀ (Power=80%) | 0.84 |
| Z₀.₁₀ (Power=90%) | 1.28 |

**R代码**:
```r
library(pwr)

# t检验样本量
pwr.t.test(d=0.5, sig.level=0.05, power=0.8)

# 两比例检验
pwr.2p.test(h=0.3, sig.level=0.05, power=0.8)

# 效应量计算
cohen.d <- function(x, y) {
  mean_diff <- mean(x) - mean(y)
  pooled_sd <- sqrt((var(x) + var(y)) / 2)
  return(mean_diff / pooled_sd)
}
```

---

## 第七部分: 眼科研究特殊考虑

### 1. 双眼数据处理

**方案对比**:

| 方案 | 方法 | 优点 | 缺点 |
|------|------|------|------|
| 1 | 单眼分析 | 简单 | 损失信息，增加I型错误 |
| 2 | 平均双眼 | 简单 | 忽略眼间差异 |
| 3 | 混合效应模型 | 正确处理相关性 | 复杂 |
| 4 | GEE | 稳健，大样本适用 | 需要大样本 |

**推荐**: 方案3 (混合效应模型)

### 2. 重复测量数据

**分析方法选择**:
- 2个时间点: 配对t检验或Wilcoxon符号秩
- 3+个时间点: 重复测量ANOVA或混合效应模型
- 非等时间间隔: 混合效应模型

### 3. 常见统计错误避免

**错误1**: p值滥用
- ✅ 报告效应量和CI
- ❌ 仅报告p<0.05

**错误2**: 多重检验未校正
- ✅ 使用Bonferroni或FDR
- ❌ 多次t检验不校正

**错误3**: 忽视双眼相关性
- ✅ 使用混合效应模型
- ❌ 将双眼视为独立

**错误4**: 样本量不足
- ✅ 研究前计算样本量
- ❌ 事后解释不显著结果

---

## 学习资源

### Wikipedia页面
- 非参数统计: https://en.wikipedia.org/wiki/Nonparametric_statistics
- Mann-Whitney U: https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
- Kruskal-Wallis: https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_test
- Wilcoxon符号秩: https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
- Bonferroni校正: https://en.wikipedia.org/wiki/Bonferroni_correction
- Tukey HSD: https://en.wikipedia.org/wiki/Tukey%27s_range_test
- 线性回归: https://en.wikipedia.org/wiki/Linear_regression
- Logistic回归: https://en.wikipedia.org/wiki/Logistic_regression
- Cox模型: https://en.wikipedia.org/wiki/Cox_proportional_hazards_model
- 混合效应模型: https://en.wikipedia.org/wiki/Mixed_model
- GEE: https://en.wikipedia.org/wiki/Generalized_estimating_equation
- 统计功效: https://en.wikipedia.org/wiki/Statistical_power

---

*最后更新: 2026-03-07*
*下次学习目标: 贝叶斯统计、机器学习在眼科的应用*
