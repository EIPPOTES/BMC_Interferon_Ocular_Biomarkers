# 高级统计学方法学习笔记

*学习来源: Wikipedia | 学习时间: 2026-03-07*
*涵盖: 贝叶斯统计、机器学习、因果推断、Meta分析、诊断试验*

---

## 第一部分: 贝叶斯统计

### 1. 贝叶斯推断 (Bayesian Inference)

**核心思想**: 使用先验知识+新数据更新对假设的信念

**贝叶斯定理**:
```
P(H|E) = P(E|H) × P(H) / P(E)
```

**术语解释**:
| 术语 | 符号 | 含义 |
|------|------|------|
| 先验概率 | P(H) | 观察数据前对假设的信念 |
| 似然 | P(E|H) | 假设成立时观察到数据的概率 |
| 后验概率 | P(H|E) | 观察数据后对假设的更新信念 |
| 证据 | P(E) | 观察数据的边际概率 |

**与频率学派对比**:

| 特征 | 频率学派 | 贝叶斯学派 |
|------|---------|-----------|
| 概率解释 | 长期频率 | 信念程度 |
| 参数 | 固定未知 | 随机变量 |
| 推断 | 点估计+置信区间 | 后验分布 |
| 先验 | 不使用 | 使用先验分布 |

**眼科应用**:
- 诊断试验的概率更新
- 临床试验的适应性设计
- 个体化治疗决策

**R代码 (使用brms包)**:
```r
library(brms)

# 贝叶斯线性回归
brm(BCVA ~ treatment + age, data=df, family=gaussian())

# 贝叶斯Logistic回归
brm(disease ~ risk_factor, data=df, family=bernoulli())
```

### 2. MCMC方法 (Markov Chain Monte Carlo)

**用途**: 从复杂后验分布中抽样

**常用算法**:
1. **Metropolis-Hastings**: 通用MCMC算法
2. **Gibbs采样**: 条件分布已知时使用
3. **Hamiltonian Monte Carlo (HMC)**: 高效采样，Stan使用
4. **No-U-Turn Sampler (NUTS)**: HMC的改进版

**关键概念**:
- **收敛诊断**: Gelman-Rubin统计量 (R̂ < 1.1)
- **有效样本量 (ESS)**: 独立样本的等效数量
- **Burn-in**: 舍弃初始不稳定的样本

**R代码 (使用rstanarm)**:
```r
library(rstanarm)

# 贝叶斯回归
stan_glm(BCVA ~ treatment, data=df, family=gaussian())

# 检查收敛
model <- stan_glm(...)
summary(model, probs=c(0.025, 0.975))
```

---

## 第二部分: 机器学习

### 1. 随机森林 (Random Forest)

**原理**: 集成学习方法，构建多棵决策树并投票

**特点**:
- 处理高维数据
- 不易过拟合
- 提供特征重要性
- 可处理缺失值

**超参数**:
- `ntree`: 树的数量 (通常500-1000)
- `mtry`: 每次分裂随机选择的特征数
- `maxdepth`: 树的最大深度

**眼科应用**:
- OCT图像分类 (正常/病变)
- 疾病风险预测
- 特征选择 (识别重要生物标志物)

**R代码**:
```r
library(randomForest)

# 分类
rf_model <- randomForest(disease ~ ., data=train_data, 
                         ntree=500, mtry=3)

# 预测
predictions <- predict(rf_model, test_data)

# 特征重要性
importance(rf_model)
varImpPlot(rf_model)
```

### 2. 支持向量机 (SVM)

**原理**: 寻找最优超平面最大化分类间隔

**核函数**:
- **线性核**: 线性可分数据
- **多项式核**: 多项式决策边界
- **RBF核**: 非线性数据 (最常用)
- **Sigmoid核**: 类似神经网络

**超参数**:
- `C`: 惩罚参数 (小C=宽间隔，大C=严格分类)
- `gamma`: 核函数系数 (RBF核)

**眼科应用**:
- 眼底图像分类
- 疾病亚型识别
- 预后预测

**R代码**:
```r
library(e1071)

# SVM分类
svm_model <- svm(disease ~ ., data=train_data, 
                 kernel="radial", cost=1, gamma=0.1)

# 交叉验证调参
tune.svm(disease ~ ., data=train_data,
         kernel="radial",
         cost=c(0.1, 1, 10),
         gamma=c(0.01, 0.1, 1))
```

### 3. 深度学习在眼科的应用

**常用架构**:
- **CNN**: 图像分类 (眼底照片、OCT)
- **U-Net**: 图像分割 (视网膜层分割)
- **ResNet**: 深层网络，解决梯度消失
- **Transformer**: 注意力机制，Vision Transformer (ViT)

**眼科AI应用**:

| 应用领域 | 任务 | 代表工作 |
|---------|------|---------|
| 糖尿病视网膜病变 | 筛查分级 | Google DeepMind |
| AMD | 进展预测 | AREDS研究 |
| 青光眼 | 视神经分析 | OCT-AI系统 |
| 白内障 | 手术规划 | 自动IOL计算 |

**Python代码示例 (PyTorch)**:
```python
import torch
import torch.nn as nn

# 简单CNN
class SimpleCNN(nn.Module):
    def __init__(self):
        super(SimpleCNN, self).__init__()
        self.conv1 = nn.Conv2d(3, 32, 3)
        self.conv2 = nn.Conv2d(32, 64, 3)
        self.fc1 = nn.Linear(64*56*56, 128)
        self.fc2 = nn.Linear(128, 2)  # 二分类
    
    def forward(self, x):
        x = self.conv1(x)
        x = nn.ReLU()(x)
        x = nn.MaxPool2d(2)(x)
        x = self.conv2(x)
        x = nn.ReLU()(x)
        x = nn.MaxPool2d(2)(x)
        x = x.view(x.size(0), -1)
        x = self.fc1(x)
        x = nn.ReLU()(x)
        x = self.fc2(x)
        return x
```

---

## 第三部分: 因果推断

### 1. 倾向评分匹配 (Propensity Score Matching, PSM)

**用途**: 观察性研究中模拟随机化

**步骤**:
1. 估计倾向评分: P(Treatment=1 | X)
2. 匹配: 为每个处理组找到对照组匹配
3. 评估平衡性
4. 在匹配样本中估计效应

**匹配方法**:
- **最近邻匹配**: 1:1或1:N匹配
- **卡尺匹配**: 限制最大距离
- **核匹配**: 加权平均

**眼科应用**:
- 比较不同治疗方案的效果
- 真实世界研究

**R代码**:
```r
library(MatchIt)

# 估计倾向评分并匹配
m.out <- matchit(treatment ~ age + gender + baseline_VA,
                 data=df, method="nearest", ratio=1)

# 检查平衡性
summary(m.out)

# 获取匹配后数据
matched_data <- match.data(m.out)

# 分析
lm(outcome ~ treatment, data=matched_data)
```

### 2. 工具变量 (Instrumental Variable, IV)

**用途**: 处理内生性问题 (遗漏变量、反向因果、测量误差)

**有效工具变量的条件**:
1. **相关性**: 与内生变量相关
2. **外生性**: 与误差项不相关
3. **排他性**: 仅通过内生变量影响结果

**常用方法**:
- **两阶段最小二乘法 (2SLS)**
- **广义矩估计 (GMM)**

**眼科应用**:
- 处理治疗选择的混杂
- 基因作为工具变量 (孟德尔随机化)

**R代码**:
```r
library(AER)

# 2SLS
ivreg(outcome ~ treatment + controls | instrument + controls, 
      data=df)

# 检验弱工具变量
summary(iv_model, diagnostics=TRUE)
```

### 3. 双重差分 (Difference-in-Differences, DID)

**用途**: 政策/干预效果评估

**基本思想**:
```
效应 = (处理后处理组 - 处理前处理组) 
     - (处理后对照组 - 处理前对照组)
```

**假设**:
- **平行趋势假设**: 无干预时，两组趋势相同

**眼科应用**:
- 评估新筛查项目的效果
- 医保政策变化对眼病治疗的影响

**R代码**:
```r
# 基本DID模型
lm(outcome ~ treatment + post + treatment:post, data=panel_data)

# 或固定效应模型
library(fixest)
feols(outcome ~ treatment:post | patient_id + time, data=panel_data)
```

---

## 第四部分: Meta分析

### 1. Meta分析基础

**定义**: 对多个独立研究的结果进行定量综合

**步骤**:
1. 制定研究问题和纳入标准
2. 系统检索文献
3. 质量评估
4. 数据提取
5. 统计分析
6. 异质性检验
7. 发表偏倚评估

### 2. 效应量计算

**连续变量**:
- **标准化均数差 (SMD)**: Cohen's d, Hedges' g
- **均数差 (MD)**: 当量纲相同时

**二分类变量**:
- **比值比 (OR)**
- **相对危险度 (RR)**
- **危险度差 (RD)**

### 3. 统计模型

**固定效应模型**:
- 假设所有研究估计同一真实效应
- 使用: Mantel-Haenszel, Inverse Variance

**随机效应模型**:
- 假设研究间效应存在差异
- 使用: DerSimonian-Laird, REML

**选择标准**:
- 异质性低 (I² < 50%): 固定效应
- 异质性高 (I² ≥ 50%): 随机效应

### 4. 异质性评估

| 指标 | 解释 |
|------|------|
| Q统计量 | 卡方检验，p<0.1表示显著异质性 |
| I²统计量 | 0-100%，表示异质性比例 |
| τ² | 研究间方差 |

**I²解释**:
- 0-25%: 低异质性
- 25-50%: 中等异质性
- 50-75%: 高异质性
- 75-100%: 极高异质性

### 5. 发表偏倚

**检测方法**:
- **漏斗图 (Funnel Plot)**: 视觉检查
- **Egger检验**: 回归检验
- **剪补法 (Trim and Fill)**: 估计缺失研究

**R代码 (使用meta包)**:
```r
library(meta)

# Meta分析
meta_result <- metagen(TE=effect_size, seTE=se, 
                       studlab=study_name,
                       data=df,
                       sm="SMD",  # 或"OR"
                       method.tau="REML",
                       hakn=TRUE)

# 结果
summary(meta_result)

# 森林图
forest(meta_result)

# 漏斗图
funnel(meta_result)

# Egger检验
metabias(meta_result, method="linreg")

# 亚组分析
metagen(..., subgroup=treatment_type)
```

---

## 第五部分: 诊断试验评价

### 1. 基本指标

**四格表**:
```
              金标准阳性    金标准阴性
检测阳性      TP (真阳性)    FP (假阳性)
检测阴性      FN (假阴性)    TN (真阴性)
```

**核心指标**:

| 指标 | 公式 | 解释 |
|------|------|------|
| **敏感性 (Sensitivity)** | TP/(TP+FN) | 有病者中检测阳性比例 |
| **特异性 (Specificity)** | TN/(TN+FP) | 无病者中检测阴性比例 |
| **阳性预测值 (PPV)** | TP/(TP+FP) | 阳性者中真有病比例 |
| **阴性预测值 (NPV)** | TN/(TN+FN) | 阴性者中真无病比例 |
| **准确率 (Accuracy)** | (TP+TN)/N | 总正确率 |

### 2. ROC分析

**ROC曲线**:
- X轴: 1-特异性 (假阳性率)
- Y轴: 敏感性 (真阳性率)
- 每个点对应一个阈值

**AUC解释**:
| AUC | 诊断准确性 |
|-----|-----------|
| 0.5 | 无诊断价值 (随机) |
| 0.6-0.7 | 较低 |
| 0.7-0.8 | 一般 |
| 0.8-0.9 | 良好 |
| 0.9-1.0 | 优秀 |

**眼科应用**:
- OCT诊断黄斑水肿
- AI模型诊断糖尿病视网膜病变
- 生物标志物预测疾病进展

**R代码**:
```r
library(pROC)

# ROC分析
roc_obj <- roc(disease_status, test_score)

# AUC和置信区间
auc(roc_obj)
ci.auc(roc_obj)

# 绘制ROC曲线
plot(roc_obj, main=paste("AUC =", round(auc(roc_obj), 3)))

# 寻找最佳截断点
coords(roc_obj, "best", ret=c("threshold", "specificity", "sensitivity"))

# 比较两个ROC曲线
roc1 <- roc(disease, test1)
roc2 <- roc(disease, test2)
roc.test(roc1, roc2)  # DeLong检验
```

### 3. 似然比 (Likelihood Ratio)

**定义**:
- **阳性似然比 (LR+)**: 敏感性/(1-特异性)
- **阴性似然比 (LR-)**: (1-敏感性)/特异性

**解释**:
| LR+ | 解释 |
|-----|------|
| 1-2 | 几乎无变化 |
| 2-5 | 小幅度增加 |
| 5-10 | 中等增加 |
| >10 | 大幅度增加 |

**贝叶斯应用**:
```
验后 odds = 验前 odds × 似然比
```

---

## 第六部分: 方法选择指南

### 1. 研究问题与方法匹配

| 研究问题 | 推荐方法 |
|---------|---------|
| 两组均数比较 | t检验/Mann-Whitney U |
| 多组均数比较 | ANOVA/Kruskal-Wallis |
| 配对比较 | 配对t检验/Wilcoxon符号秩 |
| 二分类结局 | Logistic回归 |
| 时间-事件结局 | Cox回归 |
| 相关/聚类数据 | 混合效应模型/GEE |
| 观察性因果推断 | PSM/IV/DID |
| 综合多个研究 | Meta分析 |
| 诊断试验评价 | ROC分析 |
| 预测建模 | 随机森林/SVM/深度学习 |

### 2. 眼科研究特殊考虑

**双眼数据**:
- 首选: 混合效应模型
- 备选: GEE (大样本)

**影像数据**:
- 传统特征: 手工特征+机器学习
- 端到端: 深度学习 (CNN)

**纵向数据**:
- 线性混合模型 (连续结局)
- 广义估计方程 (非正态)
- 联合模型 (生存+纵向)

---

## 学习资源

### Wikipedia页面
- 贝叶斯推断: https://en.wikipedia.org/wiki/Bayesian_inference
- MCMC: https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo
- 随机森林: https://en.wikipedia.org/wiki/Random_forest
- SVM: https://en.wikipedia.org/wiki/Support_vector_machine
- 深度学习: https://en.wikipedia.org/wiki/Deep_learning
- 倾向评分: https://en.wikipedia.org/wiki/Propensity_score_matching
- 工具变量: https://en.wikipedia.org/wiki/Instrumental_variables_estimation
- 双重差分: https://en.wikipedia.org/wiki/Difference_in_differences
- Meta分析: https://en.wikipedia.org/wiki/Meta-analysis
- ROC分析: https://en.wikipedia.org/wiki/Receiver_operating_characteristic
- 敏感性/特异性: https://en.wikipedia.org/wiki/Sensitivity_and_specificity

---

*最后更新: 2026-03-07*
