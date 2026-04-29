# 专业级统计学方法学习笔记

*学习来源: Wikipedia | 学习时间: 2026-03-07*
*涵盖: 生存分析高级主题、缺失数据、网络Meta分析、可解释AI、临床试验设计*

---

## 第一部分: 生存分析高级主题

### 1. 竞争风险 (Competing Risks)

**问题**: 当存在多种互斥终点时，传统生存分析会高估感兴趣事件的概率

**示例**:
- 研究AMD进展为晚期时，患者可能先死于其他原因
- 研究白内障手术效果时，患者可能先发生视网膜脱离

**分析方法**:

#### 1.1 特定原因风险函数 (Cause-Specific Hazard)
```
h_k(t) = lim(Δt→0) P(t ≤ T < t+Δt, K=k | T≥t) / Δt
```

#### 1.2 累积发生率函数 (Cumulative Incidence Function, CIF)
```
F_k(t) = P(T ≤ t, K = k)
```
- 考虑竞争事件的影响
- 各事件CIF之和等于总体失效概率

#### 1.3 子分布风险 (Subdistribution Hazard)
- Fine-Gray模型
- 直接建模CIF

**R代码**:
```r
library(survival)
library(cmprsk)

# 数据准备
time <- c(10, 15, 20, 25, 30)
status <- c(1, 2, 1, 0, 2)  # 1=感兴趣事件, 2=竞争事件, 0=删失
group <- c(0, 1, 0, 1, 0)

# 累积发生率估计
cif <- cuminc(time, status, group)
plot(cif)

# Fine-Gray模型
crr(time, status, cov1=matrix(group), failcode=1)
```

### 2. 多状态模型 (Multi-State Models)

**用途**: 建模疾病在多个状态间的转移

**常见结构**:
- **疾病进展模型**: 正常 → 早期 → 晚期 → 死亡
- **治疗反应模型**: 治疗 → 缓解 → 复发 → 死亡
- **移植模型**: 等待 → 移植 → 排斥/存活

**转移类型**:
- **渐进式**: 只能向前转移 (如疾病分期)
- **可逆式**: 可来回转移 (如病情反复)
- **竞争式**: 多个可能终点

**分析方法**:
- **Aalen-Johansen估计**: 非参数转移概率
- **Cox型多状态模型**: 半参数方法
- **马尔可夫模型**: 假设无记忆性

**眼科应用**:
- AMD进展: 早期AMD → 中期 → 晚期 → 失明
- 糖尿病视网膜病变分级进展
- 青光眼视野缺损进展

**R代码**:
```r
library(mstate)

# 定义转移矩阵
tmat <- transMat(x = list(c(2, 3), c(3), c()), 
                 names = c("Normal", "Early", "Late"))

# 数据准备
msdata <- msprep(time = c(NA, "time1", "time2"),
                 status = c(NA, "status1", "status2"),
                 data = df,
                 trans = tmat)

# Cox模型
coxph(Surv(Tstart, Tstop, status) ~ covariate + strata(trans),
      data = msdata)
```

---

## 第二部分: 缺失数据处理

### 1. 缺失数据机制

**Rubin分类**:

| 机制 | 缩写 | 含义 | 示例 |
|------|------|------|------|
| 完全随机缺失 | MCAR | 缺失与任何数据无关 | 随机丢失样本 |
| 随机缺失 | MAR | 缺失与观测数据有关 | 视力差的患者更易失访 |
| 非随机缺失 | MNAR | 缺失与未观测值有关 | 病情重的患者不愿复查 |

### 2. 多重插补 (Multiple Imputation)

**原理**: 创建m个完整数据集→分别分析→合并结果

**步骤**:
1. **插补**: 基于观测数据生成m个完整数据集
2. **分析**: 对每个数据集进行分析
3. **合并**: 使用Rubin规则合并结果

**Rubin规则**:
```
总方差 = 组内方差 + (1+1/m) × 组间方差
```

**插补方法**:
- **MICE**: 链式方程多重插补 (最常用)
- **EM算法**: 期望最大化
- **贝叶斯方法**: 数据增强

**R代码**:
```r
library(mice)

# 查看缺失模式
md.pattern(df)

# 多重插补
imp <- mice(df, m=5, method="pmm", seed=123)

# 查看插补值
complete(imp, 1)  # 第1个完整数据集

# 在每个数据集上分析
fit <- with(imp, lm(BCVA ~ treatment + age))

# 合并结果
pooled <- pool(fit)
summary(pooled)
```

### 3. 逆概率加权 (Inverse Probability Weighting, IPW)

**原理**: 对完整观测赋予更高权重，使其代表总体

**权重计算**:
```
权重 = 1 / P(观测完整 | 协变量)
```

**边际结构模型 (MSM)**:
- 结合IPW和GEE
- 处理时变混杂

**R代码**:
```r
library(ipw)

# 估计缺失概率
missing_model <- glm(missing ~ age + baseline_VA, 
                     family=binomial, data=df)

# 计算权重
df$weight <- 1 / predict(missing_model, type="response")

# 加权分析
library(survey)
design <- svydesign(ids=~1, weights=~weight, data=df)
svymean(~BCVA, design)
```

---

## 第三部分: 网络Meta分析

### 1. 基本概念

**定义**: 同时比较多种治疗方法，包括直接和间接比较

**优势**:
- 可以比较未直接对比的治疗
- 提高统计效能
- 对治疗方法进行排序

**网络结构**:
- **节点**: 治疗方法
- **边**: 直接比较的研究

### 2. 一致性假设

**直接比较 vs 间接比较**:
- 一致性: 直接和间接估计一致
- 不一致性: 存在差异，需调查原因

**检测方法**:
- **节点劈裂法 (Node-Splitting)**: 比较直接和间接估计
- **设计-治疗交互模型**: 检验设计效应
- **Q统计量**: 评估整体不一致性

### 3. 统计模型

**频率学派方法**:
- **一致性模型**: 假设所有证据一致
- **不一致性模型**: 允许直接/间接估计不同

**贝叶斯方法**:
- 使用MCMC估计
- 可计算治疗排序概率
- 更灵活处理复杂网络

### 4. 结果解释

**SUCRA值**:
- Surface Under the Cumulative Ranking curve
- 0-100%，越高越好
- 用于治疗方法排序

**R代码**:
```r
library(netmeta)

# 网络Meta分析
net <- netmeta(TE, seTE, treat1, treat2, studlab,
               data=df, sm="OR")

# 结果
summary(net)

# 网络图
netgraph(net)

# 森林图
forest(net)

# 一致性检验
netsplit(net)

# 排序
netrank(net)

# SUCRA
netleague(net)
```

---

## 第四部分: 可解释AI

### 1. SHAP值 (SHapley Additive exPlanations)

**原理**: 基于博弈论Shapley值，计算每个特征对预测的贡献

**特性**:
- **可加性**: 各特征贡献之和等于预测值与基准值之差
- **一致性**: 特征重要性随模型变化单调变化
- **局部准确性**: 针对单个预测的解释

**计算**:
```
φ_i = Σ_S⊆N\{i} [|S|!(|N|-|S|-1)! / |N|!] × [f(S∪{i}) - f(S)]
```

**眼科应用**:
- 解释AI诊断模型的决策依据
- 识别对预测最重要的影像特征
- 个体化解释每个病例

**Python代码**:
```python
import shap
import xgboost as xgb

# 训练模型
model = xgb.XGBClassifier()
model.fit(X_train, y_train)

# SHAP解释器
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)

# 全局特征重要性
shap.summary_plot(shap_values, X_test)

# 单个预测解释
shap.force_plot(explainer.expected_value, 
                shap_values[0], X_test.iloc[0])

# 依赖图
shap.dependence_plot("feature_name", shap_values, X_test)
```

### 2. LIME (Local Interpretable Model-agnostic Explanations)

**原理**: 在单个预测附近拟合可解释的简单模型

**步骤**:
1. 在待解释样本附近扰动生成新样本
2. 用复杂模型预测新样本
3. 根据距离加权拟合简单模型
4. 解释简单模型

**特点**:
- 模型无关: 适用于任何黑盒模型
- 局部解释: 针对单个预测
- 可解释性: 使用线性模型或决策树

**Python代码**:
```python
import lime
from lime.lime_tabular import LimeTabularExplainer

# 创建解释器
explainer = LimeTabularExplainer(
    X_train.values,
    feature_names=X_train.columns,
    class_names=['Normal', 'Disease'],
    mode='classification'
)

# 解释单个预测
exp = explainer.explain_instance(
    X_test.iloc[0].values, 
    model.predict_proba,
    num_features=5
)

# 显示解释
exp.show_in_notebook()
exp.as_list()
```

### 3. 注意力可视化

**用途**: 解释深度学习模型的关注区域

**方法**:
- **Grad-CAM**: 梯度加权类激活映射
- **Attention Map**: Transformer的注意力权重
- **Occlusion Sensitivity**: 遮挡敏感性分析

**眼科应用**:
- 显示AI诊断眼底图像时关注的区域
- 验证模型是否关注正确的病理特征
- 发现新的生物标志物

**Python代码 (Grad-CAM)**:
```python
import torch
import cv2
import numpy as np

class GradCAM:
    def __init__(self, model, target_layer):
        self.model = model
        self.target_layer = target_layer
        self.gradients = None
        self.activations = None
        
    def save_gradient(self, grad):
        self.gradients = grad
        
    def forward(self, x):
        self.activations = self.target_layer(x)
        self.activations.register_hook(self.save_gradient)
        return self.model(x)
    
    def generate(self, input_image, target_class):
        # 前向传播
        output = self.forward(input_image)
        
        # 反向传播
        self.model.zero_grad()
        output[0, target_class].backward()
        
        # 计算Grad-CAM
        pooled_gradients = torch.mean(self.gradients, dim=[0, 2, 3])
        activations = self.activations.detach()
        
        for i in range(activations.shape[1]):
            activations[:, i, :, :] *= pooled_gradients[i]
            
        heatmap = torch.mean(activations, dim=1).squeeze()
        heatmap = np.maximum(heatmap, 0)
        heatmap /= torch.max(heatmap)
        
        return heatmap.numpy()

# 使用
grad_cam = GradCAM(model, target_layer)
heatmap = grad_cam.generate(image, target_class=1)
```

---

## 第五部分: 临床试验设计

### 1. 适应性设计 (Adaptive Design)

**定义**: 根据期中分析结果调整试验设计的临床试验

**优势**:
- 提高试验效率
- 减少样本量
- 早期终止无效/优效治疗
- 适应性分配提高成功概率

**调整类型**:

| 调整类型 | 说明 | 示例 |
|---------|------|------|
| 样本量重估 | 根据观察到的效应量调整样本量 | 效应量小于预期时增加样本 |
| 治疗组淘汰 | 早期淘汰无效治疗组 | 多臂试验中淘汰无效药物 |
| 适应性随机化 | 根据疗效调整随机化比例 | 更多患者分配到有效组 |
| 人群富集 | 根据亚组分析调整入组标准 | 仅纳入获益人群 |
| 无缝II/III期 | 合并II期和III期 | 减少总体开发时间 |

**统计考虑**:
- 控制整体I型错误率
- 避免操作偏倚
- 预先规定适应性规则

### 2. 平台试验 (Platform Trial)

**定义**: 持续评估多种治疗的适应性试验平台

**特点**:
- 多治疗、多疾病
- 共享对照组
- 治疗可动态加入/退出
- 使用主协议 (Master Protocol)

**类型**:
- **伞式试验**: 一种疾病，多种治疗
- **篮式试验**: 一种治疗，多种疾病
- **平台试验**: 多疾病、多治疗

**眼科应用**:
- 评估多种抗VEGF药物治疗DME
- 比较不同手术方式治疗白内障
- 评估多种基因治疗遗传性视网膜病变

### 3. 主协议设计 (Master Protocol)

**组成**:
1. **核心设计**: 统一的入组标准、终点、统计方法
2. **治疗特定附录**: 各治疗的剂量、安全性监测
3. **适应性规则**: 预先规定的调整机制

**统计方法**:
- **贝叶斯自适应随机化**: 根据累积数据调整随机化比例
- **预测概率**: 计算试验成功的预测概率
- **响应适应性**: 根据患者响应调整治疗分配

**R代码 (模拟适应性设计)**:
```r
# 简单的适应性随机化模拟
adaptive_randomization <- function(n_patients, n_arms, 
                                   true_effects, 
                                   adapt_after = 50) {
  allocations <- rep(0, n_arms)
  responses <- vector("list", n_arms)
  
  for(i in 1:n_patients) {
    if(i <= adapt_after) {
      # 初始等概率随机化
      arm <- sample(1:n_arms, 1)
    } else {
      # 适应性随机化 (Thompson Sampling简化版)
      success_rates <- sapply(responses, function(x) {
        if(length(x) == 0) return(0.5)
        mean(x)
      })
      probs <- success_rates / sum(success_rates)
      arm <- sample(1:n_arms, 1, prob=probs)
    }
    
    # 模拟响应
    response <- rbinom(1, 1, true_effects[arm])
    allocations[arm] <- allocations[arm] + 1
    responses[[arm]] <- c(responses[[arm]], response)
  }
  
  list(allocations=allocations, responses=responses)
}

# 运行模拟
result <- adaptive_randomization(
  n_patients = 200,
  n_arms = 3,
  true_effects = c(0.3, 0.5, 0.6),
  adapt_after = 50
)

print(result$allocations)
```

---

## 第六部分: 方法选择决策树

### 1. 生存数据复杂情况

```
单一终点?
├── 是 → 标准Cox模型
└── 否 (竞争风险)
    ├── 只关心发生率 → 累积发生率函数
    └── 关心风险因素 → Fine-Gray模型

多状态?
├── 是 → 多状态模型
└── 否 → 标准方法
```

### 2. 缺失数据

```
缺失机制?
├── MCAR → 完整案例分析 (可能损失效能)
├── MAR → 多重插补
└── MNAR → 敏感性分析 + 模式混合模型

缺失比例?
├── <5% → 完整案例分析
├── 5-40% → 多重插补
└── >40% → 考虑研究设计问题
```

### 3. Meta分析

```
多种治疗比较?
├── 是 → 网络Meta分析
└── 否 → 标准Meta分析

是否有直接比较?
├── 有 → 一致性检验
└── 无 → 仅间接证据 (谨慎解释)
```

---

## 学习资源

### Wikipedia页面
- 竞争风险: https://en.wikipedia.org/wiki/Competing_risks
- 多重插补: https://en.wikipedia.org/wiki/Multiple_imputation
- 逆概率加权: https://en.wikipedia.org/wiki/Inverse_probability_weighting
- 适应性试验: https://en.wikipedia.org/wiki/Adaptive_clinical_trial
- Shapley值: https://en.wikipedia.org/wiki/Shapley_value

---

*最后更新: 2026-03-07*
