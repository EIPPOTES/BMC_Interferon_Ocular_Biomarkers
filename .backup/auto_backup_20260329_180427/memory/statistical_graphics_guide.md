# 统计学图表制作指南

*学习来源: Wikipedia | 学习时间: 2026-03-07*
*针对: 医学科研数据可视化*

---

## 第一部分: 数据可视化基础

### 1. 数据可视化的目的

**核心目标**:
- 帮助理解复杂数据
- 发现数据中的模式、趋势和异常
- 支持决策和假设生成
- 有效传达研究发现

**医学研究中的重要性**:
- 直观展示治疗效果
- 识别疾病进展模式
- 比较不同治疗方案
- 呈现患者特征分布

### 2. 图表选择决策树

```
数据类型?
├── 分类数据
│   ├── 比较各类别大小 → 条形图/柱状图
│   ├── 显示构成比例 → 饼图/堆叠条形图
│   └── 显示趋势 → 分组条形图
│
├── 连续数据
│   ├── 单变量分布 → 直方图/箱线图/密度图
│   ├── 两变量关系 → 散点图
│   ├── 时间趋势 → 折线图
│   └── 多组比较 → 箱线图/小提琴图
│
└── 生存数据
    ├── 生存曲线 → Kaplan-Meier曲线
    └── 多研究比较 → 森林图
```

---

## 第二部分: 基础统计图表

### 1. 直方图 (Histogram)

**用途**: 展示连续变量的分布

**关键要素**:
- **分箱 (Binning)**: 将数据范围分成若干区间
- **频数/频率**: 每个区间的观测数
- **面积意义**: 总面积代表总概率(密度直方图)

**与条形图的区别**:

| 特征 | 直方图 | 条形图 |
|------|--------|--------|
| 数据类型 | 连续变量 | 分类变量 |
| 柱子间距 | 无间隙 | 有间隙 |
| 柱子顺序 | 按数值排序 | 可任意排序 |
| 面积意义 | 代表频率 | 仅代表高度 |

**最佳实践**:
- 选择合适的分箱数 (Sturges公式: k = 1 + 3.322 log n)
- 标注样本量
- 添加密度曲线辅助观察

**R代码**:
```r
# 基础直方图
hist(data$BCVA, 
     main="Distribution of BCVA",
     xlab="BCVA (LogMAR)",
     ylab="Frequency",
     col="lightblue",
     border="white",
     breaks=20)

# 添加密度曲线
hist(data$BCVA, freq=FALSE, col="lightblue")
lines(density(data$BCVA), col="red", lwd=2)

# ggplot2版本
library(ggplot2)
ggplot(data, aes(x=BCVA)) +
  geom_histogram(aes(y=..density..), 
                 bins=30, fill="lightblue", color="white") +
  geom_density(color="red", size=1) +
  labs(title="Distribution of BCVA",
       x="BCVA (LogMAR)", y="Density") +
  theme_minimal()
```

**眼科应用**:
- 视力分布
- 眼压分布
- OCT厚度分布

---

### 2. 箱线图 (Box Plot)

**用途**: 展示数据分布和组间比较

**组成要素**:
```
    │                    ↑ 上 whisker (Q3 + 1.5×IQR)
    ├─●─ 异常值
    │
   ┌┴┐ 上边缘 (Q3)
   │ │
   ├┼┤ 中位数线
   │ │
   └┬┘ 下边缘 (Q1)
    │
    ├─ 下 whisker (Q1 - 1.5×IQR)
```

**五数概括**:
- 最小值 (或下whisker)
- 第一四分位数 (Q1, 25%)
- 中位数 (Q2, 50%)
- 第三四分位数 (Q3, 75%)
- 最大值 (或上whisker)

**优点**:
- 紧凑展示分布特征
- 便于比较多组数据
- 识别异常值
- 不假设分布形态

**R代码**:
```r
# 基础箱线图
boxplot(BCVA ~ treatment, data=df,
        main="BCVA by Treatment Group",
        xlab="Treatment",
        ylab="BCVA (LogMAR)",
        col=c("lightblue", "lightgreen"))

# ggplot2版本
ggplot(df, aes(x=treatment, y=BCVA, fill=treatment)) +
  geom_boxplot() +
  geom_jitter(width=0.2, alpha=0.3) +  # 添加原始点
  labs(title="BCVA by Treatment Group",
       x="Treatment", y="BCVA (LogMAR)") +
  theme_minimal() +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))

# 小提琴图 (箱线图升级版)
ggplot(df, aes(x=treatment, y=BCVA, fill=treatment)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.2) +
  theme_minimal()
```

**眼科应用**:
- 比较不同治疗组的视力
- 比较不同疾病阶段的OCT参数
- 识别异常值患者

---

### 3. 散点图 (Scatter Plot)

**用途**: 展示两个连续变量的关系

**变体**:
- **基础散点图**: 显示相关关系
- **带回归线**: 显示趋势
- **分组散点图**: 按类别着色
- **气泡图**: 点大小代表第三变量

**解读要点**:
- **方向**: 正相关/负相关
- **形式**: 线性/非线性
- **强度**: 点的分散程度
- **异常值**: 远离主体的点

**R代码**:
```r
# 基础散点图
plot(df$baseline_BCVA, df$followup_BCVA,
     main="Baseline vs Follow-up BCVA",
     xlab="Baseline BCVA",
     ylab="Follow-up BCVA",
     pch=19, col="blue")
abline(0, 1, col="red", lty=2)  # 添加对角线

# ggplot2版本
ggplot(df, aes(x=baseline_BCVA, y=followup_BCVA)) +
  geom_point(alpha=0.5, color="blue") +
  geom_smooth(method="lm", color="red", se=TRUE) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(title="Baseline vs Follow-up BCVA",
       subtitle="Red line: regression, Dashed: identity") +
  theme_minimal()

# 分组散点图
ggplot(df, aes(x=CST, y=BCVA, color=treatment)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~timepoint) +
  theme_minimal()
```

**眼科应用**:
- 视力与OCT厚度关系
- 基线与随访视力相关性
- 眼压与视野缺损关系

---

### 4. 条形图/柱状图 (Bar Chart)

**用途**: 比较不同类别的数值

**类型**:
- **垂直条形图**: 类别较少时
- **水平条形图**: 类别较多时
- **分组条形图**: 多组比较
- **堆叠条形图**: 显示构成

**注意事项**:
- Y轴必须从0开始 (避免误导)
- 按数值大小排序 (除非有自然顺序)
- 添加数值标签

**R代码**:
```r
# 基础条形图
barplot(table(df$treatment),
        main="Patient Distribution by Treatment",
        xlab="Treatment",
        ylab="Count",
        col="steelblue")

# ggplot2版本
ggplot(df_summary, aes(x=treatment, y=mean_BCVA, fill=treatment)) +
  geom_bar(stat="identity", width=0.7) +
  geom_errorbar(aes(ymin=mean_BCVA-sd_BCVA, 
                    ymax=mean_BCVA+sd_BCVA),
                width=0.2) +
  labs(title="Mean BCVA by Treatment",
       y="BCVA (LogMAR)") +
  theme_minimal()

# 分组条形图
ggplot(df_summary, aes(x=timepoint, y=mean_BCVA, fill=treatment)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=mean_BCVA-se, ymax=mean_BCVA+se),
                position=position_dodge(0.7), width=0.2) +
  labs(title="BCVA by Treatment and Time",
       x="Time Point", y="BCVA (LogMAR)") +
  theme_minimal()
```

**眼科应用**:
- 不同治疗组的患者数
- 各疾病分期的比例
- 不良事件发生率比较

---

## 第三部分: 医学研究专用图表

### 1. Kaplan-Meier生存曲线

**用途**: 展示时间-事件数据

**关键要素**:
- **X轴**: 时间
- **Y轴**: 生存概率 (0-1)
- **阶梯线**: 事件发生时的下降
- **删失标记**: 小竖线
- **置信区间**: 阴影区域

**解读**:
- 曲线位置越高越好
- 中位生存时间
- 特定时间点的生存率

**R代码**:
```r
library(survival)
library(survminer)

# 创建生存对象
surv_obj <- Surv(time=df$time_to_event, event=df$event)

# 拟合Kaplan-Meier
km_fit <- survfit(surv_obj ~ treatment, data=df)

# 基础绘图
plot(km_fit, main="Survival by Treatment",
     xlab="Time (months)", ylab="Survival Probability",
     col=c("red", "blue"), lwd=2)
legend("topright", legend=c("Treatment A", "Treatment B"),
       col=c("red", "blue"), lwd=2)

# ggsurvplot (推荐)
ggsurvplot(km_fit,
           pval=TRUE,              # 显示log-rank检验P值
           conf.int=TRUE,          # 显示置信区间
           risk.table=TRUE,        # 显示风险表
           risk.table.col="strata",
           linetype=1,
           surv.median.line="hv",  # 显示中位生存时间
           ggtheme=theme_minimal(),
           palette=c("#E7B800", "#2E9FDF"),
           title="Time to Vision Loss",
           xlab="Time (months)",
           ylab="Vision Preservation Probability")
```

**眼科应用**:
- 视力维持时间
- 疾病进展时间
- 治疗失败时间

---

### 2. 森林图 (Forest Plot)

**用途**: Meta分析或多组比较

**组成**:
- **左侧**: 研究/组别信息
- **中间**: 效应量点估计和置信区间
- **右侧**: 数值
- **底部**: 汇总效应量

**解读**:
- 菱形代表汇总效应
- 垂线代表无效线 (OR=1或HR=1)
- 点偏左=保护因素，偏右=危险因素

**R代码**:
```r
library(meta)
library(forestplot)

# 使用meta包
meta_result <- metagen(TE=effect_size, seTE=se, 
                       studlab=study_name, data=df)

# 绘制森林图
forest(meta_result,
       layout="JAMA",           # JAMA风格
       prediction=TRUE,         # 预测区间
       print.tau2=TRUE,         # 显示异质性
       leftcols=c("studlab", "TE", "seTE"),
       rightcols=c("effect", "ci"),
       xlab="Effect Size")

# 使用forestplot包 (更灵活)
library(forestplot)

table_text <- cbind(
  c("Study", "Study 1", "Study 2", "Study 3", "Overall"),
  c("N", "100", "150", "120", "370"),
  c("Effect", "0.5", "0.8", "0.3", "0.55")
)

coef <- c(NA, 0.5, 0.8, 0.3, 0.55)
low <- c(NA, 0.2, 0.5, 0.0, 0.35)
high <- c(NA, 0.8, 1.1, 0.6, 0.75)

forestplot(table_text, coef, low, high,
           zero=1, boxsize=0.3,
           col=fpColors(box="royalblue", line="darkblue"),
           xlab="Odds Ratio")
```

**眼科应用**:
- Meta分析结果展示
- 亚组分析结果
- 多中心研究各中心结果

---

### 3. 热图 (Heatmap)

**用途**: 展示矩阵数据的相关性或聚类

**特点**:
- 颜色深浅代表数值大小
- 可同时展示行和列的聚类
- 适合展示相关性矩阵

**R代码**:
```r
library(pheatmap)
library(RColorBrewer)

# 相关性矩阵热图
cor_matrix <- cor(df[, c("BCVA", "CST", "Age", "IOP")])

pheatmap(cor_matrix,
         color=colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers=TRUE,
         number_color="black",
         main="Correlation Matrix")

# 复杂热图 (ComplexHeatmap)
library(ComplexHeatmap)

Heatmap(cor_matrix,
        name="Correlation",
        col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cell_fun=function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cor_matrix[i, j]), 
                    x, y, gp=gpar(fontsize=10))
        })
```

**眼科应用**:
- 多参数相关性分析
- 基因表达数据
- 患者聚类分析

---

## 第四部分: 高级可视化

### 1. 小提琴图 (Violin Plot)

**用途**: 箱线图的升级版，展示完整分布

**优点**:
- 展示分布形状
- 保留箱线图的统计信息
- 适合大样本数据

**R代码**:
```r
ggplot(df, aes(x=treatment, y=BCVA, fill=treatment)) +
  geom_violin(alpha=0.5, trim=FALSE) +
  geom_boxplot(width=0.2, fill="white") +
  stat_summary(fun=mean, geom="point", 
               shape=23, size=3, fill="red") +
  labs(title="BCVA Distribution by Treatment",
       y="BCVA (LogMAR)") +
  theme_minimal()
```

---

### 2. 配对图 (Pairs Plot)

**用途**: 展示多变量间的关系

**R代码**:
```r
library(GGally)

# 基础配对图
pairs(df[, c("BCVA", "CST", "Age", "IOP")])

# ggplot2版本
ggpairs(df, columns=c("BCVA", "CST", "Age", "IOP"),
        aes(color=treatment),
        upper=list(continuous=wrap("cor", size=4)),
        lower=list(continuous=wrap("points", alpha=0.3)),
        diag=list(continuous=wrap("densityDiag", alpha=0.5))) +
  theme_minimal()
```

---

## 第五部分: 图表设计原则

### 1. 颜色使用

**原则**:
- 色盲友好 (避免红绿组合)
- 黑白打印友好
- 一致的颜色编码

**推荐配色方案**:
```r
# 色盲友好配色
cb_palette <- c("#999999", "#E69F00", "#56B4E9", 
                "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7")

# 使用
scale_fill_manual(values=cb_palette)
```

### 2. 字体和标签

**最佳实践**:
- 标题: 14-16pt
- 轴标签: 12pt
- 刻度标签: 10pt
- 避免倾斜标签
- 使用有意义的标题

### 3. 图表导出

**期刊要求**:
- 分辨率: 300-600 DPI
- 格式: TIFF (位图), EPS/PDF (矢量)
- 颜色模式: CMYK (印刷) 或 RGB (在线)

**R代码**:
```r
# 高分辨率导出
ggsave("figure1.tiff", plot=p, 
       width=8, height=6, 
       dpi=300, 
       compression="lzw")

# PDF矢量图
ggsave("figure1.pdf", plot=p,
       width=8, height=6)
```

---

## 第六部分: 眼科研究专用图表

### 1. 视力变化瀑布图 (Waterfall Plot)

**用途**: 展示每个患者的视力变化

**R代码**:
```r
ggplot(df, aes(x=reorder(patient_id, BCVA_change), 
               y=BCVA_change, 
               fill=BCVA_change>0)) +
  geom_bar(stat="identity", width=0.8) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_fill_manual(values=c("red", "green"),
                    labels=c("Loss", "Gain")) +
  labs(title="BCVA Change by Patient",
       x="Patient", y="BCVA Change (letters)") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
```

### 2. OCT厚度地形图

**用途**: 展示视网膜厚度分布

**R代码**:
```r
library(ggplot2)
library(reshape2)

# 假设有ETDRS网格数据
oct_matrix <- matrix(oct_data, nrow=9, ncol=9)
oct_melted <- melt(oct_matrix)

ggplot(oct_melted, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red",
                      name="Thickness (μm)") +
  coord_fixed() +
  labs(title="Retinal Thickness Map",
       x="", y="") +
  theme_minimal()
```

### 3. 视野缺损进展图

**用途**: 展示视野随时间的变化

**R代码**:
```r
ggplot(df, aes(x=time, y=MD, group=patient_id)) +
  geom_line(alpha=0.3) +
  geom_smooth(aes(group=1), method="lm", 
              color="red", size=1.5) +
  facet_wrap(~treatment) +
  labs(title="Visual Field Progression",
       x="Time (years)", y="Mean Deviation (dB)") +
  theme_minimal()
```

---

## 第七部分: 常用R包总结

| 包名 | 用途 | 推荐场景 |
|------|------|---------|
| **ggplot2** | 通用绘图 | 所有基础图表 |
| **survminer** | 生存分析 | Kaplan-Meier曲线 |
| **meta** | Meta分析 | 森林图 |
| **pheatmap** | 热图 | 相关性矩阵 |
| **GGally** | 配对图 | 多变量探索 |
| **plotly** | 交互图 | 在线展示 |
| **patchwork** | 拼图 | 多图组合 |

---

## 学习资源

- Data Visualization: https://en.wikipedia.org/wiki/Data_visualization
- Histogram: https://en.wikipedia.org/wiki/Histogram
- Box Plot: https://en.wikipedia.org/wiki/Box_plot
- Scatter Plot: https://en.wikipedia.org/wiki/Scatter_plot
- Forest Plot: https://en.wikipedia.org/wiki/Forest_plot

---

*最后更新: 2026-03-07*
