#!/usr/bin/env python3
"""
分析最终版表格数据，为制作SCI图表做准备
"""

import pandas as pd
import numpy as np
import os
import json
from datetime import datetime

print("=" * 80)
print("分析表格数据用于制作SCI图表")
print(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)

# 路径设置
final_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版"
tables_dir = os.path.join(final_dir, "03_Tables")
manuscript_dir = os.path.join(final_dir, "01_Manuscript")

# 1. 读取关键表格数据
print(f"\n📊 读取关键表格数据...")

# 描述性统计
desc_file = os.path.join(tables_dir, "Descriptive_Statistics_463eyes_20260315.xlsx")
if os.path.exists(desc_file):
    df_desc = pd.read_excel(desc_file)
    print(f"描述性统计: {df_desc.shape[0]}行 × {df_desc.shape[1]}列")
    print(f"列名: {df_desc.columns.tolist()}")
    print(f"前5个参数: {df_desc['Parameter'].head().tolist()}")
else:
    print(f"⚠️ 文件不存在: {desc_file}")

# 组间比较
group_file = os.path.join(tables_dir, "Group_Comparison_463eyes_20260315.xlsx")
if os.path.exists(group_file):
    df_group = pd.read_excel(group_file)
    print(f"\n组间比较: {df_group.shape[0]}行 × {df_group.shape[1]}列")
    print(f"列名: {df_group.columns.tolist()}")
    
    # 查看显著指标
    significant = df_group[df_group['P_value'] < 0.05]
    print(f"显著指标数 (P<0.05): {len(significant)}/{len(df_group)}")
    
    # 按效应量排序
    top_effects = df_group.sort_values('Cohens_d').head(10)
    print(f"效应量最大(负向)的10个指标:")
    for idx, row in top_effects.iterrows():
        print(f"  {row['Parameter']}: Cohen's d={row['Cohens_d']:.3f}, P={row['P_value']:.6f}")
else:
    print(f"⚠️ 文件不存在: {group_file}")

# ROC分析
roc_file = os.path.join(tables_dir, "ROC_Analysis_463eyes_20260315.xlsx")
if os.path.exists(roc_file):
    df_roc = pd.read_excel(roc_file)
    print(f"\nROC分析: {df_roc.shape[0]}行 × {df_roc.shape[1]}列")
    print(f"列名: {df_roc.columns.tolist()}")
    print(f"前5个指标: {df_roc.head().to_string()}")
else:
    print(f"⚠️ 文件不存在: {roc_file}")

# 机器学习模型性能
ml_file = os.path.join(tables_dir, "机器学习模型性能对比_20260315.xlsx")
if os.path.exists(ml_file):
    df_ml = pd.read_excel(ml_file)
    print(f"\n机器学习模型性能: {df_ml.shape[0]}行 × {df_ml.shape[1]}列")
    print(f"模型列表: {df_ml['模型'].tolist()}")
    print(f"最佳AUC: {df_ml['AUC'].max():.3f} ({df_ml.loc[df_ml['AUC'].idxmax(), '模型']})")
else:
    print(f"⚠️ 文件不存在: {ml_file}")

# 特征重要性
feat_file = os.path.join(tables_dir, "特征重要性分析_20260315.xlsx")
if os.path.exists(feat_file):
    df_feat = pd.read_excel(feat_file)
    print(f"\n特征重要性: {df_feat.shape[0]}行 × {df_feat.shape[1]}列")
    print(f"列名: {df_feat.columns.tolist()}")
    
    # 查看最重要的特征
    if 'importance' in df_feat.columns or '特征重要性' in df_feat.columns:
        importance_col = 'importance' if 'importance' in df_feat.columns else '特征重要性'
        top_features = df_feat.sort_values(importance_col, ascending=False).head(10)
        print(f"最重要的10个特征:")
        for idx, row in top_features.iterrows():
            print(f"  {row['特征'] if '特征' in df_feat.columns else 'Feature'}: {row[importance_col]:.4f}")
else:
    print(f"⚠️ 文件不存在: {feat_file}")

# 复合指标权重
weight_file = os.path.join(tables_dir, "复合指标权重_20260315.xlsx")
if os.path.exists(weight_file):
    df_weight = pd.read_excel(weight_file)
    print(f"\n复合指标权重: {df_weight.shape[0]}行 × {df_weight.shape[1]}列")
    print(f"列名: {df_weight.columns.tolist()}")
    
    # 查看权重最大的特征
    if '权重' in df_weight.columns:
        top_weights = df_weight.sort_values('权重', key=abs, ascending=False).head(10)
        print(f"权重绝对值最大的10个特征:")
        for idx, row in top_weights.iterrows():
            print(f"  {row['特征']}: {row['权重']:.4f}")
else:
    print(f"⚠️ 文件不存在: {weight_file}")

# 2. 分析需要制作的图表类型
print(f"\n🎨 分析需要制作的图表类型...")

# 基于SCI论文标准，建议制作以下图表
figures_needed = {
    "Figure 1": {
        "title": "Study flowchart",
        "type": "flowchart",
        "content": "Participant recruitment, inclusion/exclusion criteria, final sample size",
        "data_needed": "Sample sizes at each stage",
        "priority": "High"
    },
    "Figure 2": {
        "title": "Group comparison of key OCT parameters",
        "type": "forest plot / bar plot",
        "content": "Effect sizes (Cohen's d) with 95% confidence intervals for top 10-15 OCT parameters",
        "data_needed": "Group comparison results with effect sizes and confidence intervals",
        "priority": "High"
    },
    "Figure 3": {
        "title": "ROC curves for diagnostic performance",
        "type": "ROC curves",
        "content": "Traditional ROC curves for single parameters vs. machine learning composite score",
        "data_needed": "ROC data for single parameters and composite score",
        "priority": "High"
    },
    "Figure 4": {
        "title": "Feature importance and composite score weights",
        "type": "bar plot / waterfall plot",
        "content": "Machine learning feature importance and logistic regression weights for composite score",
        "data_needed": "Feature importance and composite score weights",
        "priority": "High"
    },
    "Figure 5": {
        "title": "Correlation between OCT parameters and depression severity",
        "type": "scatter plot / correlation matrix",
        "content": "Correlation of key OCT parameters with PHQ-9 scores",
        "data_needed": "Correlation analysis results",
        "priority": "Medium"
    },
    "Figure 6": {
        "title": "Subgroup analysis by sex and age",
        "type": "forest plot",
        "content": "Effect sizes stratified by sex and age groups",
        "data_needed": "Subgroup analysis results",
        "priority": "Medium"
    },
    "Supplementary Figure 1": {
        "title": "Machine learning model performance comparison",
        "type": "bar plot",
        "content": "Performance metrics (AUC, sensitivity, specificity) for 6 machine learning models",
        "data_needed": "Machine learning model performance data",
        "priority": "Medium"
    },
    "Supplementary Figure 2": {
        "title": "Simplified clinical scoring system",
        "type": "ROC curves / calibration plot",
        "content": "Performance of simplified 5-parameter scoring system",
        "data_needed": "Simplified scoring system performance",
        "priority": "Low"
    }
}

print(f"建议制作 {len(figures_needed)} 个图表:")
for fig_name, fig_info in figures_needed.items():
    print(f"  {fig_name}: {fig_info['title']} ({fig_info['priority']} priority)")

# 3. 检查数据完整性
print(f"\n🔍 检查数据完整性...")

data_availability = {
    "group_comparison": os.path.exists(group_file),
    "roc_analysis": os.path.exists(roc_file),
    "ml_performance": os.path.exists(ml_file),
    "feature_importance": os.path.exists(feat_file),
    "composite_weights": os.path.exists(weight_file),
    "correlation": os.path.exists(os.path.join(tables_dir, "相关性分析_OCT_vs_PHQ9_20260315.xlsx")),
    "subgroup": os.path.exists(os.path.join(tables_dir, "亚组分析结果_20260315.xlsx")),
    "multivariate": os.path.exists(os.path.join(tables_dir, "多变量回归_线性模型结果_20260315.xlsx"))
}

print("数据可用性:")
for data_type, available in data_availability.items():
    status = "✅" if available else "❌"
    print(f"  {status} {data_type}")

# 4. 生成图表制作计划
print(f"\n📋 生成图表制作计划...")

# 读取原始数据以获取更多信息
data_file = os.path.join(final_dir, "04_Data", "data_499eyes_20260315.xlsx")
if os.path.exists(data_file):
    print(f"原始数据文件: {data_file}")
    # 只读取列名以了解数据结构
    try:
        df_data_cols = pd.read_excel(data_file, nrows=0)
        print(f"数据列数: {len(df_data_cols.columns)}")
        
        # 统计OCT相关列
        oct_cols = [col for col in df_data_cols.columns if any(keyword in col for keyword in ['Retina', 'RNFL', 'GCL', 'Choroid', 'Disc', 'Cup', 'Rim'])]
        print(f"OCT相关列数: {len(oct_cols)}")
        
        # 查看前10个OCT列
        print(f"前10个OCT指标: {oct_cols[:10]}")
    except Exception as e:
        print(f"读取数据文件错误: {e}")
else:
    print(f"⚠️ 原始数据文件不存在: {data_file}")

# 5. 生成图表设计规范
print(f"\n🎯 生成图表设计规范...")

sci_figure_specs = {
    "尺寸": "根据期刊要求，通常:",
    "单栏宽度": "~8.5 cm (3.35 inches)",
    "双栏宽度": "~17.5 cm (6.89 inches)",
    "分辨率": "300 DPI (高分辨率)",
    "格式": "TIFF或EPS (矢量图更佳)",
    "字体": "Arial或Helvetica, 大小8-12 pt",
    "颜色": "CMYK用于印刷，RGB用于在线",
    "图例": "清晰明确，可独立理解",
    "误差线": "95%置信区间或标准差",
    "显著性标记": "*P<0.05, **P<0.01, ***P<0.001"
}

print("SCI图表规范:")
for key, value in sci_figure_specs.items():
    print(f"  {key}: {value}")

# 6. 生成制作脚本框架
print(f"\n🔧 生成制作脚本框架...")

script_template = """#!/usr/bin/env python3
\"\"\"
制作SCI论文图表 - {figure_title}
基于463眼OCT-MDD研究数据
\"\"\"

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
from scipy import stats
import matplotlib

# SCI图表设置
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1
})

# 颜色方案
colors = {
    'mdd': '#D55E00',  # 橙色
    'control': '#0072B2',  # 蓝色
    'significant': '#CC79A7',  # 粉色
    'ns': '#999999'  # 灰色
}

# 读取数据
def load_data():
    \"\"\"加载所需数据\"\"\"
    # 这里添加数据加载代码
    pass

def create_figure_{figure_num}():
    \"\"\"创建{figure_title}\"\"\"
    fig, ax = plt.subplots(figsize=(8.5/2.54, 8.5/2.54))  # 单栏宽度
    
    # 这里添加图表绘制代码
    
    # 设置标签
    ax.set_xlabel('X Label', fontsize=9)
    ax.set_ylabel('Y Label', fontsize=9)
    ax.set_title('{figure_title}', fontsize=10, fontweight='bold')
    
    # 添加图例
    ax.legend(loc='best', frameon=False)
    
    # 调整布局
    plt.tight_layout()
    
    return fig

def main():
    \"\"\"主函数\"\"\"
    print("制作{figure_title}...")
    
    # 创建图表
    fig = create_figure_{figure_num}()
    
    # 保存图表
    output_dir = "figures"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "Figure{figure_num}_{description}.tiff")
    fig.savefig(output_path, format='tiff', dpi=300)
    print(f"图表已保存: {{output_path}}")
    
    # 也保存为PDF格式（矢量图）
    pdf_path = output_path.replace('.tiff', '.pdf')
    fig.savefig(pdf_path, format='pdf')
    print(f"矢量图已保存: {{pdf_path}}")
    
    plt.close(fig)

if __name__ == "__main__":
    main()
"""

# 为每个图表生成脚本框架
for i, (fig_name, fig_info) in enumerate(figures_needed.items(), 1):
    figure_num = fig_name.split()[-1] if 'Figure' in fig_name else f"S{i}"
    description = fig_info['title'].lower().replace(' ', '_').replace('-', '_')
    
    script_content = script_template.format(
        figure_title=fig_info['title'],
        figure_num=figure_num,
        description=description
    )
    
    script_name = f"create_{description}.py"
    script_path = os.path.join("/root/.openclaw/workspace", script_name)
    
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(script_content)
    
    print(f"  ✅ 生成脚本框架: {script_name}")

# 7. 生成总体制作计划
print(f"\n📅 生成总体制作计划...")

production_plan = f"""# SCI图表制作计划
## OCT-MDD研究基于463眼版本
### 生成日期: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

---

## 📊 数据基础

### 样本信息
- **总样本**: 463眼 (有完整年龄性别数据)
- **MDD组**: 303眼
- **健康对照组**: 160眼
- **关键发现**: 最大效应量 Cohen's d=-0.498 (Retina_外环颞侧)

### 可用数据表格
1. `Descriptive_Statistics_463eyes_20260315.xlsx` - 描述性统计
2. `Group_Comparison_463eyes_20260315.xlsx` - 组间比较
3. `ROC_Analysis_463eyes_20260315.xlsx` - ROC分析
4. `机器学习模型性能对比_20260315.xlsx` - 机器学习性能
5. `特征重要性分析_20260315.xlsx` - 特征重要性
6. `复合指标权重_20260315.xlsx` - 复合指标权重
7. `相关性分析_OCT_vs_PHQ9_20260315.xlsx` - 相关性分析
8. `亚组分析结果_20260315.xlsx` - 亚组分析

---

## 🎨 图表制作清单

### 主要图表 (Figure 1-6)

| 编号 | 标题 | 类型 | 优先级 | 所需数据 | 预计完成时间 |
|------|------|------|--------|----------|--------------|
| **Figure 1** | Study flowchart | 流程图 | 高 | 样本筛选流程 | 30分钟 |
| **Figure 2** | Group comparison of key OCT parameters | 森林图 | 高 | 组间比较效应量 | 1小时 |
| **Figure 3** | ROC curves for diagnostic performance | ROC曲线 | 高 | ROC分析数据 | 1小时 |
| **Figure 4** | Feature importance and composite score weights | 条形图 | 高 | 特征重要性+权重 | 45分钟 |
| **Figure 5** | Correlation between OCT parameters and depression severity | 散点图 | 中 | 相关性数据 | 45分钟 |
| **Figure 6** | Subgroup analysis by sex and age | 森林图 | 中 | 亚组分析结果 | 45分钟 |

### 补充图表 (Supplementary Figures)

| 编号 | 标题 | 类型 | 优先级 | 所需数据 |
|------|------|------|--------|----------|
| **S1** | Machine learning model performance comparison | 条形图 | 中 | 机器学习性能 |
| **S2** | Simplified clinical scoring system | ROC曲线 | 低 | 简化评分性能 |
| **S3** | Distribution of key OCT parameters | 箱线图 | 低 | 描述性统计 |
| **S4** | Heatmap of correlations among OCT parameters | 热图 | 低 | 相关性矩阵 |

---

## 🔧 技术规范

### 图形设置
- **软件**: Python + Matplotlib + Seaborn
- **分辨率**: 300 DPI
- **格式**: TIFF (主), PDF/EPS (矢量备份)
- **字体**: Arial, 8-12 pt
- **颜色**: 色盲友好配色 (ColorBrewer Set2)

### 尺寸规范
- **单栏**: 8.5 cm (3.35 inches) 宽
- **1.5栏**: 12.7 cm (5.0 inches) 宽  
- **双栏**: 17.5 cm (6.89 inches) 宽
- **高度**: 根据内容调整，保持适当纵横比

### 学术标准
1. **误差线**: 95%置信区间或标准差
2. **显著性标记**: *P<0.05, **P<0.01, ***P<0.001
3. **图例**: 清晰独立，避免图表内过多文字
4. **坐标轴**: 明确标注单位和尺度
5. **颜色一致性**: 相同组别在不同图表中颜色一致

---

## 📈 图表内容详情

### Figure 1: 研究流程图
- 纳入标准: 首发未用药MDD患者，年龄18-60岁
- 排除标准: 眼科疾病、神经系统疾病、药物影响
- 最终样本: 463眼 (303 MDD, 160 对照)
- 数据完整性: 基于年龄性别数据完整性筛选

### Figure 2: 组间比较森林图
- 显示前15个效应量最大的OCT指标
- 包括: Cohen's d, 95% CI, P值
- 分组: 黄斑参数、RNFL参数、GCL参数、视盘参数
- 突出显示: Retina_外环颞侧 (d=-0.498)

### Figure 3: ROC曲线
- 传统ROC: 最佳单指标 (Cup-to-disc ratio, AUC=0.571)
- 机器学习: 随机森林 (AUC=0.730), 复合指标 (AUC=0.799)
- 对角线: 参考线 (AUC=0.5)
- 阴影: 95%置信区间

### Figure 4: 特征重要性
- 左侧: 随机森林特征重要性 (前20个特征)
- 右侧: 逻辑回归复合指标权重 (绝对值前20个)
- 颜色编码: 正权重(红色), 负权重(蓝色)
- 分组: RNFL、GCL、Retina、Choroid、Disc

### Figure 5: 相关性散点图
- X轴: PHQ-9评分 (0-27)
- Y轴: 关键OCT指标厚度
- 回归线: 线性拟合 + 95%置信带
- 分面: 按指标分组显示
- 标注: 相关系数(r)和P值

### Figure 6: 亚组分析森林图
- 分层: 性别(男性/女性)、年龄(年轻/年长)
- 指标: Retina_外环颞侧 (主要指标)
- 显示: 各组效应量、95% CI、P值
- 异质性检验: Q统计量和P值

---

## 🚀 实施步骤

### 第一阶段 (高优先级，今日完成)
1. **Figure 2**: 组间比较森林图
2. **Figure 3**: ROC曲线对比图
3. **Figure 4**: 特征重要性图

### 第二阶段 (中优先级，今日完成)
4. **Figure 1**: 研究流程图
5. **Figure 5**: 相关性散点图
6. **Figure 6**: 亚组分析森林图

### 第三阶段 (低优先级，可延后)
7. **Supplementary Figures**: 补充图表

### 质量控制
- 每个图表完成后检查: 分辨率、字体、颜色、标注
- 一致性检查: 不同图表间颜色和样式一致
- 学术规范: 符合目标期刊要求

---

## 📁 输出文件

### 图表文件
```
figures/
├── Figure1_Study_flowchart.tiff
├── Figure2_Group_comparison.tiff
├── Figure3_ROC_curves.tiff
├── Figure4_Feature_importance.tiff
├── Figure5_Correlation.tiff
├── Figure6_Subgroup_analysis.tiff
├── Supplementary_Figure1_ML_performance.tiff
└── Supplementary_Figure2_Simplified_scoring.tiff
```

### 矢量备份
```
figures/vector/
├── 所有图表的PDF版本
└── 所有图表的EPS版本 (如需要)
```

### 脚本文件
```
scripts/
├── create_figure1.py
├── create_figure2.py
├── create_figure3.py
├── create_figure4.py
├── create_figure5.py
└── create_figure6.py
```

---

## ⏰ 时间安排

### 今日 (2026-03-15)
- **17:40-18:30**: Figure 2, 3, 4 制作
- **18:30-19:00**: Figure 1, 5, 6 制作
- **19:00-19:30**: 质量控制和完善

### 后续
- 补充图表制作 (如时间允许)
- 图表整合到论文中
- 根据期刊要求调整格式

---

## 📞 注意事项

1. **数据准确性**: 确保所有数据基于463眼版本
2. **方法学透明**: 在图表说明中注明样本量和统计方法
3. **临床意义**: 突出关键发现和临床相关性
4. **美观与清晰**: 平衡信息密度和可读性
5. **期刊合规**: 检查目标期刊的图表要求

---
*计划生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*
*基于: 最终版/03_Tables/ 中的463眼版本数据*
*目标: 制作符合SCI发表要求的高质量图表*
"""

plan_path = os.path.join("/root/.openclaw/workspace", "figure_production_plan.md")
with open(plan_path, 'w', encoding='utf-8') as f:
    f.write(production_plan)

print(f"✅ 制作计划: {plan_path}")

print(f"\n" + "=" * 80)
print("🎯 图表制作分析完成!")
print(f"建议制作 {len(figures_needed)} 个SCI标准图表")
print("已生成制作计划和脚本框架")
print("=" * 80)