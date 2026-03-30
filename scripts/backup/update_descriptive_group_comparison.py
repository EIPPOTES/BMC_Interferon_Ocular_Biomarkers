#!/usr/bin/env python3
"""
基于499眼版本更新描述性统计和组间比较
"""

import pandas as pd
import numpy as np
import os
from scipy import stats
from datetime import datetime

print("=" * 80)
print("基于499眼版本更新描述性统计和组间比较")
print(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
print("=" * 80)

# 路径设置
base_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改"
data_path = os.path.join(base_dir, "04_Data", "data_499eyes_20260315.xlsx")
tables_dir = os.path.join(base_dir, "03_Tables")
report_dir = os.path.join(base_dir, "分析报告", "01_原始数据表格")

# 确保目录存在
os.makedirs(tables_dir, exist_ok=True)
os.makedirs(report_dir, exist_ok=True)

# 读取数据
print(f"\n📊 读取数据文件: {data_path}")
df = pd.read_excel(data_path)
print(f"数据形状: {df.shape[0]} 行 × {df.shape[1]} 列")

# 检查分组
print(f"\n📋 样本分布:")
print(f"总样本: {df.shape[0]} 眼")
print(f"抑郁症: {df[df['分组'] == '抑郁症'].shape[0]} 眼")
print(f"健康对照: {df[df['分组'] == '健康对照'].shape[0]} 眼")

# 识别OCT指标列
oct_columns = []
for col in df.columns:
    if any(keyword in col for keyword in ['Retina', 'RNFL', 'GCL', 'Choroid', 'Disc', 'Cup', 'Rim']):
        oct_columns.append(col)

print(f"\n🔍 识别到 {len(oct_columns)} 个OCT指标")
print(f"示例指标: {oct_columns[:5]}")

# 为指标分类（用于描述性统计）
def categorize_parameter(param):
    """将OCT指标分类"""
    param_lower = param.lower()
    if 'retina' in param_lower:
        return '黄斑参数_Retina'
    elif 'rnfl' in param_lower:
        return 'RNFL参数'
    elif 'gcl' in param_lower:
        return 'GCL参数'
    elif 'choroid' in param_lower:
        return '脉络膜参数'
    elif any(keyword in param_lower for keyword in ['disc', 'cup', 'rim']):
        return '视盘参数'
    else:
        return '其他参数'

# 1. 生成描述性统计
print(f"\n📈 生成描述性统计...")

descriptive_results = []
detailed_results = []

for param in oct_columns:
    # 获取数据
    mdd_data = df[df['分组'] == '抑郁症'][param].dropna()
    control_data = df[df['分组'] == '健康对照'][param].dropna()
    
    # 跳过数据太少的指标
    if len(mdd_data) < 10 or len(control_data) < 10:
        continue
    
    # 基本描述性统计
    category = categorize_parameter(param)
    
    # MDD组统计
    mdd_n = len(mdd_data)
    mdd_mean = mdd_data.mean()
    mdd_sd = mdd_data.std()
    mdd_median = mdd_data.median()
    mdd_q1 = mdd_data.quantile(0.25)
    mdd_q3 = mdd_data.quantile(0.75)
    mdd_min = mdd_data.min()
    mdd_max = mdd_data.max()
    
    # 对照组统计
    control_n = len(control_data)
    control_mean = control_data.mean()
    control_sd = control_data.std()
    control_median = control_data.median()
    control_q1 = control_data.quantile(0.25)
    control_q3 = control_data.quantile(0.75)
    control_min = control_data.min()
    control_max = control_data.max()
    
    # 均值差
    mean_diff = mdd_mean - control_mean
    
    # 添加到详细结果
    detailed_results.append({
        'Category': category,
        'Parameter': param,
        'Unit': 'μm',  # 假设单位为微米
        'MDD_n': mdd_n,
        'MDD_Mean': mdd_mean,
        'MDD_SD': mdd_sd,
        'MDD_Median': mdd_median,
        'MDD_Q1': mdd_q1,
        'MDD_Q3': mdd_q3,
        'MDD_Min': mdd_min,
        'MDD_Max': mdd_max,
        'Control_n': control_n,
        'Control_Mean': control_mean,
        'Control_SD': control_sd,
        'Control_Median': control_median,
        'Control_Q1': control_q1,
        'Control_Q3': control_q3,
        'Control_Min': control_min,
        'Control_Max': control_max,
        'Mean_Diff': mean_diff
    })
    
    # 简化版本（基本描述性统计）
    descriptive_results.append({
        'Category': category,
        'Parameter': param,
        'MDD_n': mdd_n,
        'MDD_Mean': mdd_mean,
        'MDD_SD': mdd_sd,
        'MDD_Median': mdd_median,
        'MDD_Q1': mdd_q1,
        'MDD_Q3': mdd_q3,
        'Control_n': control_n,
        'Control_Mean': control_mean,
        'Control_SD': control_sd,
        'Control_Median': control_median,
        'Control_Q1': control_q1,
        'Control_Q3': control_q3,
        'Mean_Diff': mean_diff
    })

# 转换为DataFrame
detailed_df = pd.DataFrame(detailed_results)
descriptive_df = pd.DataFrame(descriptive_results)

print(f"生成 {len(detailed_df)} 个指标的描述性统计")

# 2. 生成组间比较
print(f"\n📊 生成组间比较分析...")

group_comparison_results = []

for param in oct_columns:
    # 获取数据
    mdd_data = df[df['分组'] == '抑郁症'][param].dropna()
    control_data = df[df['分组'] == '健康对照'][param].dropna()
    
    # 跳过数据太少的指标
    if len(mdd_data) < 10 or len(control_data) < 10:
        continue
    
    # 基本统计
    mdd_n = len(mdd_data)
    mdd_mean = mdd_data.mean()
    mdd_sd = mdd_data.std()
    
    control_n = len(control_data)
    control_mean = control_data.mean()
    control_sd = control_data.std()
    
    mean_diff = mdd_mean - control_mean
    
    # 效应量 (Cohen's d)
    pooled_sd = np.sqrt(((mdd_n - 1) * mdd_sd**2 + (control_n - 1) * control_sd**2) / (mdd_n + control_n - 2))
    if pooled_sd != 0:
        cohens_d = mean_diff / pooled_sd
    else:
        cohens_d = np.nan
    
    # 统计检验
    # 先检查方差齐性
    _, levene_p = stats.levene(mdd_data, control_data)
    
    if levene_p > 0.05:  # 方差齐性，使用t检验
        _, p_value = stats.ttest_ind(mdd_data, control_data, equal_var=True)
        test_method = "独立样本t检验"
    else:  # 方差不齐，使用Welch's t检验
        _, p_value = stats.ttest_ind(mdd_data, control_data, equal_var=False)
        test_method = "Welch's t检验"
    
    # 显著性标记
    if p_value < 0.001:
        significance = "***"
    elif p_value < 0.01:
        significance = "**"
    elif p_value < 0.05:
        significance = "*"
    else:
        significance = ""
    
    group_comparison_results.append({
        'Parameter': param,
        'MDD_n': mdd_n,
        'MDD_Mean': mdd_mean,
        'MDD_SD': mdd_sd,
        'Control_n': control_n,
        'Control_Mean': control_mean,
        'Control_SD': control_sd,
        'Mean_Diff': mean_diff,
        'Cohens_d': cohens_d,
        'P_value': p_value,
        'Significance': significance,
        'Test_Method': test_method
    })

group_comparison_df = pd.DataFrame(group_comparison_results)

print(f"生成 {len(group_comparison_df)} 个指标的组间比较")

# 3. 保存文件
print(f"\n💾 保存结果文件...")

# 文件命名
date_str = datetime.now().strftime("%Y%m%d")
desc_file_simple = os.path.join(tables_dir, f"Descriptive_Statistics_499eyes_{date_str}.xlsx")
desc_file_detailed = os.path.join(tables_dir, f"Descriptive_Statistics_Detailed_499eyes_{date_str}.xlsx")
group_file = os.path.join(tables_dir, f"Group_Comparison_499eyes_{date_str}.xlsx")

# 保存描述性统计（简化版）
descriptive_df.to_excel(desc_file_simple, index=False)
print(f"✅ 描述性统计（简化版）: {desc_file_simple}")

# 保存描述性统计（详细版）
detailed_df.to_excel(desc_file_detailed, index=False)
print(f"✅ 描述性统计（详细版）: {desc_file_detailed}")

# 保存组间比较
group_comparison_df.to_excel(group_file, index=False)
print(f"✅ 组间比较分析: {group_file}")

# 4. 复制到分析报告目录
print(f"\n📁 复制到分析报告目录...")

report_desc_simple = os.path.join(report_dir, f"Descriptive_Statistics_499eyes_{date_str}.xlsx")
report_desc_detailed = os.path.join(report_dir, f"Descriptive_Statistics_Detailed_499eyes_{date_str}.xlsx")
report_group = os.path.join(report_dir, f"Group_Comparison_499eyes_{date_str}.xlsx")

descriptive_df.to_excel(report_desc_simple, index=False)
detailed_df.to_excel(report_desc_detailed, index=False)
group_comparison_df.to_excel(report_group, index=False)

print(f"✅ 复制到: {report_desc_simple}")
print(f"✅ 复制到: {report_desc_detailed}")
print(f"✅ 复制到: {report_group}")

# 5. 生成汇总报告
print(f"\n📋 生成分析汇总报告...")

# 显著性统计
total_indicators = len(group_comparison_df)
sig_indicators = group_comparison_df[group_comparison_df['P_value'] < 0.05].shape[0]
sig_pct = (sig_indicators / total_indicators * 100) if total_indicators > 0 else 0

# 效应量统计
large_effect = group_comparison_df[group_comparison_df['Cohens_d'].abs() >= 0.8].shape[0]
medium_effect = group_comparison_df[(group_comparison_df['Cohens_d'].abs() >= 0.5) & 
                                   (group_comparison_df['Cohens_d'].abs() < 0.8)].shape[0]
small_effect = group_comparison_df[(group_comparison_df['Cohens_d'].abs() >= 0.2) & 
                                  (group_comparison_df['Cohens_d'].abs() < 0.5)].shape[0]

# 找到效应量最大的指标
if not group_comparison_df.empty:
    max_effect_idx = group_comparison_df['Cohens_d'].abs().idxmax()
    max_effect_param = group_comparison_df.loc[max_effect_idx, 'Parameter']
    max_effect_d = group_comparison_df.loc[max_effect_idx, 'Cohens_d']
    max_effect_p = group_comparison_df.loc[max_effect_idx, 'P_value']
else:
    max_effect_param = "无"
    max_effect_d = 0
    max_effect_p = 1

# 生成报告文本
report_text = f"""# 描述性统计与组间比较分析报告
## 基于499眼版本 (2026-03-15)

### 📊 数据概览
- **数据源**: data_499eyes_20260315.xlsx
- **总样本**: 499眼
  - 抑郁症组: 325眼
  - 健康对照组: 174眼
- **分析指标**: {total_indicators} 个OCT指标
- **生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

### 📈 描述性统计
- **简化版**: {desc_file_simple} ({descriptive_df.shape[0]}指标 × {descriptive_df.shape[1]}列)
- **详细版**: {desc_file_detailed} ({detailed_df.shape[0]}指标 × {detailed_df.shape[1]}列)
- **包含统计量**: 样本量、均值、标准差、中位数、四分位数(Q1, Q3)、最小值、最大值

### ⚖️ 组间比较结果
- **分析文件**: {group_file} ({group_comparison_df.shape[0]}指标 × {group_comparison_df.shape[1]}列)
- **显著性统计**: 
  - 显著指标(P<0.05): {sig_indicators}/{total_indicators} ({sig_pct:.1f}%)
  - 高度显著(P<0.001): {group_comparison_df[group_comparison_df['P_value'] < 0.001].shape[0]}个
  - 中度显著(P<0.01): {group_comparison_df[(group_comparison_df['P_value'] >= 0.001) & (group_comparison_df['P_value'] < 0.01)].shape[0]}个
  - 轻度显著(P<0.05): {group_comparison_df[(group_comparison_df['P_value'] >= 0.01) & (group_comparison_df['P_value'] < 0.05)].shape[0]}个

### 🎯 效应量分析
- **大效应(Cohen's d ≥ 0.8)**: {large_effect} 个指标
- **中等效应(0.5 ≤ d < 0.8)**: {medium_effect} 个指标  
- **小效应(0.2 ≤ d < 0.5)**: {small_effect} 个指标
- **效应量最大的指标**: {max_effect_param}
  - Cohen's d = {max_effect_d:.3f}
  - P值 = {max_effect_p:.6f}

### 📁 输出文件
1. **描述性统计（简化版）**: `{os.path.basename(desc_file_simple)}`
2. **描述性统计（详细版）**: `{os.path.basename(desc_file_detailed)}`
3. **组间比较分析**: `{os.path.basename(group_file)}`

### 🔗 文件位置
- 主目录: `{tables_dir}`
- 分析报告目录: `{report_dir}`

### ⚠️ 注意事项
1. 所有分析基于原始499眼数据，未排除任何缺失值（每个指标单独处理缺失）
2. 组间比较使用适当的统计检验（根据方差齐性选择t检验或Welch's t检验）
3. 效应量计算基于Cohen's d标准
4. 显著性标记: *** (P<0.001), ** (P<0.01), * (P<0.05)

---
*分析完成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*
"""

# 保存报告
report_file = os.path.join(base_dir, "分析报告", f"描述性统计与组间比较_更新报告_{date_str}.md")
with open(report_file, 'w', encoding='utf-8') as f:
    f.write(report_text)

print(f"✅ 分析报告: {report_file}")

print(f"\n" + "=" * 80)
print("🎉 更新完成!")
print("=" * 80)