#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
计算OCT参数的Cohen's d效应量
"""

import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def cohens_d(group1, group2):
    """计算Cohen's d效应量"""
    n1, n2 = len(group1), len(group2)
    mean1, mean2 = np.mean(group1), np.mean(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    
    # 合并标准差
    pooled_sd = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
    if pooled_sd == 0:
        return 0
    return (mean1 - mean2) / pooled_sd

def calculate_all_effect_sizes(df):
    """计算所有选定参数的效应量"""
    
    # 分组 - 根据实际数据标签
    mdd_mask = df['分组'] == '抑郁症'
    control_mask = df['分组'] == '健康对照'
    
    if not control_mask.any():
        print("警告: 未找到对照组，检查分组列的值")
        print("分组列唯一值:", df['分组'].unique())
        return []
    
    print(f"MDD组标签: '抑郁症'")
    print(f"对照组标签: '健康对照'")
    
    mdd_group = df[mdd_mask]
    control_group = df[control_mask]
    
    print(f"MDD组样本量: {len(mdd_group)}")
    print(f"对照组样本量: {len(control_group)}")
    
    # 参数映射 (英文标签: 数据列名)
    param_mapping = {
        'Mean Macular Thickness': 'Retina_平均厚度',
        'Total Macular Volume': 'Retina_总体积', 
        'Outer Temporal': 'Retina_外环颞侧',
        'Inner Temporal': 'Retina_内环颞侧',
        'Superior RNFL': 'RNFL_上方',  # 或者 'RNFL_外环上方'
        'Cup Area': 'Cup Area',
        'Rim Volume': 'Rim Volume',
        'C/D Area Ratio': 'C/D Area Ratio'
    }
    
    results = []
    
    for eng_label, col_name in param_mapping.items():
        if col_name not in df.columns:
            print(f"警告: 列 '{col_name}' 不存在，跳过 {eng_label}")
            continue
        
        # 获取数据，移除NaN
        mdd_data = mdd_group[col_name].dropna()
        control_data = control_group[col_name].dropna()
        
        if len(mdd_data) < 5 or len(control_data) < 5:
            print(f"警告: {eng_label} 数据不足 (MDD: {len(mdd_data)}, Control: {len(control_data)})")
            continue
        
        # 计算效应量
        d = cohens_d(mdd_data, control_data)
        
        # 计算Mann-Whitney U检验p值（非参数检验，适合非正态分布）
        try:
            u_stat, p_value = stats.mannwhitneyu(mdd_data, control_data, alternative='two-sided')
        except:
            # 如果Mann-Whitney失败，使用t检验
            t_stat, p_value = stats.ttest_ind(mdd_data, control_data, equal_var=False, nan_policy='omit')
        
        # 计算均值
        mdd_mean = np.mean(mdd_data)
        control_mean = np.mean(control_data)
        
        results.append({
            'label': eng_label,
            'column': col_name,
            'cohens_d': d,
            'p_value': p_value,
            'mdd_mean': mdd_mean,
            'control_mean': control_mean,
            'mdd_n': len(mdd_data),
            'control_n': len(control_data)
        })
        
        print(f"{eng_label}: d = {d:.3f}, p = {p_value:.6f}, "
              f"MDD = {mdd_mean:.2f}, Control = {control_mean:.2f}")
    
    return results

def compare_with_hardcoded(results):
    """与硬编码值比较"""
    hardcoded = {
        'Mean Macular Thickness': (-0.415, 0.000003),
        'Total Macular Volume': (-0.416, 0.000003),
        'Outer Temporal': (-0.497, 0.000003),
        'Inner Temporal': (-0.375, 0.000032),
        'Superior RNFL': (-0.311, 0.002229),
        'Cup Area': (0.224, 0.022329),
        'Rim Volume': (-0.303, 0.010735),
        'C/D Area Ratio': (0.246, 0.021236)
    }
    
    print("\n与硬编码值比较:")
    print("="*60)
    print(f"{'参数':<25} {'计算d':<10} {'硬编码d':<10} {'差异':<10} {'计算p':<12} {'硬编码p':<12}")
    print("-"*60)
    
    for result in results:
        label = result['label']
        if label in hardcoded:
            h_d, h_p = hardcoded[label]
            calc_d = result['cohens_d']
            calc_p = result['p_value']
            diff_d = abs(calc_d - h_d)
            
            print(f"{label:<25} {calc_d:>9.3f} {h_d:>9.3f} {diff_d:>9.3f} "
                  f"{calc_p:>11.6f} {h_p:>11.6f}")

def main():
    data_file = '../data/raw/data.xlsx'
    
    try:
        print(f"读取数据文件: {data_file}")
        df = pd.read_excel(data_file)
        print(f"数据形状: {df.shape}")
        
        # 显示分组信息
        print(f"\n分组分布:")
        print(df['分组'].value_counts())
        
        # 计算效应量
        print(f"\n计算效应量:")
        results = calculate_all_effect_sizes(df)
        
        if results:
            compare_with_hardcoded(results)
            
            # 输出用于Figure 5的数据
            print(f"\n用于Figure 5的数据:")
            print("forest_data = [")
            for result in results:
                print(f"    ('{result['label']}', {result['cohens_d']:.3f}, {result['p_value']:.6f}),")
            print("]")
        else:
            print("未计算到有效结果")
            
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()