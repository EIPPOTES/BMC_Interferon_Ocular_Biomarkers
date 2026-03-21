#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
选项B：全面更新统计结果（基于485眼数据）
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import os

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 输出目录
OUTPUT_DIR = '/mnt/c/Users/CUI/Desktop/投稿、数据修改'

def load_data():
    """加载485眼数据"""
    df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx')
    
    # 排除年龄缺失的Control参与者
    control_ages = df[df['分组'] == '健康对照'].drop_duplicates('Patient_ID')[['Patient_ID', '年龄']]
    age_missing_ids = control_ages[control_ages['年龄'].isna()]['Patient_ID'].tolist()
    df_filtered = df[~((df['分组'] == '健康对照') & (df['Patient_ID'].isin(age_missing_ids)))]
    
    return df_filtered

def calculate_cohens_d(x, y):
    """计算Cohen's d"""
    nx, ny = len(x), len(y)
    pooled_std = np.sqrt(((nx-1)*x.var() + (ny-1)*y.var()) / (nx + ny - 2))
    return (x.mean() - y.mean()) / pooled_std

def generate_all_results():
    """生成所有统计结果"""
    
    df = load_data()
    
    print("="*70)
    print("选项B：全面统计更新（485眼数据）")
    print("="*70)
    
    # 1. 样本量
    print("\n1. 样本量")
    print("-"*70)
    
    mdd_patients = df[df['分组'] == '抑郁症']['Patient_ID'].nunique()
    control_patients = df[df['分组'] == '健康对照']['Patient_ID'].nunique()
    mdd_eyes = len(df[df['分组'] == '抑郁症'])
    control_eyes = len(df[df['分组'] == '健康对照'])
    
    print(f"  患者总数: {mdd_patients + control_patients} ({mdd_patients} MDD + {control_patients} Control)")
    print(f"  眼数总数: {mdd_eyes + control_eyes} ({mdd_eyes} MDD + {control_eyes} Control)")
    
    # 2. 主要指标统计
    print("\n2. 主要OCT指标统计")
    print("-"*70)
    
    indicators = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal Thickness'),
        ('Retina_内环颞侧', 'Inner Temporal Thickness'),
        ('RNFL_上方', 'Superior RNFL'),
    ]
    
    results = []
    for col, name in indicators:
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        # 统计检验
        stat, pvalue = stats.mannwhitneyu(mdd_data, control_data, alternative='two-sided')
        
        # 效应量
        d = calculate_cohens_d(mdd_data, control_data)
        
        results.append({
            'name': name,
            'mdd_n': len(mdd_data),
            'mdd_mean': mdd_data.mean(),
            'mdd_std': mdd_data.std(),
            'ctrl_n': len(control_data),
            'ctrl_mean': control_data.mean(),
            'ctrl_std': control_data.std(),
            'pvalue': pvalue,
            'cohens_d': d
        })
        
        print(f"\n  {name}:")
        print(f"    MDD:     n={len(mdd_data)}, {mdd_data.mean():.2f}±{mdd_data.std():.2f}")
        print(f"    Control: n={len(control_data)}, {control_data.mean():.2f}±{control_data.std():.2f}")
        print(f"    P值:     {pvalue:.6f}")
        print(f"    Cohen's d: {d:.3f}")
    
    # 保存结果
    results_df = pd.DataFrame(results)
    results_df.to_excel(f'{OUTPUT_DIR}/03_Tables/Table2_OCT_Parameters_Updated.xlsx', index=False)
    print(f"\n✅ 统计结果已保存")
    
    return results

if __name__ == "__main__":
    generate_all_results()
    print("\n" + "="*70)
    print("✅ 选项B统计更新完成!")
    print("="*70)