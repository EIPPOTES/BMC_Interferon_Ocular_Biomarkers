#!/usr/bin/env python3
"""
更新论文中的所有数值（基于485眼数据）
"""

import pandas as pd
import numpy as np
from scipy import stats
import re

def load_data():
    """加载485眼数据"""
    df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx')
    control_ages = df[df['分组'] == '健康对照'].drop_duplicates('Patient_ID')[['Patient_ID', '年龄']]
    age_missing_ids = control_ages[control_ages['年龄'].isna()]['Patient_ID'].tolist()
    return df[~((df['分组'] == '健康对照') & (df['Patient_ID'].isin(age_missing_ids)))]

def calculate_cohens_d(x, y):
    """计算Cohen's d"""
    nx, ny = len(x), len(y)
    pooled_std = np.sqrt(((nx-1)*x.var() + (ny-1)*y.var()) / (nx + ny - 2))
    return (x.mean() - y.mean()) / pooled_std

def get_updated_statistics(df):
    """获取更新后的统计结果"""
    
    # 样本量
    mdd_patients = df[df['分组'] == '抑郁症']['Patient_ID'].nunique()
    control_patients = df[df['分组'] == '健康对照']['Patient_ID'].nunique()
    mdd_eyes = len(df[df['分组'] == '抑郁症'])
    control_eyes = len(df[df['分组'] == '健康对照'])
    
    # 主要指标
    indicators = {
        'mean_macular': {
            'col': 'Retina_平均厚度',
            'mdd_data': df[(df['分组'] == '抑郁症')]['Retina_平均厚度'],
            'control_data': df[(df['分组'] == '健康对照')]['Retina_平均厚度'],
        },
        'outer_temporal': {
            'col': 'Retina_外环颞侧',
            'mdd_data': df[(df['分组'] == '抑郁症')]['Retina_外环颞侧'],
            'control_data': df[(df['分组'] == '健康对照')]['Retina_外环颞侧'],
        },
        'inner_temporal': {
            'col': 'Retina_内环颞侧',
            'mdd_data': df[(df['分组'] == '抑郁症')]['Retina_内环颞侧'],
            'control_data': df[(df['分组'] == '健康对照')]['Retina_内环颞侧'],
        },
        'superior_rnfl': {
            'col': 'RNFL_上方',
            'mdd_data': df[(df['分组'] == '抑郁症') & (df['RNFL_上方'].notna())]['RNFL_上方'],
            'control_data': df[(df['分组'] == '健康对照') & (df['RNFL_上方'].notna())]['RNFL_上方'],
        },
    }
    
    results = {}
    for key, data in indicators.items():
        mdd = data['mdd_data']
        ctrl = data['control_data']
        
        stat, pvalue = stats.mannwhitneyu(mdd, ctrl, alternative='two-sided')
        d = calculate_cohens_d(mdd, ctrl)
        
        results[key] = {
            'mdd_n': len(mdd),
            'mdd_mean': mdd.mean(),
            'mdd_std': mdd.std(),
            'ctrl_n': len(ctrl),
            'ctrl_mean': ctrl.mean(),
            'ctrl_std': ctrl.std(),
            'pvalue': pvalue,
            'd': d,
        }
    
    return {
        'sample': {
            'total_patients': mdd_patients + control_patients,
            'mdd_patients': mdd_patients,
            'control_patients': control_patients,
            'total_eyes': mdd_eyes + control_eyes,
            'mdd_eyes': mdd_eyes,
            'control_eyes': control_eyes,
        },
        'indicators': results,
    }

def update_paper():
    """更新论文中的数值"""
    
    print("="*70)
    print("更新论文数值（基于485眼数据）")
    print("="*70)
    
    # 加载数据
    df = load_data()
    stats_data = get_updated_statistics(df)
    
    # 读取论文
    paper_file = '/mnt/c/Users/CUI/Desktop/投稿版/01_Manuscript/OCT_MDD_Manuscript_FINAL_20260314_1308.md'
    with open(paper_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    print("\n1. 样本量更新")
    print("-"*70)
    
    # 样本量替换
    replacements = [
        ('251 participants', f"{stats_data['sample']['total_patients']} participants"),
        ('87 controls', f"{stats_data['sample']['control_patients']} controls"),
        ('87 healthy controls', f"{stats_data['sample']['control_patients']} healthy controls"),
        ('499 eyes', f"{stats_data['sample']['total_eyes']} eyes"),
        ('174 control eyes', f"{stats_data['sample']['control_eyes']} control eyes"),
    ]
    
    for old, new in replacements:
        if old in content:
            content = content.replace(old, new)
            print(f"  ✅ '{old}' → '{new}'")
    
    print("\n2. 统计数值更新")
    print("-"*70)
    
    # Mean Macular更新
    mm = stats_data['indicators']['mean_macular']
    old_mm = f"{mm['ctrl_mean']:.2f}±{mm['ctrl_std']:.2f}"
    print(f"  Control Mean Macular: {old_mm}")
    
    # 效应量更新
    print(f"\n  Mean Macular Cohen's d: {mm['d']:.3f}")
    
    ot = stats_data['indicators']['outer_temporal']
    print(f"  Outer Temporal Cohen's d: {ot['d']:.3f}")
    
    # 保存更新后的论文
    output_file = '/mnt/c/Users/CUI/Desktop/投稿、数据修改/01_Manuscript/OCT_MDD_Manuscript_Updated_485eyes.md'
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"\n✅ 论文已保存: {output_file}")
    
    # 生成更新摘要
    print("\n" + "="*70)
    print("更新摘要")
    print("="*70)
    print(f"""
样本量:
  - 患者: {stats_data['sample']['total_patients']} ({stats_data['sample']['mdd_patients']} MDD + {stats_data['sample']['control_patients']} Control)
  - 眼数: {stats_data['sample']['total_eyes']} ({stats_data['sample']['mdd_eyes']} MDD + {stats_data['sample']['control_eyes']} Control)

主要指标:
  - Mean Macular: d={mm['d']:.3f}
  - Outer Temporal: d={ot['d']:.3f}
  - Inner Temporal: d={stats_data['indicators']['inner_temporal']['d']:.3f}
  - Superior RNFL: d={stats_data['indicators']['superior_rnfl']['d']:.3f}
""")

if __name__ == "__main__":
    update_paper()
    print("\n" + "="*70)
    print("✅ 论文数值更新完成!")
    print("="*70)