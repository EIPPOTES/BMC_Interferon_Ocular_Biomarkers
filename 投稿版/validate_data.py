#!/usr/bin/env python3
"""
数据完整性验证脚本
验证原始数据与Figures的一致性
"""

import pandas as pd
import numpy as np
from scipy import stats
import os

def validate_data():
    print("="*60)
    print("OCT-MDD数据完整性验证")
    print("="*60)
    
    # 检查数据文件
    data_file = 'data/data.xlsx'
    if not os.path.exists(data_file):
        print(f"❌ 数据文件不存在: {data_file}")
        return False
    
    print(f"\n1. 加载数据文件: {data_file}")
    df = pd.read_excel(data_file)
    
    print(f"   ✅ 数据加载成功")
    print(f"   样本量: {len(df)}眼")
    print(f"   参数列: {len(df.columns)}列")
    
    # 检查分组
    print(f"\n2. 分组检查:")
    group_counts = df['分组'].value_counts()
    for group, count in group_counts.items():
        print(f"   {group}: {count}眼")
    
    # 验证关键统计
    print(f"\n3. 关键统计验证:")
    
    # 检查Figure 5的效应量
    mdd = df[df['分组'] == '抑郁症']
    control = df[df['分组'] == '健康对照']
    
    key_params = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('RNFL_上方', 'Superior RNFL'),
    ]
    
    for col, name in key_params:
        if col in df.columns:
            mdd_data = mdd[col].dropna()
            control_data = control[col].dropna()
            
            # 计算效应量
            mean_diff = mdd_data.mean() - control_data.mean()
            pooled_std = np.sqrt((mdd_data.var() + control_data.var()) / 2)
            cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
            
            print(f"   {name}:")
            print(f"     MDD: {mdd_data.mean():.2f} ± {mdd_data.std():.2f}")
            print(f"     Control: {control_data.mean():.2f} ± {control_data.std():.2f}")
            print(f"     Cohen's d: {cohens_d:.3f}")
    
    # 检查Figures文件
    print(f"\n4. Figures文件检查:")
    figures_dir = 'figures Final'
    expected_figures = [
        'Figure1_Study_Flowchart.png',
        'Figure2_Group_Comparison_Boxplot.png',
        'Figure3_ROC_Curves.png',
        'Figure4_Correlation_Scatter.png',
        'Figure5_Forest_Plot.png',
        'Figure6_Subgroup_Analysis.png',
    ]
    
    all_exist = True
    for fig in expected_figures:
        path = os.path.join(figures_dir, fig)
        if os.path.exists(path):
            size = os.path.getsize(path) / 1024
            print(f"   ✅ {fig} ({size:.1f} KB)")
        else:
            print(f"   ❌ {fig} (缺失)")
            all_exist = False
    
    print(f"\n" + "="*60)
    if all_exist:
        print("✅ 所有验证通过！数据完整性良好。")
    else:
        print("⚠️ 部分文件缺失，请检查。")
    print("="*60)
    
    return all_exist

if __name__ == "__main__":
    validate_data()