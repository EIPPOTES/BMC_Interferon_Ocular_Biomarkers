#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
修复Figure 2 - 确保n值与实际数据量一致
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def load_and_filter_data():
    """加载并过滤数据"""
    df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx')
    
    # 排除年龄缺失的Control参与者
    control_patients = df[df['分组'] == '健康对照']['Patient_ID'].unique()
    control_ages = df[df['分组'] == '健康对照'].drop_duplicates('Patient_ID')[['Patient_ID', '年龄']]
    age_missing_ids = control_ages[control_ages['年龄'].isna()]['Patient_ID'].tolist()
    
    # 过滤数据
    df_filtered = df[~((df['分组'] == '健康对照') & (df['Patient_ID'].isin(age_missing_ids)))]
    
    return df_filtered

def get_actual_n(df, column):
    """获取实际用于分析的n值"""
    mdd_data = df[(df['分组'] == '抑郁症') & (df[column].notna())][column]
    control_data = df[(df['分组'] == '健康对照') & (df[column].notna())][column]
    return len(mdd_data), len(control_data)

def generate_figure2_fixed():
    """生成修复后的Figure 2"""
    
    df = load_and_filter_data()
    
    # 定义指标和显示名称
    indicators = [
        ('Retina_平均厚度', 'Mean Macular\nThickness (μm)'),
        ('Retina_外环颞侧', 'Outer Temporal\nThickness (μm)'),
        ('Retina_内环颞侧', 'Inner Temporal\nThickness (μm)'),
        ('RNFL_上方', 'Superior RNFL\nThickness (μm)'),
    ]
    
    # 颜色设置
    colors = {'抑郁症': '#E69F00', '健康对照': '#56B4E9'}
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    for idx, (col, title) in enumerate(indicators):
        ax = axes[idx]
        
        # 获取数据
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        # 获取实际n值
        n_mdd = len(mdd_data)
        n_control = len(control_data)
        
        # 准备箱图数据
        data_to_plot = [mdd_data, control_data]
        
        # 绘制箱图
        bp = ax.boxplot(data_to_plot, labels=['MDD', 'Control'], patch_artist=True,
                       widths=0.6, showfliers=True)
        
        # 设置颜色
        bp['boxes'][0].set_facecolor(colors['抑郁症'])
        bp['boxes'][1].set_facecolor(colors['健康对照'])
        bp['boxes'][0].set_alpha(0.7)
        bp['boxes'][1].set_alpha(0.7)
        
        # 设置标题和标签
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_ylabel('Thickness (μm)', fontsize=10)
        
        # 添加n值标注 (在箱图上方)
        ax.text(1, ax.get_ylim()[1] * 0.98, f'n={n_mdd}', ha='center', fontsize=10, fontweight='bold')
        ax.text(2, ax.get_ylim()[1] * 0.98, f'n={n_control}', ha='center', fontsize=10, fontweight='bold')
        
        # 添加统计检验结果
        from scipy import stats
        stat, pvalue = stats.mannwhitneyu(mdd_data, control_data, alternative='two-sided')
        
        # 添加显著性标记
        y_max = max(mdd_data.max(), control_data.max())
        y_min = min(mdd_data.min(), control_data.min())
        y_range = y_max - y_min
        
        # 绘制显著性线
        ax.plot([1, 2], [y_max + y_range*0.05, y_max + y_range*0.05], 'k-', linewidth=1.5)
        
        if pvalue < 0.001:
            sig_text = '***'
        elif pvalue < 0.01:
            sig_text = '**'
        elif pvalue < 0.05:
            sig_text = '*'
        else:
            sig_text = 'ns'
        
        ax.text(1.5, y_max + y_range*0.06, sig_text, ha='center', fontsize=14, fontweight='bold')
        
        # 添加P值
        ax.text(1.5, y_max + y_range*0.12, f'P={pvalue:.4f}', ha='center', fontsize=9)
    
    # 总标题
    fig.suptitle('Figure 2. Comparison of OCT Parameters Between MDD and Control Groups', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    # 图注
    note_text = (
        "Note: Box plots show median (horizontal line), interquartile range (box), "
        "and whiskers (1.5×IQR). Outliers shown as individual points. "
        "n values indicate actual number of eyes with valid data for each parameter. "
        "Statistical significance: *P<0.05, **P<0.01, ***P<0.001."
    )
    fig.text(0.5, 0.02, note_text, ha='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    
    # 保存
    output_path = '/mnt/c/Users/CUI/Desktop/投稿版/Figures_修订版/Figure2_Group_Comparison_修正版.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"✅ Figure 2已保存: {output_path}")
    print("\n各指标实际n值:")
    for col, title in indicators:
        n_mdd, n_control = get_actual_n(df, col)
        print(f"  {title.split('(')[0].strip()}: MDD n={n_mdd}, Control n={n_control}")

if __name__ == "__main__":
    print("="*60)
    print("生成修复后的Figure 2")
    print("="*60)
    generate_figure2_fixed()
    print("\n" + "="*60)
    print("✅ Figure 2修复完成!")
    print("="*60)