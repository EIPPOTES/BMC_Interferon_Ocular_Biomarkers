#!/usr/bin/env python3
"""
生成增强版Figure 3和Figure 5（更多指标）
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from scipy import stats
import os

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

OUTPUT_DIR = '/mnt/c/Users/CUI/Desktop/投稿、数据修改/02_Figures'

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

def generate_figure3_enhanced():
    """生成增强版Figure 3 - 更多ROC曲线"""
    
    df = load_data()
    
    # 选择更多指标（8个）
    indicators = [
        ('Retina_平均厚度', 'Mean Macular', '#E69F00'),
        ('Retina_外环颞侧', 'Outer Temporal', '#56B4E9'),
        ('Retina_内环颞侧', 'Inner Temporal', '#009E73'),
        ('Retina_外环上方', 'Outer Superior', '#F0E442'),
        ('RNFL_上方', 'Superior RNFL', '#0072B2'),
        ('RNFL_颞侧', 'Temporal RNFL', '#D55E00'),
        ('Cup Area', 'Cup Area', '#CC79A7'),
        ('C/D Area Ratio', 'C/D Area Ratio', '#999999'),
    ]
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    results = []
    for col, name, color in indicators:
        # 获取数据
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        if len(mdd_data) < 10 or len(control_data) < 10:
            continue
        
        # 创建标签和数值
        y_true = [1] * len(mdd_data) + [0] * len(control_data)
        y_scores = list(mdd_data) + list(control_data)
        
        # 计算ROC
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        results.append({'name': name, 'auc': roc_auc, 'color': color})
        
        # 绘制ROC曲线
        ax.plot(fpr, tpr, color=color, lw=2.5, 
               label=f'{name} (AUC={roc_auc:.3f})')
    
    # 对角线
    ax.plot([0, 1], [0, 1], 'k--', lw=1.5, label='Chance level (AUC=0.500)')
    
    ax.set_xlabel('False Positive Rate', fontsize=13)
    ax.set_ylabel('True Positive Rate', fontsize=13)
    ax.set_title('Figure 3. ROC Curves for OCT Parameters (n=8 indicators)', 
                fontsize=15, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure3_ROC_Curves_Enhanced.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Figure 3增强版已生成 (8条ROC曲线)")
    print("\nAUC排名:")
    for r in sorted(results, key=lambda x: x['auc'], reverse=True):
        print(f"  {r['name']:20s}: AUC={r['auc']:.3f}")

def generate_figure5_enhanced():
    """生成增强版Figure 5 - 更多效应量"""
    
    df = load_data()
    
    # 选择更多指标（12个）
    indicators = [
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_中心厚度', 'Central Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_内环鼻侧', 'Inner Nasal'),
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
        ('RNFL_上方', 'Superior RNFL'),
        ('RNFL_颞侧', 'Temporal RNFL'),
        ('RNFL_鼻侧', 'Nasal RNFL'),
        ('Cup Area', 'Cup Area'),
        ('Rim Volume', 'Rim Volume'),
        ('C/D Area Ratio', 'C/D Area Ratio'),
    ]
    
    effect_sizes = []
    for col, name in indicators:
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        if len(mdd_data) < 10 or len(control_data) < 10:
            continue
        
        # Cohen's d
        d = calculate_cohens_d(mdd_data, control_data)
        
        # 95% CI
        nx, ny = len(mdd_data), len(control_data)
        se = np.sqrt((nx + ny) / (nx * ny) + d**2 / (2 * (nx + ny)))
        ci_lower = d - 1.96 * se
        ci_upper = d + 1.96 * se
        
        effect_sizes.append({
            'name': name,
            'd': d,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'abs_d': abs(d)
        })
    
    # 按效应量绝对值排序
    effect_sizes.sort(key=lambda x: x['abs_d'], reverse=True)
    
    # 绘制森林图
    fig, ax = plt.subplots(figsize=(14, 10))
    
    y_pos = np.arange(len(effect_sizes))
    
    for i, es in enumerate(effect_sizes):
        color = '#D55E00' if es['d'] < 0 else '#0072B2'
        # 置信区间线
        ax.plot([es['ci_lower'], es['ci_upper']], [i, i], 'k-', linewidth=2)
        # 效应量点
        ax.plot(es['d'], i, 'o', markersize=12, color=color, markeredgecolor='black', markeredgewidth=1.5)
        # 数值标签
        ax.text(es['d'], i+0.25, f"d={es['d']:.2f}", ha='center', fontsize=9, fontweight='bold')
    
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([es['name'] for es in effect_sizes], fontsize=11)
    ax.set_xlabel("Cohen's d (Effect Size)", fontsize=13)
    ax.set_title('Figure 5. Effect Sizes for OCT Parameters (n=12 indicators)', 
                fontsize=15, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x', linestyle='--')
    
    # 图例
    ax.text(0.02, 0.98, 'Red: MDD < Control (thinner)\nBlue: MDD > Control (thicker)', 
           transform=ax.transAxes, fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5, pad=0.5))
    
    # 效应量大小区间标注
    ax.axvspan(-0.2, 0.2, alpha=0.1, color='gray', label='Negligible')
    ax.axvspan(0.2, 0.5, alpha=0.1, color='yellow')
    ax.axvspan(-0.5, -0.2, alpha=0.1, color='yellow')
    ax.axvspan(0.5, 0.8, alpha=0.1, color='orange')
    ax.axvspan(-0.8, -0.5, alpha=0.1, color='orange')
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure5_Forest_Plot_Enhanced.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("\n✅ Figure 5增强版已生成 (12个效应量)")
    print("\n效应量排名:")
    for es in effect_sizes[:5]:
        direction = "MDD thinner" if es['d'] < 0 else "MDD thicker"
        print(f"  {es['name']:20s}: d={es['d']:+.3f} ({direction})")

if __name__ == "__main__":
    print("="*70)
    print("生成增强版Figure 3和Figure 5")
    print("="*70)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    generate_figure3_enhanced()
    print()
    generate_figure5_enhanced()
    
    print("\n" + "="*70)
    print("✅ 增强版图表已生成!")
    print("="*70)