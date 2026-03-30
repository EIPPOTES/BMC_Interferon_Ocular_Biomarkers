#!/usr/bin/env python3
"""
生成选项B的剩余图表（Figure 3-6）
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

def generate_figure3():
    """生成Figure 3 - ROC曲线"""
    
    df = load_data()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 准备数据
    indicators = [
        ('Retina_平均厚度', 'Mean Macular', '#E69F00'),
        ('Retina_外环颞侧', 'Outer Temporal', '#56B4E9'),
    ]
    
    for col, name, color in indicators:
        # 获取数据
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        # 创建标签和数值
        y_true = [1] * len(mdd_data) + [0] * len(control_data)
        y_scores = list(mdd_data) + list(control_data)
        
        # 计算ROC
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        # 绘制ROC曲线
        ax.plot(fpr, tpr, color=color, lw=2, 
               label=f'{name} (AUC={roc_auc:.3f})')
    
    # 对角线
    ax.plot([0, 1], [0, 1], 'k--', lw=1, label='Chance level (AUC=0.500)')
    
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title('Figure 3. ROC Curves for OCT Parameters', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure3_ROC_Curves_485eyes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Figure 3已生成 (ROC曲线)")

def generate_figure4():
    """生成Figure 4 - 相关性分析"""
    
    df = load_data()
    
    # MDD组数据
    mdd_df = df[df['分组'] == '抑郁症'].copy()
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    indicators = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal Thickness'),
        ('RNFL_上方', 'Superior RNFL'),
    ]
    
    for idx, (col, title) in enumerate(indicators):
        ax = axes[idx]
        
        # 获取有效数据
        valid_data = mdd_df[(mdd_df[col].notna()) & (mdd_df['PHQ-9'].notna())]
        x = valid_data['PHQ-9']
        y = valid_data[col]
        
        # 散点图
        ax.scatter(x, y, alpha=0.5, s=30)
        
        # 回归线
        if len(x) > 1:
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            ax.plot(x, p(x), "r--", alpha=0.8, linewidth=2)
            
            # 计算相关系数
            r, pvalue = stats.pearsonr(x, y)
            ax.text(0.05, 0.95, f'r={r:.3f}, P={pvalue:.3f}', 
                   transform=ax.transAxes, fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlabel('PHQ-9 Score', fontsize=10)
        ax.set_ylabel(f'{title} (μm)', fontsize=10)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    fig.suptitle('Figure 4. Correlation Between OCT Parameters and Depression Severity', 
                fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure4_Correlation_485eyes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Figure 4已生成 (相关性分析)")

def generate_figure5():
    """生成Figure 5 - 森林图"""
    
    df = load_data()
    
    # 计算各指标的效应量
    indicators = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal Thickness'),
        ('Retina_内环颞侧', 'Inner Temporal Thickness'),
        ('RNFL_上方', 'Superior RNFL'),
    ]
    
    effect_sizes = []
    for col, name in indicators:
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        # Cohen's d
        nx, ny = len(mdd_data), len(control_data)
        pooled_std = np.sqrt(((nx-1)*mdd_data.var() + (ny-1)*control_data.var()) / (nx + ny - 2))
        d = (mdd_data.mean() - control_data.mean()) / pooled_std
        
        # 95% CI
        se = np.sqrt((nx + ny) / (nx * ny) + d**2 / (2 * (nx + ny)))
        ci_lower = d - 1.96 * se
        ci_upper = d + 1.96 * se
        
        effect_sizes.append({
            'name': name,
            'd': d,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper
        })
    
    # 绘制森林图
    fig, ax = plt.subplots(figsize=(12, 8))
    
    y_pos = np.arange(len(effect_sizes))
    
    for i, es in enumerate(effect_sizes):
        color = '#D55E00' if es['d'] < 0 else '#0072B2'
        ax.plot([es['ci_lower'], es['ci_upper']], [i, i], 'k-', linewidth=2)
        ax.plot(es['d'], i, 'o', markersize=10, color=color)
        ax.text(es['d'], i+0.2, f"d={es['d']:.2f}", ha='center', fontsize=9)
    
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([es['name'] for es in effect_sizes])
    ax.set_xlabel("Cohen's d (Effect Size)", fontsize=12)
    ax.set_title('Figure 5. Effect Sizes for OCT Parameters', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # 图例
    ax.text(0.02, 0.98, 'Red: MDD < Control\nBlue: MDD > Control', 
           transform=ax.transAxes, fontsize=9, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure5_Forest_Plot_485eyes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Figure 5已生成 (森林图)")

def generate_figure6():
    """生成Figure 6 - 亚组分析"""
    
    df = load_data()
    
    # 按年龄分组
    mdd_df = df[df['分组'] == '抑郁症'].copy()
    mdd_df['年龄组'] = pd.cut(mdd_df['年龄'], bins=[0, 30, 50, 100], labels=['<30', '30-50', '>50'])
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    age_groups = ['<30', '30-50', '>50']
    x_pos = np.arange(len(age_groups))
    width = 0.35
    
    means_mdd = []
    stds_mdd = []
    means_control = []
    
    control_data = df[df['分组'] == '健康对照']['Retina_平均厚度']
    control_mean = control_data.mean()
    
    for age_group in age_groups:
        group_data = mdd_df[mdd_df['年龄组'] == age_group]['Retina_平均厚度']
        means_mdd.append(group_data.mean())
        stds_mdd.append(group_data.std())
        means_control.append(control_mean)
    
    # 绘制条形图
    bars1 = ax.bar(x_pos - width/2, means_mdd, width, label='MDD', 
                  color='#E69F00', alpha=0.7, yerr=stds_mdd, capsize=5)
    bars2 = ax.bar(x_pos + width/2, means_control, width, label='Control', 
                  color='#56B4E9', alpha=0.7)
    
    ax.set_xlabel('Age Group', fontsize=12)
    ax.set_ylabel('Mean Macular Thickness (μm)', fontsize=12)
    ax.set_title('Figure 6. Subgroup Analysis by Age', fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(age_groups)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure6_Subgroup_Analysis_485eyes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Figure 6已生成 (亚组分析)")

if __name__ == "__main__":
    print("="*70)
    print("生成选项B剩余图表（Figure 3-6）")
    print("="*70)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    generate_figure3()
    generate_figure4()
    generate_figure5()
    generate_figure6()
    
    print("\n" + "="*70)
    print("✅ 所有图表已生成!")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*70)