#!/usr/bin/env python3
"""
修复Figure 5 - 添加杯盘比(C/D Ratio)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11

# 色盲友好配色
cb_palette = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999']

output_dir = '/mnt/c/Users/CUI/Desktop/投稿版/Figures_修订版'

def load_data():
    data_file = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx'
    df = pd.read_excel(data_file)
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    
    # 计算杯盘比
    df_patient['C/D_Area_Ratio'] = df_patient['Cup Area'] / df_patient['Disc Area']
    
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    control = df_patient[df_patient['分组'] == '健康对照']
    return df_patient, mdd, control

def fix_figure5_with_cd_ratio(df_patient, mdd, control):
    """Figure 5: 添加杯盘比"""
    print("修复 Figure 5: 添加杯盘比(C/D Ratio)...")
    
    def cohens_d(x, y):
        nx, ny = len(x), len(y)
        dof = nx + ny - 2
        return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1)**2 + (ny-1)*np.std(y, ddof=1)**2) / dof)
    
    def cohens_d_ci(d, nx, ny):
        se = np.sqrt((nx + ny) / (nx * ny) + d**2 / (2 * (nx + ny)))
        ci_lower = d - 1.96 * se
        ci_upper = d + 1.96 * se
        return ci_lower, ci_upper
    
    # 包含杯盘比的参数列表
    effect_params = [
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
        ('RNFL_Total', 'Total RNFL'),
        ('Cup Area', 'Cup Area'),
        ('Rim Volume', 'Rim Volume'),
        ('C/D_Area_Ratio', 'C/D Area Ratio'),  # 添加杯盘比
    ]
    
    effects = []
    ci_lowers = []
    ci_uppers = []
    labels = []
    p_values = []
    
    print("  效应量计算:")
    for col, name in effect_params:
        if col in df_patient.columns:
            mdd_data = mdd[col].dropna()
            ctrl_data = control[col].dropna()
            
            if len(mdd_data) > 0 and len(ctrl_data) > 0:
                d = cohens_d(mdd_data, ctrl_data)
                ci_lower, ci_upper = cohens_d_ci(d, len(mdd_data), len(ctrl_data))
                _, p = stats.mannwhitneyu(mdd_data, ctrl_data, alternative='two-sided')
                
                effects.append(d)
                ci_lowers.append(ci_lower)
                ci_uppers.append(ci_upper)
                labels.append(name)
                p_values.append(p)
                
                sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
                print(f"    {name}: d={d:.2f} {sig}, 95%CI=({ci_lower:.2f}, {ci_upper:.2f})")
    
    # 绘制森林图
    fig, ax = plt.subplots(figsize=(13, 10))
    
    y_pos = np.arange(len(labels))
    
    # 绘制每个效应量
    for i, (effect, ci_lower, ci_upper, p, label) in enumerate(
        zip(effects, ci_lowers, ci_uppers, p_values, labels)):
        
        # 根据效应量方向选择颜色
        # 负值(减少)用红色，正值(增加)用蓝色
        if effect < 0:
            color = cb_palette[0]  # 红色 - Reduced in MDD
        else:
            color = cb_palette[1]  # 蓝色 - Increased in MDD
        
        # 杯盘比用粗线突出
        if 'C/D' in label:
            linewidth = 3
            markersize = 120
        else:
            linewidth = 2
            markersize = 100
        
        # CI线
        ax.plot([ci_lower, ci_upper], [i, i], color=color, linewidth=linewidth, alpha=0.7)
        
        # 效应量点
        sig_marker = '*' if p < 0.05 else ''
        ax.scatter(effect, i, color=color, s=markersize, zorder=3, 
                  edgecolors='black', linewidth=1.5)
        
        # 数值标签
        ax.text(effect, i+0.2, f'{effect:.2f}{sig_marker}', 
               ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # 参考线
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1.5)
    ax.axvline(x=-0.2, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(x=0.2, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    # 区域标注
    ax.text(-0.6, len(labels)+0.8, 'MDD < Control\n(Reduced in MDD)', 
           ha='center', fontsize=10, color=cb_palette[0], fontweight='bold')
    ax.text(0.6, len(labels)+0.8, 'MDD > Control\n(Increased in MDD)', 
           ha='center', fontsize=10, color=cb_palette[1], fontweight='bold')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel("Cohen's d (Effect Size)", fontsize=13, fontweight='bold')
    ax.set_xlim(-1, 1)
    ax.set_title('Figure 5. Effect Sizes of Retinal Changes in MDD Patients\n'
                'Including Cup-to-Disc Ratio (C/D Area Ratio). Horizontal lines represent 95% CI.',
                fontsize=13, fontweight='bold', pad=20)
    ax.grid(True, axis='x', alpha=0.3)
    
    # 图例
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=cb_palette[0], 
               markersize=10, label='Reduced in MDD (negative d)', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=cb_palette[1], 
               markersize=10, label='Increased in MDD (positive d)', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=cb_palette[6], 
               markersize=10, label='Optic Disc Parameter (C/D Ratio)', markeredgecolor='black'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
    
    # 添加杯盘比说明
    cd_text = (
        "Note: C/D Area Ratio = Cup Area / Disc Area. "
        "Higher values indicate larger cup relative to disc. "
        "MDD patients show increased C/D ratio, suggesting optic nerve changes."
    )
    fig.text(0.5, 0.02, cd_text, ha='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))
    
    # 样本量说明
    fig.text(0.02, 0.98, f'MDD n={len(mdd)}, Control n={len(control)}', 
            fontsize=9, verticalalignment='top')
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(f'{output_dir}/Figure5_Forest_Plot_含杯盘比版.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✅ 已保存: Figure5_Forest_Plot_含杯盘比版.png")
    print(f"     共{len(labels)}个参数，包括C/D Area Ratio")

def main():
    print("="*70)
    print("修复 Figure 5 - 添加杯盘比(C/D Ratio)")
    print("="*70)
    
    df_patient, mdd, control = load_data()
    
    # 检查杯盘比数据
    print(f"\n杯盘比统计:")
    print(f"  MDD: {mdd['C/D_Area_Ratio'].mean():.3f} ± {mdd['C/D_Area_Ratio'].std():.3f}")
    print(f"  Control: {control['C/D_Area_Ratio'].mean():.3f} ± {control['C/D_Area_Ratio'].std():.3f}")
    
    fix_figure5_with_cd_ratio(df_patient, mdd, control)
    
    print("\n" + "="*70)
    print("✅ Figure 5 修复完成!")
    print("="*70)

if __name__ == "__main__":
    main()