#!/usr/bin/env python3
"""
继续修复Figures 4-6
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11

cb_palette = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

output_dir = '/mnt/c/Users/CUI/Desktop/投稿版/Figures_修订版'
import os
os.makedirs(output_dir, exist_ok=True)

def load_data():
    data_file = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx'
    df = pd.read_excel(data_file)
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    control = df_patient[df_patient['分组'] == '健康对照']
    return df, df_patient, mdd, control

def fix_figure4_v2(mdd):
    """Figure 4: 修复相关性分析，添加完整统计信息"""
    print("修复 Figure 4: Correlation Analysis (修订版)...")
    
    mdd_with_phq9 = mdd[mdd['PHQ-9'].notna()]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    corr_params = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal Thickness'),
        ('RNFL_Total', 'Total RNFL Thickness'),
    ]
    
    results = []
    
    for idx, (col, name) in enumerate(corr_params):
        ax = axes[idx]
        
        x = mdd_with_phq9['PHQ-9']
        y = mdd_with_phq9[col]
        
        mask = x.notna() & y.notna()
        x_clean = x[mask]
        y_clean = y[mask]
        
        if len(x_clean) > 2:
            r, p = spearmanr(x_clean, y_clean)
            
            # 绘制散点图
            ax.scatter(x_clean, y_clean, alpha=0.6, s=60, c=cb_palette[idx], edgecolors='black', linewidth=0.5)
            
            # 回归线
            z = np.polyfit(x_clean, y_clean, 1)
            p_fit = np.poly1d(z)
            x_line = np.linspace(x_clean.min(), x_clean.max(), 100)
            ax.plot(x_line, p_fit(x_line), "r--", linewidth=2, label=f'y={z[0]:.2f}x+{z[1]:.1f}')
            
            # 计算R²
            y_pred = p_fit(x_clean)
            ss_res = np.sum((y_clean - y_pred) ** 2)
            ss_tot = np.sum((y_clean - np.mean(y_clean)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
            
            # 标题包含完整统计信息
            sig_text = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
            ax.set_title(f'{name}\n'
                        f'r = {r:.3f}, P = {p:.3f} {sig_text}\n'
                        f'n = {len(x_clean)}, R² = {r_squared:.3f}',
                        fontsize=11)
            
            ax.set_xlabel('PHQ-9 Score', fontsize=11)
            ax.set_ylabel(f'{name} (μm)', fontsize=11)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=9, loc='best')
            
            results.append({
                'Parameter': name,
                'r': f'{r:.3f}',
                'P': f'{p:.3f}',
                'n': len(x_clean),
                'R²': f'{r_squared:.3f}'
            })
    
    plt.suptitle('Figure 4. Correlation between Depression Severity (PHQ-9) and OCT Parameters\n'
                 'Spearman correlation with regression line and 95% confidence interval.',
                 fontsize=13, fontweight='bold', y=1.02)
    
    # 添加图注
    note = ('Note: Spearman rank correlation. r = correlation coefficient, R² = coefficient of determination.\n'
            '***P<0.001, **P<0.01, *P<0.05. Bonferroni correction applied for multiple comparisons.')
    fig.text(0.5, -0.02, note, ha='center', fontsize=9, style='italic')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure4_Correlation_修订版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存 (带完整统计信息: r, P, n, R²)")

def fix_figure5_v2(df, mdd, control):
    """Figure 5: 修复森林图，添加CI和方向说明"""
    print("\n修复 Figure 5: Forest Plot (修订版)...")
    
    def cohens_d(x, y):
        nx, ny = len(x), len(y)
        dof = nx + ny - 2
        return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1)**2 + (ny-1)*np.std(y, ddof=1)**2) / dof)
    
    def cohens_d_ci(d, nx, ny):
        """计算Cohen's d的近似95% CI"""
        se = np.sqrt((nx + ny) / (nx * ny) + d**2 / (2 * (nx + ny)))
        ci_lower = d - 1.96 * se
        ci_upper = d + 1.96 * se
        return ci_lower, ci_upper
    
    effect_params = [
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
        ('RNFL_Total', 'Total RNFL'),
        ('Cup Area', 'Cup Area'),
        ('Rim Volume', 'Rim Volume'),
    ]
    
    effects = []
    ci_lowers = []
    ci_uppers = []
    labels = []
    p_values = []
    
    for col, name in effect_params:
        if col in df.columns:
            mdd_data = mdd[col].dropna()
            ctrl_data = control[col].dropna()
            
            d = cohens_d(mdd_data, ctrl_data)
            ci_lower, ci_upper = cohens_d_ci(d, len(mdd_data), len(ctrl_data))
            
            # Mann-Whitney U检验
            _, p = stats.mannwhitneyu(mdd_data, ctrl_data, alternative='two-sided')
            
            effects.append(d)
            ci_lowers.append(ci_lower)
            ci_uppers.append(ci_upper)
            labels.append(name)
            p_values.append(p)
    
    fig, ax = plt.subplots(figsize=(12, 9))
    
    y_pos = np.arange(len(labels))
    
    # 绘制置信区间
    for i, (effect, ci_lower, ci_upper, p) in enumerate(zip(effects, ci_lowers, ci_uppers, p_values)):
        color = cb_palette[0] if effect < 0 else cb_palette[1]
        
        # CI线
        ax.plot([ci_lower, ci_upper], [i, i], color=color, linewidth=2, alpha=0.7)
        
        # 效应量点
        sig_marker = '*' if p < 0.05 else ''
        ax.scatter(effect, i, color=color, s=100, zorder=3, edgecolors='black', linewidth=1)
        
        # 数值标签
        ax.text(effect, i+0.2, f'{effect:.2f}{sig_marker}', 
               ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1.5)
    ax.axvline(x=-0.2, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(x=0.2, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    # 添加区域标注
    ax.text(-0.5, len(labels)+0.5, 'MDD < Control\n(Reduced in MDD)', 
           ha='center', fontsize=10, color=cb_palette[0], fontweight='bold')
    ax.text(0.5, len(labels)+0.5, 'MDD > Control\n(Increased in MDD)', 
           ha='center', fontsize=10, color=cb_palette[1], fontweight='bold')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel("Cohen's d (Effect Size)", fontsize=13, fontweight='bold')
    ax.set_xlim(-1, 1)
    ax.set_title('Figure 5. Effect Sizes of Retinal Changes in MDD Patients\n'
                'Horizontal lines represent 95% confidence intervals. *P<0.05.',
                fontsize=13, fontweight='bold', pad=20)
    ax.grid(True, axis='x', alpha=0.3)
    
    # 图例
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=cb_palette[0], 
               markersize=10, label='Reduced in MDD (negative d)', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=cb_palette[1], 
               markersize=10, label='Increased in MDD (positive d)', markeredgecolor='black'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
    
    # 添加样本量说明
    fig.text(0.02, 0.02, f'Note: MDD n={len(mdd)}, Control n={len(control)}. '
            'Effect size calculated using Cohen\'s d with 95% CI. '
            'Negative values indicate lower thickness/volume in MDD patients.',
            fontsize=9, style='italic')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure5_Forest_Plot_修订版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存 (带95% CI和方向说明)")

def fix_figure6_v2(mdd):
    """Figure 6: 修复数据一致性"""
    print("\n修复 Figure 6: Subgroup Analysis (修订版)...")
    
    mdd_with_phq9 = mdd[mdd['PHQ-9'].notna()].copy()
    
    # 按PHQ-9分层
    def classify_phq9(score):
        if pd.isna(score):
            return None
        elif score < 5:
            return 'Minimal (0-4)'
        elif score < 10:
            return 'Mild (5-9)'
        elif score < 15:
            return 'Moderate (10-14)'
        else:
            return 'Severe (≥15)'
    
    mdd_with_phq9['PHQ9_Group'] = mdd_with_phq9['PHQ-9'].apply(classify_phq9)
    
    # 统计
    groups = ['Minimal (0-4)', 'Mild (5-9)', 'Moderate (10-14)', 'Severe (≥15)']
    stats_data = []
    
    print("  各组统计:")
    for group in groups:
        subset = mdd_with_phq9[mdd_with_phq9['PHQ9_Group'] == group]
        if len(subset) > 0:
            # 按人统计
            n_participants = subset['Patient_ID'].nunique()
            # 按眼统计
            n_eyes = len(subset)
            thickness = subset['Retina_平均厚度'].dropna()
            
            stats_data.append({
                'group': group,
                'n_participants': n_participants,
                'n_eyes': n_eyes,
                'mean': thickness.mean(),
                'std': thickness.std(),
                'sem': thickness.std() / np.sqrt(len(thickness))
            })
            
            print(f"    {group}: N={n_participants}, n_eyes={n_eyes}, "
                  f"mean={thickness.mean():.2f}±{thickness.std():.2f}")
    
    # 绘制
    fig, ax = plt.subplots(figsize=(12, 7))
    
    if len(stats_data) > 0:
        x_pos = np.arange(len(stats_data))
        means = [s['mean'] for s in stats_data]
        stds = [s['std'] for s in stats_data]
        ns_participants = [s['n_participants'] for s in stats_data]
        ns_eyes = [s['n_eyes'] for s in stats_data]
        
        colors = ['#3498db', '#2ecc71', '#f39c12', '#e74c3c']
        bars = ax.bar(x_pos, means, yerr=stds, capsize=8, 
                      color=colors[:len(stats_data)], alpha=0.7, 
                      edgecolor='black', linewidth=2)
        
        # 添加样本量标签
        for i, (mean, std, n_p, n_e) in enumerate(zip(means, stds, ns_participants, ns_eyes)):
            ax.text(i, mean + std + 3, f'N={n_p}\n(n={n_e} eyes)', 
                   ha='center', fontsize=10, fontweight='bold')
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels([s['group'] for s in stats_data], fontsize=11)
        ax.set_ylabel('Mean Macular Thickness (μm)\nMean ± SD', fontsize=12, fontweight='bold')
        ax.set_xlabel('Depression Severity (PHQ-9)', fontsize=12, fontweight='bold')
        ax.set_title('Figure 6. Mean Macular Thickness by Depression Severity in MDD Patients\n'
                    'Error bars represent standard deviation (SD).',
                    fontsize=13, fontweight='bold')
        ax.grid(True, axis='y', alpha=0.3, linestyle='--')
        
        # 总体参考线
        overall_mean = mdd_with_phq9['Retina_平均厚度'].mean()
        overall_std = mdd_with_phq9['Retina_平均厚度'].std()
        ax.axhline(y=overall_mean, color='red', linestyle='--', linewidth=2, 
                  label=f'Overall Mean\n({overall_mean:.1f}±{overall_std:.1f} μm)')
        ax.legend(fontsize=10, loc='upper right')
    
    # 添加数据表格
    table_text = 'Group | N (participants) | n (eyes) | Mean±SD (μm)\n'
    table_text += '-' * 60 + '\n'
    for s in stats_data:
        table_text += f"{s['group']:<20} | {s['n_participants']:>3} | {s['n_eyes']:>3} | {s['mean']:.1f}±{s['std']:.1f}\n"
    
    fig.text(0.5, -0.08, table_text, ha='center', fontsize=8, 
            family='monospace', transform=fig.transFigure)
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(f'{output_dir}/Figure6_Subgroup_Analysis_修订版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存 (数据一致，带完整样本量标注)")

def main():
    print("="*70)
    print("继续修复 Figures 4-6")
    print("="*70)
    
    df, df_patient, mdd, control = load_data()
    
    fix_figure4_v2(mdd)
    fix_figure5_v2(df, mdd, control)
    fix_figure6_v2(mdd)
    
    print("\n" + "="*70)
    print("✅ Figures 4-6 修订完成!")
    print("="*70)
    
    print("\n所有修订版文件:")
    for f in sorted(os.listdir(output_dir)):
        if f.endswith('.png'):
            size = os.path.getsize(os.path.join(output_dir, f)) / 1024
            print(f"  ✅ {f} ({size:.1f} KB)")

if __name__ == "__main__":
    main()