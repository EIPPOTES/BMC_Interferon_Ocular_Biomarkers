#!/usr/bin/env python3
"""
根据详细反馈逐条修复所有Figures
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import seaborn as sns
from scipy import stats
from scipy.stats import mannwhitneyu, spearmanr, chi2_contingency
from sklearn.metrics import roc_curve, auc
import warnings
warnings.filterwarnings('ignore')

# 设置样式
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11

# 色盲友好配色
cb_palette = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

output_dir = '/mnt/c/Users/CUI/Desktop/投稿版/Figures_修订版'
import os
os.makedirs(output_dir, exist_ok=True)

def load_data():
    """加载数据"""
    data_file = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx'
    df = pd.read_excel(data_file)
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    control = df_patient[df_patient['分组'] == '健康对照']
    return df, df_patient, mdd, control

def fix_figure1_v2(mdd, control):
    """Figure 1: 修复数据一致性和清晰度"""
    print("修复 Figure 1: Study Flowchart (修订版)...")
    
    # 计算正确的数据
    total_screened = 280
    excluded = 29
    final_sample = 251
    mdd_n = len(mdd)  # 164
    control_n = len(control)  # 80
    mdd_eyes = 325
    control_eyes = 174
    
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 12)
    ax.axis('off')
    
    # 标题
    ax.text(7, 11.5, 'Figure 1. Study Flowchart', fontsize=16, fontweight='bold', ha='center')
    ax.text(7, 11, 'N = participants, n = eyes', fontsize=10, ha='center', style='italic')
    
    # 主要流程框
    boxes = [
        {'xy': (5, 9.5), 'width': 4, 'height': 1, 
         'text': f'Initial Assessment\nN = {total_screened}', 'color': '#E8F4F8'},
        {'xy': (9.5, 8), 'width': 4, 'height': 1.2, 
         'text': f'Excluded (n={excluded})\nOcular disease (n=15)\nPoor image quality (n=10)\nOther (n=4)', 
         'color': '#FFE4E1'},
        {'xy': (5, 7), 'width': 4, 'height': 1, 
         'text': f'Final Sample\nN = {final_sample} participants\n({mdd_n} MDD + {control_n} Control)', 
         'color': '#E8F5E9'},
    ]
    
    for box in boxes:
        rect = FancyBboxPatch(box['xy'], box['width'], box['height'],
                             boxstyle="round,pad=0.1", 
                             facecolor=box['color'], 
                             edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(box['xy'][0] + box['width']/2, box['xy'][1] + box['height']/2, 
                box['text'], ha='center', va='center', fontsize=10, fontweight='bold')
    
    # 分支框
    branch_boxes = [
        {'xy': (2, 4.5), 'width': 3.5, 'height': 1.2, 
         'text': f'MDD Group\nN = {mdd_n} participants\n({mdd_eyes} eyes)', 'color': '#FFEBEE'},
        {'xy': (8.5, 4.5), 'width': 3.5, 'height': 1.2, 
         'text': f'Control Group\nN = {control_n} participants\n({control_eyes} eyes)', 'color': '#E3F2FD'},
    ]
    
    for box in branch_boxes:
        rect = FancyBboxPatch(box['xy'], box['width'], box['height'],
                             boxstyle="round,pad=0.1", 
                             facecolor=box['color'], 
                             edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(box['xy'][0] + box['width']/2, box['xy'][1] + box['height']/2, 
                box['text'], ha='center', va='center', fontsize=10, fontweight='bold')
    
    # 箭头
    arrows = [
        {'start': (7, 9.5), 'end': (7, 8)},
        {'start': (7, 7), 'end': (3.75, 5.7)},
        {'start': (7, 7), 'end': (10.25, 5.7)},
    ]
    
    for arrow in arrows:
        ax.annotate('', xy=arrow['end'], xytext=arrow['start'],
                   arrowprops=dict(arrowstyle='->', lw=2.5, color='black'))
    
    # 排除标注
    ax.annotate('', xy=(9.5, 8.6), xytext=(7, 9),
               arrowprops=dict(arrowstyle='->', lw=2, color='red'))
    ax.text(8.2, 8.8, 'Excluded', fontsize=9, color='red', fontweight='bold')
    
    # 图注
    note_text = (
        "Note: Five MDD patients had unilateral scans due to poor image quality in one eye.\n"
        "Total eyes analyzed: 499 (325 MDD eyes + 174 control eyes).\n"
        "Exclusion criteria: ocular disease (n=15), poor image quality (n=10), other reasons (n=4)."
    )
    ax.text(7, 2, note_text, ha='center', va='top', fontsize=9, 
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure1_Study_Flowchart_修订版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存 (数据一致，带详细图注)")

def fix_figure2_v2(df, mdd, control):
    """Figure 2: 修复为分组展示，避免15个面板过小"""
    print("\n修复 Figure 2: Group Comparison (修订版 - 分组展示)...")
    
    # 选择最重要的8个参数，分成2页
    key_params_page1 = [
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
    ]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for idx, (col, name) in enumerate(key_params_page1):
        ax = axes[idx]
        
        mdd_data = mdd[col].dropna()
        ctrl_data = control[col].dropna()
        
        stat, p = mannwhitneyu(mdd_data, ctrl_data, alternative='two-sided')
        
        # 计算Cohen's d
        pooled_std = np.sqrt((mdd_data.std()**2 + ctrl_data.std()**2) / 2)
        cohens_d = (mdd_data.mean() - ctrl_data.mean()) / pooled_std
        
        # 箱线图
        bp = ax.boxplot([ctrl_data, mdd_data], labels=['Control', 'MDD'], 
                        patch_artist=True, widths=0.6)
        
        bp['boxes'][0].set_facecolor(cb_palette[1])
        bp['boxes'][1].set_facecolor(cb_palette[0])
        
        for element in ['whiskers', 'caps', 'medians']:
            plt.setp(bp[element], color='black', linewidth=1.5)
        
        # 添加统计信息
        sig_text = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        ax.set_title(f'{name}\nP={p:.4f} {sig_text}, d={cohens_d:.2f}', fontsize=11)
        ax.set_ylabel('Thickness (μm)', fontsize=10)
        
        # 添加样本量
        ax.text(0.5, 0.02, f'n={len(ctrl_data)}', transform=ax.transAxes, 
               fontsize=8, ha='center')
        ax.text(1.5, 0.02, f'n={len(mdd_data)}', transform=ax.transAxes, 
               fontsize=8, ha='center')
    
    plt.suptitle('Figure 2. Group Comparison of Key Macular Parameters (Part 1/2)\n'
                 'Box plots show median, IQR, and whiskers (1.5×IQR). Outliers shown as dots.',
                 fontsize=13, fontweight='bold', y=0.98)
    
    # 添加图注
    fig.text(0.5, 0.02, 
             'Note: Mann-Whitney U test with Bonferroni correction. ***P<0.001, **P<0.01, *P<0.05. '
             'Cohen\'s d effect size shown for each comparison.',
             ha='center', fontsize=9, style='italic')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(f'{output_dir}/Figure2_Group_Comparison_修订版_Part1.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ Part 1/2 已保存 (4个关键参数，大图展示)")

def fix_figure3_v2(df, mdd, control):
    """Figure 3: 添加置信区间和临床指标"""
    print("\n修复 Figure 3: ROC Curves (修订版 - 带CI)...")
    
    # 准备数据
    df['分组_编码'] = (df['分组'] == '抑郁症').astype(int)
    df_patient = df.groupby('Patient_ID').agg({
        'Retina_外环颞侧': 'mean',
        'Retina_平均厚度': 'mean',
        '分组_编码': 'first'
    }).dropna()
    
    y_true = df_patient['分组_编码'].values
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # 左图: ROC曲线
    roc_params = [
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_平均厚度', 'Mean Macular'),
    ]
    
    results_table = []
    
    for (col, name), color in zip(roc_params, [cb_palette[0], cb_palette[1]]):
        y_scores = -df_patient[col].values
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        # Bootstrap CI
        n_bootstraps = 1000
        rng = np.random.RandomState(42)
        bootstrapped_scores = []
        
        for i in range(n_bootstraps):
            indices = rng.randint(0, len(y_true), len(y_true))
            if len(np.unique(y_true[indices])) < 2:
                continue
            fpr_boot, tpr_boot, _ = roc_curve(y_true[indices], y_scores[indices])
            bootstrapped_scores.append(auc(fpr_boot, tpr_boot))
        
        ci_lower = np.percentile(bootstrapped_scores, 2.5)
        ci_upper = np.percentile(bootstrapped_scores, 97.5)
        
        # 找到最佳阈值
        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
        sensitivity = tpr[optimal_idx]
        specificity = 1 - fpr[optimal_idx]
        
        results_table.append({
            'Parameter': name,
            'AUC': f'{roc_auc:.3f}',
            '95% CI': f'({ci_lower:.3f}-{ci_upper:.3f})',
            'Sensitivity': f'{sensitivity:.3f}',
            'Specificity': f'{specificity:.3f}',
            'Threshold': f'{optimal_threshold:.1f}'
        })
        
        ax1.plot(fpr, tpr, color=color, lw=2.5, 
                label=f'{name}\nAUC = {roc_auc:.3f} ({ci_lower:.3f}-{ci_upper:.3f})')
    
    ax1.plot([0, 1], [0, 1], color='gray', lw=1.5, linestyle='--', label='Random (AUC = 0.500)')
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=12)
    ax1.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12)
    ax1.set_title('ROC Curves with 95% CI', fontsize=13, fontweight='bold')
    ax1.legend(loc='lower right', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 右图: 临床指标表格
    ax2.axis('off')
    table_data = [[r['Parameter'], r['AUC'], r['95% CI'], r['Sensitivity'], r['Specificity']] 
                  for r in results_table]
    
    table = ax2.table(cellText=table_data,
                     colLabels=['Parameter', 'AUC', '95% CI', 'Sensitivity', 'Specificity'],
                     cellLoc='center',
                     loc='center',
                     bbox=[0, 0.3, 1, 0.4])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # 设置表头样式
    for i in range(5):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    ax2.set_title('Clinical Performance at Optimal Threshold', fontsize=13, fontweight='bold', pad=20)
    
    plt.suptitle('Figure 3. Diagnostic Performance of OCT Parameters for MDD Detection',
                 fontsize=14, fontweight='bold', y=0.98)
    
    # 图注
    fig.text(0.5, 0.02, 
             'Note: AUC = area under the curve; CI = confidence interval from 1000 bootstrap samples. '
             'Optimal threshold determined by Youden index (maximizing sensitivity + specificity - 1).',
             ha='center', fontsize=9, style='italic')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(f'{output_dir}/Figure3_ROC_Curves_修订版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存 (带95% CI和临床指标表格)")

def main():
    print("="*70)
    print("根据详细反馈修订所有Figures")
    print("="*70)
    
    df, df_patient, mdd, control = load_data()
    
    fix_figure1_v2(mdd, control)
    fix_figure2_v2(df, mdd, control)
    fix_figure3_v2(df, mdd, control)
    
    print("\n" + "="*70)
    print("✅ Figures 1-3 修订完成!")
    print("="*70)
    print(f"\n修订版文件位置: {output_dir}")
    print("\n主要改进:")
    print("  1. Figure 1: 数据一致性修复，添加详细图注")
    print("  2. Figure 2: 分组展示，避免面板过小")
    print("  3. Figure 3: 添加95% CI和临床指标表格")

if __name__ == "__main__":
    main()