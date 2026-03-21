#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
生成选项B的所有图表（基于485眼数据）
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
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

def generate_figure1():
    """生成Figure 1（485眼版本）"""
    
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 12)
    ax.axis('off')
    
    # 485眼数据
    total_screened = 280
    excluded = 36
    final_sample = 244
    mdd_n = 164
    control_n = 80
    mdd_eyes = 325
    control_eyes = 160
    total_eyes = 485
    
    # 标题
    ax.text(7, 11.5, 'Figure 1. Study Flowchart', fontsize=16, fontweight='bold', ha='center')
    ax.text(7, 11, f'N = participants, n = eyes (Total: {final_sample} participants, {total_eyes} eyes)', 
            fontsize=10, ha='center', style='italic')
    
    # 流程框
    boxes = [
        {'xy': (5, 9.5), 'width': 4, 'height': 1, 
         'text': f'Initial Assessment\nN = {total_screened}', 'color': '#E8F4F8'},
        {'xy': (9.5, 8), 'width': 4, 'height': 1.2, 
         'text': f'Excluded (n={excluded})\n• Ocular disease (n=15)\n• Poor image quality (n=10)\n• Missing age (n=7)\n• Other (n=4)', 
         'color': '#FFE4E1'},
        {'xy': (5, 7), 'width': 4, 'height': 1, 
         'text': f'Final Sample\nN = {final_sample} ({mdd_n} MDD + {control_n} Control)\n{total_eyes} eyes ({mdd_eyes} MDD + {control_eyes} Control)', 
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
        {'xy': (1.5, 4.5), 'width': 3.5, 'height': 1.2, 
         'text': f'MDD Group\nN = {mdd_n} participants\n{mdd_eyes} eyes', 'color': '#FFEBEE'},
        {'xy': (8.5, 4.5), 'width': 3.5, 'height': 1.2, 
         'text': f'Control Group\nN = {control_n} participants\n{control_eyes} eyes', 'color': '#E3F2FD'},
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
        {'start': (7, 7), 'end': (3.25, 5.7)},
        {'start': (7, 7), 'end': (10.25, 5.7)},
    ]
    
    for arrow in arrows:
        ax.annotate('', xy=arrow['end'], xytext=arrow['start'],
                   arrowprops=dict(arrowstyle='->', lw=2.5, color='black'))
    
    # 图注
    note_text = (
        f"Note: Seven control participants were excluded due to missing age data. "
        f"Total eyes analyzed: {total_eyes} ({mdd_eyes} MDD eyes + {control_eyes} control eyes)."
    )
    ax.text(7, 2, note_text, ha='center', va='top', fontsize=9, 
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure1_Study_Flowchart_485eyes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Figure 1已生成 (485眼版本)")

def generate_figure2():
    """生成Figure 2（485眼版本）"""
    
    df = load_data()
    
    indicators = [
        ('Retina_平均厚度', 'Mean Macular\nThickness', 'μm'),
        ('Retina_外环颞侧', 'Outer Temporal\nThickness', 'μm'),
        ('Retina_内环颞侧', 'Inner Temporal\nThickness', 'μm'),
        ('RNFL_上方', 'Superior RNFL\nThickness', 'μm'),
    ]
    
    colors = {'抑郁症': '#E69F00', '健康对照': '#56B4E9'}
    n_mdd = 325
    n_control = 160
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
    
    for idx, (col, title, unit) in enumerate(indicators):
        ax = axes[idx]
        
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        bp = ax.boxplot([mdd_data, control_data], 
                       tick_labels=['MDD', 'Control'],
                       patch_artist=True, widths=0.5, showfliers=True)
        
        bp['boxes'][0].set_facecolor(colors['抑郁症'])
        bp['boxes'][1].set_facecolor(colors['健康对照'])
        bp['boxes'][0].set_alpha(0.7)
        bp['boxes'][1].set_alpha(0.7)
        
        ax.set_title(f'{title} ({unit})', fontsize=11, fontweight='bold', pad=10)
        ax.set_ylabel(f'Thickness ({unit})', fontsize=10)
        
        y_max = ax.get_ylim()[1]
        y_range = y_max - ax.get_ylim()[0]
        
        ax.text(1, y_max - y_range*0.02, f'n={n_mdd}', ha='center', fontsize=9, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        ax.text(2, y_max - y_range*0.02, f'n={n_control}', ha='center', fontsize=9, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        stat, pvalue = stats.mannwhitneyu(mdd_data, control_data, alternative='two-sided')
        
        data_max = max(mdd_data.max(), control_data.max())
        line_y = data_max + y_range * 0.08
        text_y = data_max + y_range * 0.12
        
        ax.plot([1, 1.5, 2], [line_y, line_y, line_y], 'k-', linewidth=1)
        
        if pvalue < 0.001:
            sig_text = '***'
        elif pvalue < 0.01:
            sig_text = '**'
        elif pvalue < 0.05:
            sig_text = '*'
        else:
            sig_text = 'ns'
        
        ax.text(1.5, text_y, sig_text, ha='center', fontsize=12, fontweight='bold')
        
        if pvalue < 0.001:
            p_text = 'P<0.001'
        else:
            p_text = f'P={pvalue:.3f}'
        ax.text(1.5, text_y + y_range*0.06, p_text, ha='center', fontsize=8)
        
        ax.set_ylim(ax.get_ylim()[0], y_max + y_range*0.25)
    
    fig.suptitle('Figure 2. Comparison of OCT Parameters Between MDD and Control Groups', 
                 fontsize=13, fontweight='bold', y=0.98)
    
    note_text = (
        "Box plots show median (horizontal line), interquartile range (box), and whiskers (1.5×IQR). "
        "n=325 eyes for MDD, n=160 eyes for Control. "
        "Statistical significance: *P<0.05, **P<0.01, ***P<0.001 (Mann-Whitney U test)."
    )
    fig.text(0.5, 0.02, note_text, ha='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3, pad=0.5))
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(f'{OUTPUT_DIR}/Figure2_Group_Comparison_485eyes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Figure 2已生成 (485眼版本)")

if __name__ == "__main__":
    print("="*70)
    print("生成选项B图表（485眼版本）")
    print("="*70)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    generate_figure1()
    generate_figure2()
    
    print("\n" + "="*70)
    print("✅ 所有图表已生成!")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*70)