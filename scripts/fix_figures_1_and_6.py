#!/usr/bin/env python3
# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
修复Figure 1和Figure 6
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

def load_data():
    """加载数据"""
    data_file = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    df = pd.read_excel(data_file)
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    control = df_patient[df_patient['分组'] == '健康对照']
    return df, df_patient, mdd, control

def fix_figure1(mdd, control):
    """修复Figure 1 - 添加箭头"""
    print("修复 Figure 1: 添加流程箭头...")
    
    output_dir = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # 标题
    ax.text(5, 9.5, 'Figure 1. Study Flowchart', fontsize=16, fontweight='bold', ha='center')
    
    # 定义框的位置和大小
    boxes = [
        {'xy': (3.5, 8), 'width': 3, 'height': 0.8, 'text': 'Initial Assessment\nN = 280', 'color': 'lightblue'},
        {'xy': (3.5, 6.5), 'width': 3, 'height': 0.8, 'text': 'Excluded (n=29)\n• Ocular disease (n=15)\n• Poor image quality (n=10)\n• Other reasons (n=4)', 'color': 'lightyellow'},
        {'xy': (3.5, 5), 'width': 3, 'height': 0.8, 'text': f'Final Sample\nN = 251 ({len(mdd)} MDD + {len(control)} Control)', 'color': 'lightgreen'},
        {'xy': (1.5, 3), 'width': 2.5, 'height': 0.8, 'text': f'MDD Group\nN = {len(mdd)}', 'color': 'lightcoral'},
        {'xy': (6, 3), 'width': 2.5, 'height': 0.8, 'text': f'Control Group\nN = {len(control)}', 'color': 'lightcyan'},
        {'xy': (1.5, 1.5), 'width': 2.5, 'height': 0.6, 'text': 'OCT Scans\n325 eyes', 'color': 'wheat'},
        {'xy': (6, 1.5), 'width': 2.5, 'height': 0.6, 'text': 'OCT Scans\n174 eyes', 'color': 'wheat'},
    ]
    
    # 绘制框
    for box in boxes:
        rect = mpatches.FancyBboxPatch(box['xy'], box['width'], box['height'],
                                       boxstyle="round,pad=0.1", 
                                       facecolor=box['color'], 
                                       edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(box['xy'][0] + box['width']/2, box['xy'][1] + box['height']/2, 
                box['text'], ha='center', va='center', fontsize=10, fontweight='bold')
    
    # 添加箭头
    arrows = [
        {'start': (5, 8), 'end': (5, 6.5)},  # Initial -> Excluded (side)
        {'start': (5, 6.5), 'end': (5, 5.8)},  # Excluded -> Final
        {'start': (5, 5), 'end': (2.75, 3.8)},  # Final -> MDD
        {'start': (5, 5), 'end': (7.25, 3.8)},  # Final -> Control
        {'start': (2.75, 3), 'end': (2.75, 2.1)},  # MDD -> OCT
        {'start': (7.25, 3), 'end': (7.25, 2.1)},  # Control -> OCT
    ]
    
    for arrow in arrows:
        ax.annotate('', xy=arrow['end'], xytext=arrow['start'],
                   arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # 添加排除标记
    ax.text(6.8, 7.25, '→ Excluded', fontsize=9, color='red', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure1_Study_Flowchart_修复版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存: Figure1_Study_Flowchart_修复版.png (带箭头)")

def fix_figure6(mdd):
    """修复Figure 6 - 检查并修正数据"""
    print("\n修复 Figure 6: 检查数据...")
    
    output_dir = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    
    # 检查PHQ-9数据
    mdd_with_phq9 = mdd[mdd['PHQ-9'].notna()]
    print(f"  有PHQ-9数据的MDD患者: {len(mdd_with_phq9)}人")
    
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
    
    mdd_analysis = mdd_with_phq9.copy()
    mdd_analysis['PHQ9_Group'] = mdd_analysis['PHQ-9'].apply(classify_phq9)
    
    # 统计各组
    groups = ['Minimal (0-4)', 'Mild (5-9)', 'Moderate (10-14)', 'Severe (≥15)']
    stats = []
    
    print(f"\n  各组统计:")
    for group in groups:
        subset = mdd_analysis[mdd_analysis['PHQ9_Group'] == group]
        if len(subset) > 0:
            thickness = subset['Retina_平均厚度'].dropna()
            stats.append({
                'group': group,
                'n': len(thickness),
                'mean': thickness.mean(),
                'std': thickness.std()
            })
            print(f"    {group}: n={len(thickness)}, mean={thickness.mean():.2f}, std={thickness.std():.2f}")
        else:
            print(f"    {group}: n=0 (无数据)")
    
    # 绘制修正后的Figure 6
    fig, ax = plt.subplots(figsize=(12, 7))
    
    if len(stats) > 0:
        x_pos = np.arange(len(stats))
        means = [s['mean'] for s in stats]
        stds = [s['std'] for s in stats]
        ns = [s['n'] for s in stats]
        
        colors = ['#3498db', '#2ecc71', '#f39c12', '#e74c3c']
        bars = ax.bar(x_pos, means, yerr=stds, capsize=8, 
                      color=colors[:len(stats)], alpha=0.7, 
                      edgecolor='black', linewidth=2)
        
        # 添加样本量标签
        for i, (mean, std, n) in enumerate(zip(means, stds, ns)):
            ax.text(i, mean + std + 3, f'n={n}', ha='center', fontsize=11, fontweight='bold')
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels([s['group'] for s in stats], fontsize=11)
        ax.set_ylabel('Mean Macular Thickness (μm)', fontsize=13, fontweight='bold')
        ax.set_xlabel('Depression Severity (PHQ-9)', fontsize=13, fontweight='bold')
        ax.set_title('Figure 6. Mean Macular Thickness by Depression Severity\n(MDD Patients Only)', 
                     fontsize=14, fontweight='bold')
        ax.grid(True, axis='y', alpha=0.3, linestyle='--')
        
        # 添加均值线
        overall_mean = mdd_analysis['Retina_平均厚度'].mean()
        ax.axhline(y=overall_mean, color='red', linestyle='--', linewidth=2, 
                  label=f'Overall Mean ({overall_mean:.1f} μm)')
        ax.legend(fontsize=11)
    else:
        ax.text(0.5, 0.5, 'No PHQ-9 data available', ha='center', va='center', 
               transform=ax.transAxes, fontsize=14)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure6_Subgroup_Analysis_修复版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存: Figure6_Subgroup_Analysis_修复版.png")

def main():
    print("="*70)
    print("修复 Figure 1 和 Figure 6")
    print("="*70)
    
    df, df_patient, mdd, control = load_data()
    
    fix_figure1(mdd, control)
    fix_figure6(mdd)
    
    print("\n" + "="*70)
    print("✅ 修复完成!")
    print("="*70)
    
    print("""
修复内容:
  1. Figure 1: 添加了流程箭头，显示研究流程方向
  2. Figure 6: 检查了PHQ-9数据分布，修正了显示

修复后的文件:
  - Figure1_Study_Flowchart_修复版.png
  - Figure6_Subgroup_Analysis_修复版.png
    """)

if __name__ == "__main__":
    main()