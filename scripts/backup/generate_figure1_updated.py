#!/usr/bin/env python3
"""
重新生成Figure 1 - 更新眼数数据
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import os

def generate_figure1_updated():
    """生成更新后的Figure 1"""
    
    output_dir = '/mnt/c/Users/CUI/Desktop/投稿版/Figures_修订版'
    
    # 更新后的数据
    total_screened = 280
    excluded = 36  # 29 + 7
    final_sample = 244
    mdd_n = 164
    control_n = 80
    mdd_eyes = 325
    control_eyes = 160
    total_eyes = 485
    
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 12)
    ax.axis('off')
    
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
    
    # 排除标注
    ax.annotate('', xy=(9.5, 8.6), xytext=(7, 9),
               arrowprops=dict(arrowstyle='->', lw=2, color='red'))
    ax.text(8.2, 8.8, 'Excluded', fontsize=9, color='red', fontweight='bold')
    
    # 图注
    note_text = (
        f"Note: Five MDD patients had unilateral scans due to poor image quality in one eye. "
        f"Seven control participants were excluded due to missing age data. "
        f"Total eyes analyzed: {total_eyes} ({mdd_eyes} MDD eyes + {control_eyes} control eyes)."
    )
    ax.text(7, 2, note_text, ha='center', va='top', fontsize=9, 
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    output_path = f'{output_dir}/Figure1_Study_Flowchart_更新版.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Figure 1已保存: {output_path}")
    print(f"\n更新数据:")
    print(f"  人数: {final_sample} ({mdd_n} MDD + {control_n} Control)")
    print(f"  眼数: {total_eyes} ({mdd_eyes} MDD + {control_eyes} Control)")

if __name__ == "__main__":
    print("="*60)
    print("生成更新后的Figure 1")
    print("="*60)
    generate_figure1_updated()
    print("\n" + "="*60)
    print("✅ Figure 1更新完成!")
    print("="*60)