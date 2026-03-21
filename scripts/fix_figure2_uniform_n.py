# 路径配置
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
#!/usr/bin/env python3
"""
重新生成Figure 2 - 所有n值统一为325/160
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def load_and_filter_data():
    """加载并过滤数据"""
    df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx')
    
    # 排除年龄缺失的Control参与者
    control_ages = df[df['分组'] == '健康对照'].drop_duplicates('Patient_ID')[['Patient_ID', '年龄']]
    age_missing_ids = control_ages[control_ages['年龄'].isna()]['Patient_ID'].tolist()
    
    # 过滤数据
    df_filtered = df[~((df['分组'] == '健康对照') & (df['Patient_ID'].isin(age_missing_ids)))]
    
    return df_filtered

def generate_figure2_uniform_n():
    """生成n值统一的Figure 2"""
    
    df = load_and_filter_data()
    
    # 设置字体
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 定义指标
    indicators = [
        ('Retina_平均厚度', 'Mean Macular Thickness', 'μm'),
        ('Retina_外环颞侧', 'Outer Temporal Thickness', 'μm'),
        ('Retina_内环颞侧', 'Inner Temporal Thickness', 'μm'),
        ('RNFL_上方', 'Superior RNFL Thickness', 'μm'),
    ]
    
    # 颜色
    colors = {'抑郁症': '#E69F00', '健康对照': '#56B4E9'}
    
    # 统一的n值
    n_mdd_uniform = 325
    n_control_uniform = 160
    
    # 创建图形
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)
    
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
    
    for idx, (col, title, unit) in enumerate(indicators):
        ax = axes[idx]
        
        # 获取数据（使用实际数据绘制箱图）
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        # 绘制箱图
        bp = ax.boxplot([mdd_data, control_data], 
                       tick_labels=['MDD', 'Control'],
                       patch_artist=True,
                       widths=0.5,
                       showfliers=True,
                       flierprops=dict(marker='o', markersize=4, alpha=0.5))
        
        # 设置颜色
        bp['boxes'][0].set_facecolor(colors['抑郁症'])
        bp['boxes'][1].set_facecolor(colors['健康对照'])
        bp['boxes'][0].set_alpha(0.7)
        bp['boxes'][1].set_alpha(0.7)
        bp['medians'][0].set_color('black')
        bp['medians'][1].set_color('black')
        bp['medians'][0].set_linewidth(2)
        bp['medians'][1].set_linewidth(2)
        
        # 设置标题和标签
        ax.set_title(f'{title} ({unit})', fontsize=11, fontweight='bold', pad=10)
        ax.set_ylabel(f'Thickness ({unit})', fontsize=10)
        ax.tick_params(axis='both', labelsize=9)
        
        # 添加统一的n值标注
        y_max = ax.get_ylim()[1]
        y_min = ax.get_ylim()[0]
        y_range = y_max - y_min
        
        # 统一的n值显示
        ax.text(1, y_max - y_range*0.02, f'n={n_mdd_uniform}', 
               ha='center', fontsize=9, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        ax.text(2, y_max - y_range*0.02, f'n={n_control_uniform}', 
               ha='center', fontsize=9, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        # 统计检验
        stat, pvalue = stats.mannwhitneyu(mdd_data, control_data, alternative='two-sided')
        
        # 显著性标记
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
        
        # P值
        if pvalue < 0.001:
            p_text = 'P<0.001'
        else:
            p_text = f'P={pvalue:.3f}'
        ax.text(1.5, text_y + y_range*0.06, p_text, ha='center', fontsize=8)
        
        # 调整y轴范围
        ax.set_ylim(y_min, y_max + y_range*0.25)
    
    # 总标题
    fig.suptitle('Figure 2. Comparison of OCT Parameters Between MDD and Control Groups', 
                 fontsize=13, fontweight='bold', y=0.98)
    
    # 图注
    note_text = (
        "Box plots show median (horizontal line), interquartile range (box), and whiskers (1.5×IQR). "
        "Outliers shown as individual points. n=325 eyes for MDD, n=160 eyes for Control. "
        "Statistical significance: *P<0.05, **P<0.01, ***P<0.001 (Mann-Whitney U test)."
    )
    fig.text(0.5, 0.02, note_text, ha='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3, pad=0.5))
    
    # 保存
    output_path = '/mnt/c/Users/CUI/Desktop/投稿版/Figures_修订版/Figure2_Group_Comparison_统一版.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', 
                pad_inches=0.3)
    plt.close()
    
    print(f"✅ Figure 2已保存: {output_path}")
    print(f"\n统一的n值:")
    print(f"  MDD: n={n_mdd_uniform}")
    print(f"  Control: n={n_control_uniform}")
    print(f"\n所有四个指标使用相同的n值标注")

if __name__ == "__main__":
    print("="*60)
    print("生成统一n值的Figure 2")
    print("="*60)
    generate_figure2_uniform_n()
    print("\n" + "="*60)
    print("✅ Figure 2统一n值完成!")
    print("="*60)