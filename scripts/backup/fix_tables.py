import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import numpy as np

# 强制使用DejaVu Sans字体（无中文）
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 读取数据
df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/OCT数据_完整整合.xlsx')

# Table 1: Baseline Characteristics
def create_table1():
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.axis('off')
    
    # 计算基线特征
    participant_df = df.drop_duplicates(subset=['ID'])
    mdd_df = participant_df[participant_df['Group'] == 'MDD']
    control_df = participant_df[participant_df['Group'] == 'Control']
    
    table_data = [
        ['Characteristic', 'MDD (n={})'.format(len(mdd_df)), 'Control (n={})'.format(len(control_df)), 'P-value'],
        ['Age, years', '{:.1f} ± {:.1f}'.format(mdd_df['Age'].mean(), mdd_df['Age'].std()), 
         '{:.1f} ± {:.1f}'.format(control_df['Age'].mean(), control_df['Age'].std()), '0.XXX'],
        ['Female, n (%)', '{} ({:.1f})'.format((mdd_df['Gender']=='F').sum(), (mdd_df['Gender']=='F').mean()*100),
         '{} ({:.1f})'.format((control_df['Gender']=='F').sum(), (control_df['Gender']=='F').mean()*100), '0.XXX'],
    ]
    
    table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                     colWidths=[0.3, 0.25, 0.25, 0.2])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # 设置表头样式（三线表）
    for i in range(4):
        table[(0, i)].set_facecolor('white')
        table[(0, i)].set_text_props(fontweight='bold')
        table[(0, i)].set_edgecolor('black')
        table[(0, i)].set_linewidth(1.5)
    
    # 只保留顶部和底部横线
    for key, cell in table.get_celld().items():
        cell.set_edgecolor('none')
    
    # 添加顶部和底部线
    ax.plot([0.1, 0.9], [0.85, 0.85], 'k-', linewidth=1.5, transform=ax.transAxes)
    ax.plot([0.1, 0.9], [0.15, 0.15], 'k-', linewidth=1.5, transform=ax.transAxes)
    
    plt.title('Table 1. Baseline Characteristics of Study Participants', fontsize=12, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig('/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/Table1_Baseline_Characteristics_Fixed.png', 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print('Table 1 saved successfully!')

create_table1()
print('All tables regenerated with English-only labels!')
