import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import numpy as np

# 强制重新初始化字体管理器
fm.fontManager.__init__()

# 适配中文字体
cjk_list = ['CJK', 'Han', 'CN', 'TW', 'JP']
cjk_fonts = [f.name for f in fm.fontManager.ttflist if any(s.lower() in f.name.lower() for s in cjk_list)]

plt.rcParams['font.family'] = ['DejaVu Sans'] + cjk_fonts
plt.rcParams['axes.unicode_minus'] = False

# Table 3 数据 - 视盘参数
optic_disc_data = {
    'Parameter': [
        'RNFL Total',
        'RNFL Superior',
        'RNFL Temporal', 
        'RNFL Nasal',
        'RNFL Inferior',
        'Disc Area',
        'Cup Area',
        'Rim Area',
        'Cup Volume',
        'Rim Volume',
        'C/D Area Ratio',
        'Linear C/D Ratio',
        'Vertical C/D Ratio'
    ],
    'MDD Group': [
        '105.89±17.43',
        '86.32±24.51',
        '132.59±27.04',
        '67.94±20.31',
        '136.27±29.43',
        '2.06±1.05',
        '0.67±0.66',
        '1.39±0.64',
        '0.12±0.14',
        '0.25±0.16',
        '0.30±0.19',
        '0.51±0.21',
        '0.48±0.20'
    ],
    'Control Group': [
        '109.67±13.28',
        '93.53±20.44',
        '137.80±21.80',
        '65.09±19.43',
        '141.61±25.65',
        '1.97±0.83',
        '0.54±0.54',
        '1.43±0.58',
        '0.10±0.14',
        '0.30±0.24',
        '0.25±0.18',
        '0.46±0.21',
        '0.43±0.21'
    ],
    'P-value': [
        '0.047',
        '<0.001',
        '0.027',
        '0.065',
        '0.094',
        '0.627',
        '0.010',
        '0.202',
        '0.024',
        '0.002',
        '0.008',
        '0.008',
        '0.003'
    ],
    "Cohen's d": [
        '-0.23',
        '-0.31',
        '-0.21',
        '0.14',
        '-0.19',
        '0.10',
        '0.22',
        '-0.07',
        '0.15',
        '-0.30',
        '0.25',
        '0.24',
        '0.26'
    ]
}

df = pd.DataFrame(optic_disc_data)

# 创建图形
fig, ax = plt.subplots(figsize=(14, 8))
ax.axis('off')

# 标题
fig.suptitle('Table 3. Comparison of Optic Disc Parameters Between Groups', 
             fontsize=14, fontweight='bold', y=0.96)

# 准备表格数据
table_data = [df.columns.tolist()] + df.values.tolist()

# 创建表格
table = ax.table(
    cellText=df.values,
    colLabels=df.columns,
    cellLoc='center',
    loc='center',
    bbox=[0.05, 0.15, 0.9, 0.75]
)

# 设置表格样式 - 三线表
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)

# 设置表头样式
for i in range(len(df.columns)):
    cell = table[(0, i)]
    cell.set_facecolor('white')
    cell.set_text_props(fontweight='bold')
    cell.set_edgecolor('black')
    cell.set_linewidth(1.5)

# 设置数据行样式
for i in range(1, len(df) + 1):
    for j in range(len(df.columns)):
        cell = table[(i, j)]
        cell.set_facecolor('white')
        cell.set_edgecolor('none')

# 添加底部横线
for j in range(len(df.columns)):
    table[(len(df), j)].set_edgecolor('black')
    table[(len(df), j)].set_linewidth(1.5)

# 添加注释
fig.text(0.5, 0.08, 
         'Data are mean ± SD. P-values from Mann-Whitney U tests with FDR correction.',
         ha='center', fontsize=9, style='italic')

plt.tight_layout()
plt.savefig('/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/Table3_Optic_Disc_ThreeLine.png', 
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('/root/.openclaw/workspace/tables/Table3_Optic_Disc_ThreeLine.png', 
            dpi=300, bbox_inches='tight', facecolor='white')
print("Table 3 generated successfully!")
