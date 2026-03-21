import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import matplotlib.patches as patches
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

# 创建图形 - 使用text方式绘制标准三线表
fig, ax = plt.subplots(figsize=(14, 10))
ax.set_xlim(0, 14)
ax.set_ylim(0, 16)
ax.axis('off')

# 标题
ax.text(7, 15, 'Table 3. Comparison of Optic Disc Parameters Between Groups', 
        ha='center', va='center', fontsize=14, fontweight='bold')

# 列宽设置
col_x = [0.5, 4.0, 7.0, 10.0, 12.5]
col_widths = [3.0, 2.5, 2.5, 1.8, 1.5]

# 表头
headers = ['Parameter', 'MDD Group', 'Control Group', 'P-value', "Cohen's d"]
for i, (header, x) in enumerate(zip(headers, col_x)):
    ax.text(x + col_widths[i]/2, 13.5, header, ha='center', va='center', 
            fontsize=11, fontweight='bold')

# 顶部横线
ax.plot([0.3, 13.7], [14, 14], 'k-', linewidth=1.5)
# 表头下方横线
ax.plot([0.3, 13.7], [13, 13], 'k-', linewidth=1.0)

# 数据行
row_height = 0.8
start_y = 12.2
for idx, row in df.iterrows():
    y = start_y - idx * row_height
    for i, (col, x) in enumerate(zip(df.columns, col_x)):
        ax.text(x + col_widths[i]/2, y, row[col], ha='center', va='center', fontsize=10)

# 底部横线
end_y = start_y - (len(df) - 1) * row_height - 0.4
ax.plot([0.3, 13.7], [end_y, end_y], 'k-', linewidth=1.5)

# 添加注释
ax.text(7, end_y - 0.8, 
        'Data are mean ± SD. P-values from Mann-Whitney U tests with FDR correction.',
        ha='center', va='center', fontsize=9, style='italic')

plt.tight_layout()
plt.savefig('/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/Table3_Optic_Disc_ThreeLine.png', 
            dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.3)
plt.savefig('/root/.openclaw/workspace/tables/Table3_Optic_Disc_ThreeLine.png', 
            dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.3)
print("Table 3 (standard three-line format) generated successfully!")
