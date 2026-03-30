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

# Table 5 数据 - ROC分析 (Top 12)
roc_data = {
    'Parameter': [
        'Outer Temporal',
        'Inner Temporal',
        'Outer Superior',
        'Total Volume',
        'Mean Thickness',
        'Inner Inferior',
        'Inner Superior',
        'GCL+ Inner Temporal',
        'RNFL Superior',
        'GCL++ Inner Temporal',
        'GCL++ Inner Superior',
        'GCL+ Inner Superior'
    ],
    'AUC': [
        '0.650',
        '0.638',
        '0.632',
        '0.629',
        '0.628',
        '0.624',
        '0.623',
        '0.622',
        '0.621',
        '0.616',
        '0.612',
        '0.609'
    ],
    'Sensitivity': [
        '0.360',
        '0.524',
        '0.713',
        '0.720',
        '0.720',
        '0.671',
        '0.518',
        '0.463',
        '0.531',
        '0.439',
        '0.384',
        '0.707'
    ],
    'Specificity': [
        '0.885',
        '0.713',
        '0.552',
        '0.552',
        '0.552',
        '0.552',
        '0.701',
        '0.793',
        '0.721',
        '0.793',
        '0.828',
        '0.506'
    ],
    'Cut-off': [
        '265.5',
        '304.7',
        '257.8',
        '7.9',
        '278.9',
        '306.1',
        '290.3',
        '115.7',
        '84.4',
        '87.7',
        '83.1',
        '109.2'
    ]
}

df = pd.DataFrame(roc_data)

# 创建图形 - 使用text方式绘制标准三线表
fig, ax = plt.subplots(figsize=(14, 10))
ax.set_xlim(0, 14)
ax.set_ylim(0, 16)
ax.axis('off')

# 标题
ax.text(7, 15, 'Table 5. ROC Analysis for Diagnostic Performance (Top 12)', 
        ha='center', va='center', fontsize=14, fontweight='bold')
ax.text(7, 14.3, 'Diagnostic Performance of OCT Parameters', 
        ha='center', va='center', fontsize=11, style='italic')

# 列宽设置
col_x = [0.5, 4.5, 7.0, 9.5, 12.0]
col_widths = [3.5, 2.0, 2.0, 2.0, 1.5]

# 表头
headers = ['Parameter', 'AUC', 'Sensitivity', 'Specificity', 'Cut-off']
for i, (header, x) in enumerate(zip(headers, col_x)):
    ax.text(x + col_widths[i]/2, 13.0, header, ha='center', va='center', 
            fontsize=11, fontweight='bold')

# 顶部横线
ax.plot([0.3, 13.7], [13.5, 13.5], 'k-', linewidth=1.5)
# 表头下方横线
ax.plot([0.3, 13.7], [12.5, 12.5], 'k-', linewidth=1.0)

# 数据行
row_height = 0.7
start_y = 11.7
for idx, row in df.iterrows():
    y = start_y - idx * row_height
    for i, (col, x) in enumerate(zip(df.columns, col_x)):
        ax.text(x + col_widths[i]/2, y, row[col], ha='center', va='center', fontsize=10)

# 底部横线
end_y = start_y - (len(df) - 1) * row_height - 0.35
ax.plot([0.3, 13.7], [end_y, end_y], 'k-', linewidth=1.5)

# 添加注释
ax.text(7, end_y - 0.6, 
        'AUC: Area Under the Curve. Cut-off values in μm for thickness measurements.',
        ha='center', va='center', fontsize=9, style='italic')

plt.tight_layout()
plt.savefig('/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/Table5_ROC_Analysis_ThreeLine.png', 
            dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.3)
plt.savefig('/root/.openclaw/workspace/tables/Table5_ROC_Analysis_ThreeLine.png', 
            dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.3)
print("Table 5 (corrected) generated successfully!")
