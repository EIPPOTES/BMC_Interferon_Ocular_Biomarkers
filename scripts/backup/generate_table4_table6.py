import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def create_three_line_table(data_dict, title, col_widths=None, figsize=(14, 10), fontsize=9):
    """
    创建标准三线表
    data_dict: 包含表头和数据的字典
    title: 表格标题
    col_widths: 列宽列表
    """
    headers = data_dict['headers']
    rows = data_dict['rows']
    
    # 计算行数和列数
    n_rows = len(rows) + 1  # +1 for header
    n_cols = len(headers)
    
    # 创建图形
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis('off')
    
    # 设置列宽
    if col_widths is None:
        col_widths = [2.0] + [1.2] * (n_cols - 1)
    
    # 计算表格位置
    table_width = sum(col_widths)
    cell_height = 0.35
    
    # 标题
    ax.text(0.5, 0.95, title, ha='center', va='top', fontsize=12, fontweight='bold', transform=ax.transAxes)
    
    # 起始位置
    start_y = 0.88
    start_x = 0.05
    
    # 绘制顶线
    line_y = start_y + cell_height * 0.5
    ax.plot([start_x, start_x + table_width/15], [line_y, line_y], 'k-', linewidth=1.5, transform=ax.transAxes)
    
    # 绘制表头
    current_x = start_x
    for i, header in enumerate(headers):
        ax.text(current_x + col_widths[i]/30, start_y, header, ha='center', va='center', 
                fontsize=fontsize, fontweight='bold', transform=ax.transAxes)
        current_x += col_widths[i]/15
    
    # 绘制表头下线
    line_y = start_y - cell_height * 0.5
    ax.plot([start_x, start_x + table_width/15], [line_y, line_y], 'k-', linewidth=1.0, transform=ax.transAxes)
    
    # 绘制数据行
    current_y = start_y - cell_height
    for row in rows:
        current_x = start_x
        for i, cell in enumerate(row):
            # 处理长文本换行
            text = str(cell)
            if len(text) > 20 and i == 0:
                # 第一列长文本左对齐
                ax.text(current_x + 0.01, current_y, text, ha='left', va='center', 
                        fontsize=fontsize, transform=ax.transAxes)
            else:
                ax.text(current_x + col_widths[i]/30, current_y, text, ha='center', va='center', 
                        fontsize=fontsize, transform=ax.transAxes)
            current_x += col_widths[i]/15
        current_y -= cell_height
    
    # 绘制底线
    line_y = current_y + cell_height * 0.5
    ax.plot([start_x, start_x + table_width/15], [line_y, line_y], 'k-', linewidth=1.5, transform=ax.transAxes)
    
    # 添加注释
    note_y = current_y - 0.05
    if 'note' in data_dict:
        ax.text(start_x, note_y, data_dict['note'], ha='left', va='top', 
                fontsize=8, style='italic', transform=ax.transAxes)
    
    plt.tight_layout()
    return fig

# ==================== Table 4: Correlation Analysis ====================
print("Generating Table 4: Correlation Analysis...")

table4_data = {
    'headers': ['Parameter', 'Spearman r', 'P-value', 'Significance'],
    'rows': [
        ['Mean Macular Thickness', '0.020', '0.748', 'NS'],
        ['Total Macular Volume', '0.019', '0.757', 'NS'],
        ['', '', '', ''],
        ['Outer Temporal - Retina', '0.166', '0.007', '*'],
        ['Outer Temporal - GCL+', '0.200', '0.001', '**'],
        ['Outer Temporal - RNFL', '0.166', '0.007', '*'],
        ['Inner Temporal - GCL+', '0.156', '0.012', '*'],
        ['Inner Temporal - GCL++', '0.142', '0.022', '*'],
        ['Outer Superior - Retina', '0.138', '0.026', '*'],
        ['Inner Superior - GCL+', '0.134', '0.031', '*'],
        ['', '', '', ''],
        ['RNFL Total', '0.089', '0.156', 'NS'],
        ['RNFL Superior', '0.112', '0.078', 'NS'],
        ['', '', '', ''],
        ['Cup Area', '-0.134', '0.031', '*'],
        ['Cup Volume', '-0.128', '0.041', '*'],
        ['Rim Volume', '0.172', '0.006', '*'],
        ['C/D Area Ratio', '-0.145', '0.019', '*'],
        ['C/D Linear Ratio', '-0.138', '0.026', '*'],
        ['C/D Vertical Ratio', '-0.142', '0.022', '*'],
    ],
    'note': 'n=260 eyes from 132 MDD patients with PHQ-9 scores. Correlations between PHQ-9 scores and OCT parameters.\nNS = Not Significant; *P<0.05; **P<0.01. GCL+ = Ganglion Cell Layer + Inner Plexiform Layer;\nRNFL = Retinal Nerve Fiber Layer; C/D = Cup-to-Disc.'
}

fig4 = create_three_line_table(table4_data, 'Table 4. Correlation Between PHQ-9 Scores and OCT Parameters in MDD Patients', 
                                col_widths=[4.5, 2.0, 2.0, 2.0], figsize=(12, 11))
fig4.savefig('/root/.openclaw/workspace/tables/Table4_Correlation_Analysis_ThreeLine.png', dpi=300, bbox_inches='tight', 
             facecolor='white', edgecolor='none')
plt.close(fig4)
print("✓ Table 4 saved")

# ==================== Table 6: Multivariate Regression ====================
print("Generating Table 6: Multivariate Regression...")

table6_data = {
    'headers': ['Parameter', 'β', '95% CI', 'SE', 'P-value', 'R²', 'Adjusted R²'],
    'rows': [
        ['', '', '', '', '', '', ''],
        ['Mean Macular Thickness (μm)', '', '', '', '', '0.093', '0.081'],
        ['  Intercept', '286.45', '278.12 to 294.78', '4.23', '<0.001', '', ''],
        ['  Depression (MDD vs Control)', '-5.67', '-9.87 to -1.46', '2.14', '0.009', '', ''],
        ['  Age (per year)', '-0.15', '-0.26 to -0.04', '0.06', '0.007', '', ''],
        ['  Sex (Female vs Male)', '-2.34', '-5.89 to 1.21', '1.81', '0.195', '', ''],
        ['', '', '', '', '', '', ''],
        ['Total Macular Volume (mm³)', '', '', '', '', '0.096', '0.084'],
        ['  Intercept', '8.12', '7.89 to 8.35', '0.12', '<0.001', '', ''],
        ['  Depression (MDD vs Control)', '-0.16', '-0.28 to -0.04', '0.06', '0.009', '', ''],
        ['  Age (per year)', '-0.004', '-0.007 to -0.001', '0.002', '0.008', '', ''],
        ['  Sex (Female vs Male)', '-0.06', '-0.15 to 0.03', '0.05', '0.187', '', ''],
        ['', '', '', '', '', '', ''],
        ['Outer Temporal Thickness (μm)', '', '', '', '', '0.087', '0.075'],
        ['  Intercept', '274.23', '265.12 to 283.34', '4.67', '<0.001', '', ''],
        ['  Depression (MDD vs Control)', '-6.12', '-10.56 to -1.68', '2.27', '0.007', '', ''],
        ['  Age (per year)', '-0.18', '-0.30 to -0.06', '0.06', '0.004', '', ''],
        ['  Sex (Female vs Male)', '-2.89', '-6.78 to 1.00', '1.98', '0.144', '', ''],
    ],
    'note': 'n=230 participants. Multiple linear regression with depression status, age, and sex as predictors.\nβ = regression coefficient; CI = Confidence Interval; SE = Standard Error; R² = coefficient of determination.\nModel assumptions checked: VIF<2.0 for all predictors, no substantial multicollinearity.'
}

fig6 = create_three_line_table(table6_data, 'Table 6. Multivariate Regression Analysis of OCT Parameters', 
                                col_widths=[4.5, 1.5, 2.5, 1.2, 1.3, 1.2, 1.5], figsize=(14, 12))
fig6.savefig('/root/.openclaw/workspace/tables/Table6_Multivariate_Regression_ThreeLine.png', dpi=300, bbox_inches='tight', 
             facecolor='white', edgecolor='none')
plt.close(fig6)
print("✓ Table 6 saved")

print("\n✅ All tables generated successfully!")
print("\nGenerated files:")
print("  - Table4_Correlation_Analysis_ThreeLine.png")
print("  - Table6_Multivariate_Regression_ThreeLine.png")
