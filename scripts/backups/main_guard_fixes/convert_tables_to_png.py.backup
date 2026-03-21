import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm
import warnings
warnings.filterwarnings('ignore')

# 强制重新初始化字体管理器
fm.fontManager.__init__()

# 适配跨平台中文字体
cjk_list = ['CJK', 'Han', 'CN', 'TW', 'JP']
cjk_fonts = [f.name for f in fm.fontManager.ttflist if any(s.lower() in f.name.lower() for s in cjk_list)]

plt.rcParams['font.family'] = ['DejaVu Sans'] + cjk_fonts
plt.rcParams['axes.unicode_minus'] = False

print("=" * 80)
print("将Excel表格转换为PNG图片")
print("=" * 80)

input_dir = '/mnt/c/Users/CUI/Desktop/最终修改'
output_dir = '/mnt/c/Users/CUI/Desktop/论文及图表'

# 定义要转换的表格
tables = [
    ('Table1_基线特征表.xlsx', 'Table1_基线特征表.png', 'Table 1. Participant Characteristics'),
    ('Table2_黄斑所有层_Statsmodels.xlsx', 'Table2_黄斑所有层.png', 'Table 2. Macular Layer Analysis'),
    ('Table3_黄斑综合指标_Statsmodels.xlsx', 'Table3_黄斑综合指标.png', 'Table 3. Macular Comprehensive Metrics'),
    ('Table4_视盘所有指标组间比较.xlsx', 'Table4_视盘指标.png', 'Table 4. Optic Disc Parameters'),
    ('Table5_全面相关分析.xlsx', 'Table5_相关分析.png', 'Table 5. Correlation Analysis'),
    ('Table6_全面ROC分析.xlsx', 'Table6_ROC分析.png', 'Table 6. ROC Analysis'),
    ('Table7_亚组分析_严重程度.xlsx', 'Table7_亚组分析.png', 'Table 7. Subgroup Analysis by Severity'),
    ('Table8_多因素回归_Statsmodels.xlsx', 'Table8_多因素回归.png', 'Table 8. Multivariate Regression'),
    ('Table9_双眼一致性分析.xlsx', 'Table9_双眼一致性.png', 'Table 9. Inter-eye Consistency'),
    ('Table10_机器学习模型比较.xlsx', 'Table10_机器学习模型.png', 'Table 10. Machine Learning Model Comparison'),
    ('Table11_特征重要性.xlsx', 'Table11_特征重要性.png', 'Table 11. Feature Importance'),
]

def excel_to_png(excel_path, png_path, title):
    """将Excel表格转换为PNG图片"""
    try:
        # 读取Excel
        df = pd.read_excel(excel_path)
        
        # 确定图片大小（根据列数和行数）
        n_rows = len(df) + 1  # +1 for header
        n_cols = len(df.columns)
        
        # 计算图片尺寸
        fig_width = max(12, n_cols * 2.5)
        fig_height = max(6, n_rows * 0.4 + 1)
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ax.axis('tight')
        ax.axis('off')
        
        # 准备表格数据
        cell_text = df.values.tolist()
        col_labels = df.columns.tolist()
        
        # 格式化数值（保留3位小数）
        formatted_cells = []
        for row in cell_text:
            formatted_row = []
            for cell in row:
                if isinstance(cell, float):
                    if abs(cell) < 0.001 and cell != 0:
                        formatted_row.append(f'{cell:.2e}')
                    else:
                        formatted_row.append(f'{cell:.4f}' if abs(cell) < 1 else f'{cell:.2f}')
                else:
                    formatted_row.append(str(cell))
            formatted_cells.append(formatted_row)
        
        # 创建表格
        table = ax.table(
            cellText=formatted_cells,
            colLabels=col_labels,
            loc='center',
            cellLoc='center'
        )
        
        # 设置表格样式
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)
        
        # 设置表头样式
        for i in range(len(col_labels)):
            table[(0, i)].set_facecolor('#4472C4')
            table[(0, i)].set_text_props(color='white', fontweight='bold')
        
        # 设置交替行颜色
        for i in range(1, len(formatted_cells) + 1):
            for j in range(len(col_labels)):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#E7E6E6')
                else:
                    table[(i, j)].set_facecolor('#FFFFFF')
        
        # 添加标题
        plt.title(title, fontsize=14, fontweight='bold', pad=20)
        
        # 保存图片
        plt.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.5)
        plt.close()
        
        return True
    except Exception as e:
        print(f"  错误: {e}")
        return False

# 转换所有表格
print(f"\n输入目录: {input_dir}")
print(f"输出目录: {output_dir}\n")

success_count = 0
for excel_file, png_file, title in tables:
    excel_path = f'{input_dir}/{excel_file}'
    png_path = f'{output_dir}/{png_file}'
    
    print(f"转换: {excel_file} → {png_file}")
    if excel_to_png(excel_path, png_path, title):
        print(f"  ✓ 成功")
        success_count += 1
    else:
        print(f"  ✗ 失败")

print(f"\n{'=' * 80}")
print(f"转换完成: {success_count}/{len(tables)} 个表格")
print(f"{'=' * 80}")

# 列出输出文件
print(f"\n输出文件列表:")
import os
for f in sorted(os.listdir(output_dir)):
    if f.endswith('.png'):
        size = os.path.getsize(f'{output_dir}/{f}') / 1024
        print(f"  {f} ({size:.1f} KB)")
