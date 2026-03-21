# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm
import warnings

def main():
    """主函数，包装原有执行代码"""
    warnings.filterwarnings('ignore')

    print("=" * 80)
    print("修复中文字体问题，重新生成表格图片")
    print("=" * 80)

    # 查找系统中文字体
    fm.fontManager.__init__()

    # 列出所有可用字体
    all_fonts = [f.name for f in fm.fontManager.ttflist]
    print(f"\n系统中共有 {len(all_fonts)} 个字体")

    # 查找中文字体
    chinese_fonts = []
    for font in all_fonts:
        if any(keyword in font.lower() for keyword in ['cjk', 'chinese', 'hei', 'song', 'ming', 'kai', 'fang', 'yuan', 'noto', 'source', 'wqy', 'micro', 'yahei', 'simsun', 'simhei', 'arial unicode']):
            chinese_fonts.append(font)

    print(f"找到 {len(chinese_fonts)} 个可能的中文字体:")
    for f in chinese_fonts[:20]:
        print(f"  - {f}")

    # 尝试设置中文字体
    if chinese_fonts:
        # 优先使用的中文字体
        preferred = ['WenQuanYi Micro Hei', 'Noto Sans CJK SC', 'Source Han Sans CN', 'SimHei', 'Microsoft YaHei', 'SimSun', 'Arial Unicode MS']
        selected_font = None
        for font in preferred:
            if font in chinese_fonts:
                selected_font = font
                break
        if not selected_font:
            selected_font = chinese_fonts[0]

        print(f"\n选择字体: {selected_font}")
        plt.rcParams['font.family'] = [selected_font, 'DejaVu Sans']
    else:
        print("\n警告: 未找到中文字体，使用默认字体")
        plt.rcParams['font.family'] = ['DejaVu Sans']

    plt.rcParams['axes.unicode_minus'] = False

    input_dir = '/mnt/c/Users/CUI/Desktop/最终修改'
    output_dir = 'str(/mnt/c/Users/CUI/Desktop/论文及图表)'

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
            fig_width = max(14, n_cols * 2.2)
            fig_height = max(6, n_rows * 0.35 + 1.5)

            # 创建图形
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            ax.axis('tight')
            ax.axis('off')

            # 准备表格数据
            cell_text = df.values.tolist()
            col_labels = df.columns.tolist()

            # 格式化数值
            formatted_cells = []
            for row in cell_text:
                formatted_row = []
                for cell in row:
                    if isinstance(cell, float):
                        if abs(cell) < 0.001 and cell != 0:
                            formatted_row.append(f'{cell:.2e}')
                        elif abs(cell) < 0.01:
                            formatted_row.append(f'{cell:.4f}')
                        elif abs(cell) < 1:
                            formatted_row.append(f'{cell:.3f}')
                        else:
                            formatted_row.append(f'{cell:.2f}')
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
            table.set_fontsize(8)
            table.scale(1, 1.8)

            # 设置表头样式
            for i in range(len(col_labels)):
                table[(0, i)].set_facecolor('#4472C4')
                table[(0, i)].set_text_props(color='white', fontweight='bold', fontsize=9)

            # 设置交替行颜色
            for i in range(1, len(formatted_cells) + 1):
                for j in range(len(col_labels)):
                    if i % 2 == 0:
                        table[(i, j)].set_facecolor('#E7E6E6')
                    else:
                        table[(i, j)].set_facecolor('#FFFFFF')

            # 添加标题
            plt.title(title, fontsize=12, fontweight='bold', pad=15)

            # 保存图片
            plt.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.3)
            plt.close()

            return True
        except Exception as e:
            print(f"  错误: {e}")
            import traceback
            traceback.print_exc()
            return False

    # 转换所有表格
    print(f"\n{'=' * 80}")
    print(f"开始转换表格...")
    print(f"{'=' * 80}")

    success_count = 0
    for excel_file, png_file, title in tables:
        excel_path = f'{input_dir}/{excel_file}'
        png_path = f'{output_dir}/{png_file}'

        print(f"\n转换: {excel_file}")
        if excel_to_png(excel_path, png_path, title):
            print(f"  ✓ 成功 → {png_file}")
            success_count += 1
        else:
            print(f"  ✗ 失败")

    print(f"\n{'=' * 80}")
    print(f"转换完成: {success_count}/{len(tables)} 个表格")
    print(f"{'=' * 80}")



if __name__ == "__main__":
    main()