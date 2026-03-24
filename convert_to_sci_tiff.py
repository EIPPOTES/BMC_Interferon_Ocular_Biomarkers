#!/usr/bin/env python3
"""
SCI期刊标准TIFF转换脚本
基于高分辨率PNG生成符合期刊要求的TIFF文件
"""

import os
import sys
from PIL import Image

def convert_png_to_sci_tiff(png_path, output_path, target_dpi=300):
    """将PNG转换为SCI标准TIFF"""
    try:
        # 打开PNG
        img = Image.open(png_path)
        
        # 转换为RGB模式(去除Alpha通道，更适合印刷)
        if img.mode in ('RGBA', 'LA', 'P'):
            # 创建白色背景
            background = Image.new('RGB', img.size, (255, 255, 255))
            if img.mode == 'P':
                img = img.convert('RGBA')
            if img.mode in ('RGBA', 'LA'):
                background.paste(img, mask=img.split()[-1] if img.mode in ('RGBA', 'LA') else None)
                img = background
            else:
                img = img.convert('RGB')
        elif img.mode != 'RGB':
            img = img.convert('RGB')
        
        # 计算目标尺寸以保持高DPI
        # 假设单栏图8.5cm，双栏17.5cm
        current_width_px = img.size[0]
        target_width_inches = 6.69  # 17cm = 6.69 inches
        target_width_px = int(target_width_inches * target_dpi)
        
        # 如果当前分辨率足够高，保持原尺寸
        if current_width_px >= target_width_px:
            # 只调整DPI元数据，不缩放
            img.save(
                output_path,
                format='TIFF',
                dpi=(target_dpi, target_dpi),
                compression='tiff_lzw'  # LZW压缩，期刊常用
            )
            print(f"  ✅ 已转换 (保持原尺寸，设置DPI为{target_dpi})")
        else:
            # 需要放大
            ratio = target_width_px / current_width_px
            new_height = int(img.size[1] * ratio)
            img_resized = img.resize((target_width_px, new_height), Image.Resampling.LANCZOS)
            img_resized.save(
                output_path,
                format='TIFF',
                dpi=(target_dpi, target_dpi),
                compression='tiff_lzw'
            )
            print(f"  ✅ 已转换 (放大至{target_dpi} DPI)")
        
        return True
    except Exception as e:
        print(f"  ❌ 转换失败: {e}")
        return False

def batch_convert_figures(figures_dir):
    """批量转换所有图表"""
    print("="*60)
    print("SCI期刊图表批量转换工具")
    print("="*60)
    print()
    
    # 需要转换的图表
    figure_mappings = [
        ("Figure1_Study_Flowchart_463eyes.png", "Figure1_Study_Flowchart_463eyes_SCI.tiff"),
        ("Figure2_Group_comparison_forest_plot.png", "Figure2_Group_comparison_forest_plot_SCI.tiff"),
        ("Figure3_ROC_curves_comparison.png", "Figure3_ROC_curves_comparison_SCI.tiff"),
        ("Figure4_Correlation_scatter_plots.png", "Figure4_Correlation_scatter_plots_SCI.tiff"),
        ("Figure5_Effect_size_forest_plot.png", "Figure5_Effect_size_forest_plot_SCI.tiff"),
        ("Figure6_Subgroup_analysis_forest_plot.png", "Figure6_Subgroup_analysis_forest_plot_SCI.tiff"),
    ]
    
    converted = []
    failed = []
    
    for png_file, tiff_file in figure_mappings:
        png_path = os.path.join(figures_dir, png_file)
        tiff_path = os.path.join(figures_dir, tiff_file)
        
        if os.path.exists(png_path):
            print(f"📊 转换 {png_file}...")
            if convert_png_to_sci_tiff(png_path, tiff_path, target_dpi=300):
                converted.append(tiff_file)
            else:
                failed.append(png_file)
        else:
            print(f"⚠️  跳过 {png_file} (文件不存在)")
            failed.append(png_file)
    
    print()
    print("="*60)
    print("【转换结果】")
    print(f"✅ 成功: {len(converted)} 个文件")
    for f in converted:
        print(f"   - {f}")
    if failed:
        print(f"❌ 失败: {len(failed)} 个文件")
        for f in failed:
            print(f"   - {f}")
    
    print()
    print("【使用说明】")
    print("1. 新生成的 *_SCI.tiff 文件符合SCI期刊要求")
    print("2. DPI设置为300，适合期刊印刷")
    print("3. 使用LZW压缩，平衡质量和文件大小")
    print("4. 投稿时优先使用这些文件")
    print("="*60)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        figures_dir = "/mnt/c/Users/CUI/Desktop/final/02_Figures"
    else:
        figures_dir = sys.argv[1]
    
    if os.path.exists(figures_dir):
        batch_convert_figures(figures_dir)
    else:
        print(f"❌ 目录不存在: {figures_dir}")