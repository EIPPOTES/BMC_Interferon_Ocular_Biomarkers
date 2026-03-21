#!/usr/bin/env python3
# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
检查所有Figures是否符合期刊标准
"""

from PIL import Image
import os
import sys

def check_figure(image_path):
    """检查单个图像文件"""
    try:
        img = Image.open(image_path)
        width, height = img.size
        
        # 期刊标准检查
        issues = []
        
        # 1. 分辨率检查（DPI）
        dpi = img.info.get('dpi', (72, 72))
        if dpi[0] < 300 or dpi[1] < 300:
            issues.append(f"DPI不足: {dpi[0]}×{dpi[1]} (要求: ≥300)")
        
        # 2. 尺寸检查（像素）
        min_pixels = 2000 * 1500  # 最小3百万像素
        pixels = width * height
        if pixels < min_pixels:
            issues.append(f"尺寸过小: {width}×{height} = {pixels:,}像素")
        
        # 3. 格式检查
        if img.format != 'PNG':
            issues.append(f"格式: {img.format} (推荐: PNG)")
        
        # 4. 模式检查
        if img.mode not in ['RGB', 'RGBA']:
            issues.append(f"颜色模式: {img.mode} (推荐: RGB/RGBA)")
        
        # 文件大小
        file_size = os.path.getsize(image_path) / 1024  # KB
        
        return {
            'size': (width, height),
            'dpi': dpi,
            'format': img.format,
            'mode': img.mode,
            'file_size_kb': file_size,
            'issues': issues,
            'pixels': pixels
        }
        
    except Exception as e:
        return {'error': str(e)}

def main():
    figures = [
        ("Figure 1", "/root/.openclaw/workspace/Figure1_Study_Flowchart.png"),
        ("Figure 2", "/root/.openclaw/workspace/Figure2_Group_Comparison_Boxplot.png"),
        ("Figure 4", "/root/.openclaw/workspace/Figure4_Correlation_Scatter.png"),
        ("Figure 5", "/root/.openclaw/workspace/Figure5_Forest_Plot.png"),
        ("Figure 6", "/root/.openclaw/workspace/Figure6_Subgroup_Analysis.png"),
    ]
    
    print("图表期刊标准检查")
    print("=" * 60)
    
    all_ok = True
    
    for name, path in figures:
        if not os.path.exists(path):
            print(f"{name}: ❌ 文件不存在 - {path}")
            all_ok = False
            continue
            
        print(f"\n{name}: {os.path.basename(path)}")
        print("-" * 40)
        
        result = check_figure(path)
        
        if 'error' in result:
            print(f"  错误: {result['error']}")
            all_ok = False
            continue
            
        # 显示基本信息
        width, height = result['size']
        print(f"  尺寸: {width}×{height}像素")
        print(f"  DPI: {result['dpi'][0]}×{result['dpi'][1]}")
        print(f"  格式: {result['format']}, 模式: {result['mode']}")
        print(f"  文件大小: {result['file_size_kb']:.1f} KB")
        print(f"  总像素: {result['pixels']:,}")
        
        # 检查问题
        if result['issues']:
            print(f"  ⚠️ 发现问题:")
            for issue in result['issues']:
                print(f"    • {issue}")
            all_ok = False
        else:
            print(f"  ✅ 符合期刊标准")
    
    print(f"\n{'='*60}")
    if all_ok:
        print("✅ 所有图表符合期刊标准")
    else:
        print("⚠️ 部分图表需要优化")
        
    # 期刊标准要求总结
    print(f"\n期刊标准要求总结:")
    print(f"• 分辨率: ≥300 DPI")
    print(f"• 尺寸: ≥2000×1500像素 (3百万像素)")
    print(f"• 格式: PNG (推荐)")
    print(f"• 颜色模式: RGB 或 RGBA")
    print(f"• 文件大小: 通常<10MB")

if __name__ == "__main__":
    main()