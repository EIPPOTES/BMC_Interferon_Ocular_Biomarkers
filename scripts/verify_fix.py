#!/usr/bin/env python3
# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
验证Figure 1修复效果
"""

from PIL import Image
import numpy as np
import os

def check_fix_result():
    original_path = '/root/.openclaw/workspace/Figure1_Study_Flowchart.png'
    fixed_path = '/root/.openclaw/workspace/Figure1_Study_Flowchart_fixed.png'
    
    print("="*60)
    print("Figure 1 修复效果验证")
    print("="*60)
    
    # 检查文件
    if not os.path.exists(fixed_path):
        print(f"❌ 修复后的文件不存在: {fixed_path}")
        return
    
    print(f"\n1. 文件信息:")
    print(f"   原始文件: {os.path.getsize(original_path)/1024:.1f} KB")
    print(f"   修复文件: {os.path.getsize(fixed_path)/1024:.1f} KB")
    
    # 加载修复后的图像
    img = Image.open(fixed_path)
    img_array = np.array(img)
    
    print(f"\n2. 图像属性:")
    print(f"   尺寸: {img.size}")
    print(f"   模式: {img.mode}")
    
    # 分析透明度变化
    if img.mode == 'RGBA' and len(img_array.shape) == 3 and img_array.shape[2] == 4:
        alpha = img_array[:, :, 3]
        
        print(f"\n3. 透明度分析:")
        print(f"   完全透明像素 (α=0): {np.sum(alpha == 0)}")
        print(f"   半透明像素 (α<128): {np.sum(alpha < 128)}")
        print(f"   低透明度像素 (α<200): {np.sum(alpha < 200)}")
        print(f"   不透明像素 (α=255): {np.sum(alpha == 255)}")
        
        # 检查问题区域的透明度
        center_x, center_y = 1480, 1006
        region_size = 200
        
        y_start = max(0, center_y - region_size)
        y_end = min(img.size[1], center_y + region_size)
        x_start = max(0, center_x - region_size)
        x_end = min(img.size[0], center_x + region_size)
        
        region_alpha = alpha[y_start:y_end, x_start:x_end]
        
        print(f"\n4. 问题区域透明度 (中心1480,1006, 范围200像素):")
        print(f"   区域大小: {region_alpha.shape}")
        print(f"   平均透明度: {np.mean(region_alpha):.1f}")
        print(f"   最小透明度: {np.min(region_alpha)}")
        print(f"   半透明像素 (α<200): {np.sum(region_alpha < 200)}")
        
        if np.sum(region_alpha < 200) > 1000:
            print(f"   ✅ 问题区域有大量半透明像素，箭头透明度已调整")
        else:
            print(f"   ⚠️ 问题区域透明度变化不明显")
    
    # 检查DPI
    dpi = img.info.get('dpi', (0, 0))
    print(f"\n5. DPI信息:")
    print(f"   DPI: {dpi[0]}×{dpi[1]}")
    if dpi[0] >= 300 and dpi[1] >= 300:
        print(f"   ✅ 符合期刊标准 (≥300)")
    else:
        print(f"   ⚠️ 不符合期刊标准 (<300)")
    
    print(f"\n" + "="*60)
    print("验证完成")
    print("="*60)
    
    print(f"\n修复文件位置:")
    print(f"  {fixed_path}")
    print(f"\n建议:")
    print(f"1. 查看修复后的图像确认效果")
    print(f"2. 如自动修复效果不理想，使用GIMP手动修复:")
    print(f"   gimp {original_path}")
    print(f"3. 手动修复步骤:")
    print(f"   - 定位到中心区域 (1480, 1006)")
    print(f"   - 使用图层蒙版或橡皮擦工具调整箭头")
    print(f"   - 保存修复结果")

if __name__ == "__main__":
    check_fix_result()