#!/usr/bin/env python3
"""
检查修复后Figure 5的质量
"""

from PIL import Image
import os

def check_figure5():
    fig_path = '/root/.openclaw/workspace/figures/Figure5_Forest_Plot.png'
    
    if not os.path.exists(fig_path):
        print(f"错误: 文件不存在: {fig_path}")
        return False
    
    try:
        img = Image.open(fig_path)
        width, height = img.size
        dpi = img.info.get('dpi', (72, 72))
        
        print(f"检查修复后的Figure 5:")
        print(f"  文件: {fig_path}")
        print(f"  尺寸: {width}×{height}像素")
        print(f"  DPI: {dpi[0]}×{dpi[1]}")
        print(f"  格式: {img.format}")
        print(f"  模式: {img.mode}")
        
        # 检查是否满足期刊标准
        issues = []
        if dpi[0] < 300 or dpi[1] < 300:
            issues.append(f"DPI不足: {dpi[0]}×{dpi[1]} (要求: ≥300)")
        
        min_pixels = 2000 * 1500
        pixels = width * height
        if pixels < min_pixels:
            issues.append(f"尺寸过小: {pixels:,}像素 (要求: ≥{min_pixels:,})")
        
        if issues:
            print(f"  ⚠️ 发现问题:")
            for issue in issues:
                print(f"    • {issue}")
            return False
        else:
            print(f"  ✅ 符合期刊标准")
            return True
            
    except Exception as e:
        print(f"检查出错: {e}")
        return False

def compare_with_original():
    """与原始Figure 5比较"""
    new_path = '/root/.openclaw/workspace/figures/Figure5_Forest_Plot.png'
    old_path = '/root/.openclaw/workspace/Figure5_Forest_Plot.png'
    
    if not os.path.exists(old_path):
        print("原始Figure 5不存在，跳过比较")
        return
    
    try:
        new_img = Image.open(new_path)
        old_img = Image.open(old_path)
        
        print(f"\n与原始Figure 5比较:")
        print(f"  新文件: {new_img.size[0]}×{new_img.size[1]}, DPI: {new_img.info.get('dpi', (0,0))}")
        print(f"  旧文件: {old_img.size[0]}×{old_img.size[1]}, DPI: {old_img.info.get('dpi', (0,0))}")
        
        # 检查文件大小
        new_size = os.path.getsize(new_path) / 1024
        old_size = os.path.getsize(old_path) / 1024
        
        print(f"  新文件大小: {new_size:.1f} KB")
        print(f"  旧文件大小: {old_size:.1f} KB")
        
        if abs(new_img.size[0] - old_img.size[0]) > 10 or abs(new_img.size[1] - old_img.size[1]) > 10:
            print(f"  ⚠️ 尺寸有明显差异")
        else:
            print(f"  ✅ 尺寸基本一致")
            
    except Exception as e:
        print(f"比较出错: {e}")

if __name__ == "__main__":
    print("=" * 60)
    print("Figure 5修复后质量检查")
    print("=" * 60)
    
    check_figure5()
    compare_with_original()