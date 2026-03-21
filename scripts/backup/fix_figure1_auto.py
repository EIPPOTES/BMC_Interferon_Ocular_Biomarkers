#!/usr/bin/env python3
"""
自动修复Figure 1箭头重叠问题
尝试使用图像处理技术自动检测和修复箭头与文字的重叠
"""

from PIL import Image, ImageDraw, ImageFilter
import numpy as np
import os

def analyze_and_fix_figure1():
    """分析并尝试修复Figure 1"""
    
    input_path = '/root/.openclaw/workspace/Figure1_Study_Flowchart.png'
    output_path = '/root/.openclaw/workspace/Figure1_Study_Flowchart_fixed.png'
    
    print("="*60)
    print("Figure 1 箭头重叠自动修复")
    print("="*60)
    
    # 加载图像
    print(f"\n1. 加载图像: {input_path}")
    img = Image.open(input_path)
    img_array = np.array(img)
    
    print(f"   图像尺寸: {img.size}")
    print(f"   图像模式: {img.mode}")
    
    # 分析重叠区域
    print(f"\n2. 分析重叠区域...")
    
    # 转换为RGBA（如果不是）
    if img.mode != 'RGBA':
        img = img.convert('RGBA')
    
    # 获取图像数据
    width, height = img.size
    pixels = img.load()
    
    # 定义问题区域（根据之前的分析，中心在1480, 1006）
    center_x, center_y = 1480, 1006
    region_size = 200
    
    # 提取问题区域
    left = max(0, center_x - region_size)
    upper = max(0, center_y - region_size)
    right = min(width, center_x + region_size)
    lower = min(height, center_y + region_size)
    
    print(f"   问题区域: ({left}, {upper}) - ({right}, {lower})")
    
    # 创建修复后的图像副本
    fixed_img = img.copy()
    draw = ImageDraw.Draw(fixed_img)
    
    # 方法1: 尝试检测并淡化箭头
    print(f"\n3. 尝试自动修复...")
    
    # 在问题区域，尝试找到灰色箭头并调整其透明度
    fixed_pixels = 0
    
    for y in range(upper, lower):
        for x in range(left, right):
            r, g, b, a = pixels[x, y]
            
            # 检测灰色箭头（R≈G≈B，且不是白色或黑色）
            if abs(int(r) - int(g)) < 30 and abs(int(r) - int(b)) < 30 and abs(int(g) - int(b)) < 30:
                brightness = (int(r) + int(g) + int(b)) / 3
                # 如果是中等亮度的灰色（可能是箭头）
                if 80 < brightness < 180 and a > 100:
                    # 降低箭头的不透明度（使其更透明）
                    new_a = int(a * 0.3)  # 降低到30%不透明度
                    fixed_img.putpixel((x, y), (r, g, b, new_a))
                    fixed_pixels += 1
    
    print(f"   修复像素数: {fixed_pixels}")
    
    if fixed_pixels > 0:
        print(f"   ✅ 已调整箭头透明度")
    else:
        print(f"   ⚠️ 未检测到可修复的箭头像素")
    
    # 方法2: 在文本区域添加白色背景层（如果方法1效果不佳）
    print(f"\n4. 添加文本背景增强...")
    
    # 创建白色背景覆盖层
    overlay = Image.new('RGBA', img.size, (255, 255, 255, 0))
    overlay_draw = ImageDraw.Draw(overlay)
    
    # 在关键文本区域添加半透明白色背景
    text_regions = [
        (center_x - 150, center_y - 50, center_x + 150, center_y + 50),
    ]
    
    for region in text_regions:
        overlay_draw.rectangle(region, fill=(255, 255, 255, 180))
    
    # 合并图层
    fixed_img = Image.alpha_composite(fixed_img, overlay)
    
    # 保存修复后的图像
    print(f"\n5. 保存修复后的图像...")
    fixed_img.save(output_path, 'PNG', dpi=(300, 300))
    print(f"   输出: {output_path}")
    
    # 验证修复
    print(f"\n6. 验证修复效果...")
    original_size = os.path.getsize(input_path) / 1024
    fixed_size = os.path.getsize(output_path) / 1024
    
    print(f"   原文件: {original_size:.1f} KB")
    print(f"   修复后: {fixed_size:.1f} KB")
    
    # 重新分析修复后的图像
    fixed_array = np.array(fixed_img)
    if len(fixed_array.shape) == 3 and fixed_array.shape[2] == 4:
        # 检查透明度变化
        alpha_channel = fixed_array[:, :, 3]
        transparent_pixels = np.sum(alpha_channel < 200)
        print(f"   透明/半透明像素: {transparent_pixels}")
    
    print(f"\n" + "="*60)
    print("自动修复完成")
    print("="*60)
    print(f"\n注意: 自动修复可能不够完美。建议:")
    print(f"1. 检查修复后的图像: {output_path}")
    print(f"2. 如效果不佳，使用GIMP进行手动精细修复")
    print(f"3. GIMP修复步骤:")
    print(f"   - 打开: gimp {input_path}")
    print(f"   - 定位到中心区域 (1480, 1006)")
    print(f"   - 使用橡皮擦工具或图层调整修复箭头")
    print(f"   - 保存为: {output_path}")
    
    return True

if __name__ == "__main__":
    try:
        analyze_and_fix_figure1()
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()