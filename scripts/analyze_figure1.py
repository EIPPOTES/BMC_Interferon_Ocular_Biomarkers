#!/usr/bin/env python3
"""
分析Figure 1图像，检测箭头与文字重叠问题
"""

from PIL import Image
import numpy as np
import sys
import os

def analyze_image(image_path):
    """分析图像，检测潜在问题"""
    try:
        # 打开图像
        img = Image.open(image_path)
        print(f"图像信息:")
        print(f"  尺寸: {img.size} (宽×高)")
        print(f"  模式: {img.mode}")
        print(f"  格式: {img.format}")
        
        # 转换为numpy数组
        img_np = np.array(img)
        
        # 检查图像质量
        print(f"\n图像质量检查:")
        
        # 检查图像基本信息
        height, width = img_np.shape[0], img_np.shape[1]
        print(f"  NumPy数组形状: {img_np.shape}")
        
        # 转换为灰度图
        if len(img_np.shape) == 3:
            # RGB图像
            gray_img = np.mean(img_np[:, :, :3], axis=2).astype(np.uint8)
        else:
            # 已经是灰度图
            gray_img = img_np
            
        # 检查图像中心区域（通常重叠问题出现在中心）
        center_y, center_x = height // 2, width // 2
        check_size = min(100, height//4, width//4)  # 检查区域大小
        
        print(f"\n中心区域检查 (x={center_x-check_size}:{center_x+check_size}, y={center_y-check_size}:{center_y+check_size}):")
        
        # 检查中心区域的颜色变化
        y_start = max(0, center_y - check_size)
        y_end = min(height, center_y + check_size)
        x_start = max(0, center_x - check_size)
        x_end = min(width, center_x + check_size)
        
        center_region = gray_img[y_start:y_end, x_start:x_end]
        if center_region.size > 0:
            print(f"  中心区域平均亮度: {np.mean(center_region):.1f}")
            print(f"  中心区域亮度标准差: {np.std(center_region):.1f}")
            
            # 高标准差可能表示有重叠
            if np.std(center_region) > 40:
                print(f"  ⚠️ 中心区域亮度变化大，可能存在重叠")
            else:
                print(f"  ✓ 中心区域亮度变化正常")
        
        # 简单重叠检测：检查是否有低亮度区域覆盖在高亮度文本上
        # 这是简化版本，实际需要更复杂的检测
        
        print(f"\n建议:")
        print(f"  1. 人工检查图像中心区域是否有箭头覆盖文字")
        print(f"  2. 如果需要修复，建议使用图像编辑软件（如GIMP）")
        print(f"  3. 可能的修复方法:")
        print(f"     - 缩短箭头长度")
        print(f"     - 调整箭头位置")
        print(f"     - 重新排列文本框")
        
        return True
        
    except Exception as e:
        print(f"分析图像时出错: {e}")
        return False

if __name__ == "__main__":
    image_path = "/root/.openclaw/workspace/Figure1_Study_Flowchart.png"
    
    if not os.path.exists(image_path):
        print(f"错误: 图像文件不存在: {image_path}")
        sys.exit(1)
    
    print(f"分析图像: {image_path}")
    print("=" * 50)
    
    analyze_image(image_path)