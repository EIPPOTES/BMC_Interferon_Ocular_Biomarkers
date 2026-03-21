#!/usr/bin/env python3
"""
检测Figure 1中的箭头并识别重叠问题
"""

from PIL import Image
import numpy as np
import sys
import os

def detect_arrows(image_path):
    """检测图像中的箭头元素"""
    try:
        img = Image.open(image_path)
        img_np = np.array(img)
        
        print("图像分析:")
        print(f"尺寸: {img.size}")
        
        # 提取RGBA通道
        if img_np.shape[2] == 4:  # RGBA
            r, g, b, a = img_np[:,:,0], img_np[:,:,1], img_np[:,:,2], img_np[:,:,3]
        else:  # RGB
            r, g, b = img_np[:,:,0], img_np[:,:,1], img_np[:,:,2]
            a = np.ones_like(r) * 255
            
        # 寻找可能的箭头颜色（灰色）
        # 箭头通常是灰色（R≈G≈B）
        gray_mask = (np.abs(r - g) < 20) & (np.abs(r - b) < 20) & (np.abs(g - b) < 20)
        
        # 排除白色背景（高亮度）
        brightness = (r.astype(float) + g + b) / 3
        not_white = brightness < 200
        
        # 组合掩码
        arrow_candidate_mask = gray_mask & not_white
        
        # 计算箭头区域
        arrow_pixels = np.sum(arrow_candidate_mask)
        total_pixels = img_np.shape[0] * img_np.shape[1]
        
        print(f"\n箭头检测:")
        print(f"  疑似箭头像素: {arrow_pixels} ({arrow_pixels/total_pixels*100:.2f}%)")
        
        # 寻找箭头集中的区域
        if arrow_pixels > 0:
            # 计算箭头像素的坐标
            arrow_coords = np.argwhere(arrow_candidate_mask)
            
            # 计算中心
            center_y = np.mean(arrow_coords[:, 0])
            center_x = np.mean(arrow_coords[:, 1])
            
            print(f"  箭头区域中心: ({center_x:.1f}, {center_y:.1f})")
            
            # 检查中心区域是否与文字重叠
            # 文字通常是黑色或深色
            text_mask = (r < 100) & (g < 100) & (b < 100)
            
            # 检查箭头和文字的重叠
            overlap_mask = arrow_candidate_mask & text_mask
            overlap_pixels = np.sum(overlap_mask)
            
            print(f"\n重叠检测:")
            print(f"  箭头与文字重叠像素: {overlap_pixels}")
            
            if overlap_pixels > 100:  # 阈值
                print(f"  ⚠️ 检测到显著重叠 ({overlap_pixels}像素)")
                
                # 找到重叠区域
                overlap_coords = np.argwhere(overlap_mask)
                if len(overlap_coords) > 0:
                    overlap_center_y = np.mean(overlap_coords[:, 0])
                    overlap_center_x = np.mean(overlap_coords[:, 1])
                    print(f"  重叠区域中心: ({overlap_center_x:.1f}, {overlap_center_y:.1f})")
                    
                    # 提供修复建议
                    print(f"\n修复建议:")
                    print(f"  1. 在图像编辑软件中，定位到坐标 ({int(overlap_center_x)}, {int(overlap_center_y)}) 附近")
                    print(f"  2. 检查灰色箭头是否覆盖黑色文字")
                    print(f"  3. 缩短箭头或调整位置以避免重叠")
            else:
                print(f"  ✓ 未检测到显著重叠")
        
        return True
        
    except Exception as e:
        print(f"分析出错: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    image_path = "/root/.openclaw/workspace/Figure1_Study_Flowchart.png"
    
    if not os.path.exists(image_path):
        print(f"文件不存在: {image_path}")
        sys.exit(1)
    
    print(f"分析: {image_path}")
    print("=" * 50)
    
    detect_arrows(image_path)