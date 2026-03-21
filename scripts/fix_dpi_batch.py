#!/usr/bin/env python3
"""
批量修复所有Figures的DPI至300
使用PIL/Pillow修改图像元数据
"""

from PIL import Image
import os
import glob

def fix_dpi(image_path, output_path=None):
    """修复单个图像的DPI"""
    try:
        img = Image.open(image_path)
        
        # 保存时设置DPI为300
        if output_path is None:
            output_path = image_path
        
        # 保留原始图像的所有元数据
        img.save(output_path, 'PNG', dpi=(300, 300))
        
        return True
    except Exception as e:
        print(f"  错误: {e}")
        return False

def batch_fix_dpi():
    """批量修复所有Figures"""
    
    print("="*60)
    print("批量修复所有Figures的DPI至300")
    print("="*60)
    
    # 定义需要修复的figures
    figures = [
        ('Figure1', '/root/.openclaw/workspace/Figure1_Study_Flowchart_fixed.png'),
        ('Figure2', '/root/.openclaw/workspace/Figure2_Group_Comparison_Boxplot.png'),
        ('Figure3', '/root/.openclaw/workspace/Figure3_ROC_Curves.png'),
        ('Figure4', '/root/.openclaw/workspace/Figure4_Correlation_Scatter.png'),
        ('Figure5', '/root/.openclaw/workspace/Figure5_Forest_Plot.png'),
        ('Figure6', '/root/.openclaw/workspace/Figure6_Subgroup_Analysis.png'),
    ]
    
    # 同时修复figures目录中的版本
    figures_dir = '/root/.openclaw/workspace/figures'
    if os.path.exists(figures_dir):
        figures.extend([
            ('Figure1 (alt)', f'{figures_dir}/Figure1_Study_Flowchart.png'),
            ('Figure2 (alt)', f'{figures_dir}/Figure2_Group_Comparison.png'),
            ('Figure3 (alt)', f'{figures_dir}/Figure3_ROC_Curves.png'),
            ('Figure4 (alt)', f'{figures_dir}/Figure4_Correlation_Scatter.png'),
            ('Figure5 (alt)', f'{figures_dir}/Figure5_Forest_Plot.png'),
            ('Figure6 (alt)', f'{figures_dir}/Figure6_Subgroup_Analysis.png'),
        ])
    
    results = []
    
    for name, path in figures:
        if not os.path.exists(path):
            print(f"\n⚠️ {name}: 文件不存在 - {path}")
            continue
        
        print(f"\n处理: {name}")
        print(f"  文件: {path}")
        
        # 获取原始信息
        img = Image.open(path)
        original_dpi = img.info.get('dpi', (0, 0))
        original_size = os.path.getsize(path) / 1024
        
        print(f"  原始DPI: {original_dpi[0]:.1f}×{original_dpi[1]:.1f}")
        print(f"  原始大小: {original_size:.1f} KB")
        
        # 修复DPI
        if fix_dpi(path):
            # 验证修复结果
            img_fixed = Image.open(path)
            fixed_dpi = img_fixed.info.get('dpi', (0, 0))
            fixed_size = os.path.getsize(path) / 1024
            
            print(f"  修复后DPI: {fixed_dpi[0]:.1f}×{fixed_dpi[1]:.1f}")
            print(f"  修复后大小: {fixed_size:.1f} KB")
            
            if fixed_dpi[0] >= 300 and fixed_dpi[1] >= 300:
                print(f"  ✅ DPI修复成功")
                results.append((name, True))
            else:
                print(f"  ⚠️ DPI可能未正确设置")
                results.append((name, False))
        else:
            print(f"  ❌ 修复失败")
            results.append((name, False))
    
    # 总结
    print("\n" + "="*60)
    print("批量修复完成")
    print("="*60)
    
    success_count = sum(1 for _, success in results if success)
    total_count = len(results)
    
    print(f"\n修复统计:")
    print(f"  成功: {success_count}/{total_count}")
    print(f"  失败: {total_count - success_count}/{total_count}")
    
    if success_count == total_count:
        print(f"\n🎉 所有Figures的DPI已成功修复至300！")
    else:
        print(f"\n⚠️ 部分Figures修复失败，请检查错误信息")
    
    return success_count == total_count

if __name__ == "__main__":
    batch_fix_dpi()