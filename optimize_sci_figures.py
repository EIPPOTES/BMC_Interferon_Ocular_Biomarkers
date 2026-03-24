#!/usr/bin/env python3
"""
SCI期刊图表优化脚本
针对Journal of Affective Disorders要求优化
作者: OpenClaw 眼科科研助手
日期: 2026-03-22
"""

import os
import sys
from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def check_image_quality(filepath):
    """检查图像质量参数"""
    try:
        img = Image.open(filepath)
        width, height = img.size
        file_size = os.path.getsize(filepath)
        
        # 估算DPI (假设标准期刊尺寸)
        estimated_dpi_width = width / 6.69  # 17cm in inches
        estimated_dpi_height = height / 6.69
        estimated_dpi = min(estimated_dpi_width, estimated_dpi_height)
        
        return {
            'filename': os.path.basename(filepath),
            'width_px': width,
            'height_px': height,
            'estimated_dpi': estimated_dpi,
            'file_size_mb': file_size / (1024*1024),
            'format': img.format,
            'mode': img.mode
        }
    except Exception as e:
        return {'filename': os.path.basename(filepath), 'error': str(e)}

def generate_quality_report(figures_dir):
    """生成图表质量评估报告"""
    report = []
    report.append("="*60)
    report.append("SCI期刊图表质量评估报告")
    report.append("期刊: Journal of Affective Disorders")
    report.append("="*60)
    report.append("")
    
    # 标准requirement
    report.append("【期刊标准要求】")
    report.append("- 线图分辨率: ≥1200 DPI")
    report.append("- 照片/混合图: ≥300 DPI")
    report.append("- 推荐格式: TIFF (LZW压缩)")
    report.append("- 字体大小: ≥8 pt (印刷后)")
    report.append("- 线条粗细: ≥0.5 pt")
    report.append("")
    
    # 检查所有图表
    report.append("【当前图表评估】")
    report.append("-"*60)
    
    figure_files = sorted([f for f in os.listdir(figures_dir) 
                          if f.startswith('Figure') and f.endswith(('.png', '.tiff', '.pdf'))])
    
    issues_found = []
    
    for filename in figure_files:
        filepath = os.path.join(figures_dir, filename)
        quality = check_image_quality(filepath)
        
        if 'error' in quality:
            report.append(f"❌ {filename}: 检查失败 - {quality['error']}")
            continue
        
        report.append(f"\n📊 {quality['filename']}")
        report.append(f"   格式: {quality['format']} | 模式: {quality['mode']}")
        report.append(f"   尺寸: {quality['width_px']} x {quality['height_px']} px")
        report.append(f"   估算DPI: {quality['estimated_dpi']:.0f}")
        report.append(f"   文件大小: {quality['file_size_mb']:.2f} MB")
        
        # 评估
        if quality['estimated_dpi'] < 300:
            report.append(f"   ⚠️  分辨率可能不足 (建议≥300 DPI)")
            issues_found.append(f"{filename}: 分辨率低 ({quality['estimated_dpi']:.0f} DPI)")
        elif quality['estimated_dpi'] < 1200 and 'line' in filename.lower():
            report.append(f"   ⚠️  线图建议≥1200 DPI")
            issues_found.append(f"{filename}: 线图分辨率建议提高")
        else:
            report.append(f"   ✅ 分辨率充足")
    
    report.append("")
    report.append("="*60)
    report.append("【优化建议】")
    
    if issues_found:
        report.append("\n发现的问题:")
        for issue in issues_found:
            report.append(f"  - {issue}")
        report.append("\n建议操作:")
        report.append("1. 使用矢量格式(PDF/EPS)重新导出图表")
        report.append("2. 确保字体嵌入且≥8pt")
        report.append("3. 线条粗细设置为≥0.5pt")
        report.append("4. 使用CMYK颜色模式(印刷)")
    else:
        report.append("✅ 所有图表基本符合要求")
    
    return "\n".join(report)

def create_optimized_figure_template():
    """创建符合SCI标准的图表模板"""
    
    # 设置SCI标准参数
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,  # 8pt for printed figure
        'axes.linewidth': 0.5,
        'lines.linewidth': 1.0,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.format': 'tiff',
        'savefig.bbox': 'tight'
    })
    
    # 创建示例森林图模板
    fig, ax = plt.subplots(figsize=(8.5/2.54, 6/2.54))  # 8.5cm width
    
    # 示例数据
    studies = ['Study 1', 'Study 2', 'Study 3', 'Overall']
    effect_sizes = [0.5, 0.3, 0.7, 0.45]
    ci_lower = [0.2, 0.0, 0.4, 0.25]
    ci_upper = [0.8, 0.6, 1.0, 0.65]
    
    y_pos = np.arange(len(studies))
    
    # 绘制森林图
    for i, (study, es, cil, ciu) in enumerate(zip(studies, effect_sizes, ci_lower, ci_upper)):
        ax.plot([cil, ciu], [i, i], 'k-', linewidth=1.5)
        ax.plot(es, i, 's', markersize=8, color='black' if study == 'Overall' else 'gray')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(studies)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Effect Size (95% CI)', fontsize=8)
    ax.set_title('Forest Plot Template (SCI Standard)', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    
    return fig

if __name__ == "__main__":
    if len(sys.argv) < 2:
        figures_dir = "/mnt/c/Users/CUI/Desktop/final/02_Figures"
    else:
        figures_dir = sys.argv[1]
    
    print("🔬 SCI期刊图表优化工具")
    print("="*60)
    
    if os.path.exists(figures_dir):
        report = generate_quality_report(figures_dir)
        print(report)
        
        # 保存报告
        report_path = os.path.join(figures_dir, "SCI_Quality_Report.txt")
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"\n📄 报告已保存: {report_path}")
    else:
        print(f"❌ 目录不存在: {figures_dir}")
        print("用法: python optimize_sci_figures.py <figures_directory>")
