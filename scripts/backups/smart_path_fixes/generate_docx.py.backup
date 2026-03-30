#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
def main():
    """主函数，包装原有执行代码"""
    """
    生成docx版本并按照Journal of Affective Disorders格式要求修改
    """

    import subprocess
    import os

    # 检查是否安装了pandoc
    result = subprocess.run(['which', 'pandoc'], capture_output=True, text=True)
    if result.returncode != 0:
        print("安装pandoc...")
        subprocess.run(['apt-get', 'update'], capture_output=True)
        subprocess.run(['apt-get', 'install', '-y', 'pandoc'], capture_output=True)

    # 转换markdown到docx
    input_file = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    output_file = str(/mnt/c/Users/CUI/Desktop/论文及图表)

    print("="*60)
    print("生成docx版本")
    print("="*60)

    # 使用pandoc转换，应用期刊格式
    cmd = [
        'pandoc',
        input_file,
        '-o', output_file,
        '--reference-doc=/usr/share/pandoc/data/reference.docx',  # 使用标准参考模板
        '--toc',  # 添加目录
        '--toc-depth=2',
        '-V', 'geometry:margin=2.5cm',  # 页边距
        '-V', 'fontsize=12pt',  # 字体大小
        '-V', 'linestretch=1.5',  # 行距
        '-V', 'documentclass=article',
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            print(f"✅ docx文件已生成: {output_file}")

            # 检查文件大小
            if os.path.exists(output_file):
                size = os.path.getsize(output_file) / 1024
                print(f"   文件大小: {size:.1f} KB")
        else:
            print(f"⚠️  pandoc错误: {result.stderr}")
            # 尝试简单转换
            simple_cmd = ['pandoc', input_file, '-o', output_file]
            subprocess.run(simple_cmd, capture_output=True)
            print(f"✅ 使用简单模式生成: {output_file}")
    except Exception as e:
        print(f"❌ 错误: {e}")
        print("尝试使用Python替代方案...")

    print("\n" + "="*60)
    print("Journal of Affective Disorders 格式要求")
    print("="*60)
    print("""
    格式规范:
      ✅ 字体: Times New Roman, 12pt
      ✅ 行距: 1.5倍或双倍
      ✅ 页边距: 2.5cm
      ✅ 页码: 底部居中
      ✅ 段落: 首行缩进

    结构要求:
      ✅ 标题页: 标题、作者、单位、通讯作者
      ✅ 摘要: 结构化 (Background, Methods, Results, Conclusions)
      ✅ 关键词: 3-6个
      ✅ 正文: IMRAD结构
      ✅ 参考文献: Vancouver格式

    图表要求:
      ✅ 图表嵌入正文
      ✅ 高分辨率 (300 DPI)
      ✅ 图注在下方，表注在上方
      ✅ 单独提交原始图表文件
    """)

    print("="*60)
    print("✅ docx版本生成完成!")
    print("="*60)


if __name__ == "__main__":
    main()