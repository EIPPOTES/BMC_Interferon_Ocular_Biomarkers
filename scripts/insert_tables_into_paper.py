#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *

def main():
    """主函数，包装原有执行代码"""
    """
    在论文中插入更新后的表格数据
    """

    import os
    import re
    from datetime import datetime

    print("=" * 80)
    print("在论文中插入更新后的表格数据")
    print(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)

    # 文件路径
    paper_path = "/root/.openclaw/workspace/revised_paper/manuscript_final_integrated_463eyes_ml.md"
    tables_md_path = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/01_Manuscript/Updated_Tables_Markdown_20260315.md"
    output_path = "/root/.openclaw/workspace/revised_paper/manuscript_with_tables_463eyes_ml.md"

    # 读取论文和表格
    print(f"📄 读取论文文件: {paper_path}")
    with open(paper_path, 'r', encoding='utf-8') as f:
        paper_content = f.read()

    print(f"📋 读取表格数据: {tables_md_path}")
    with open(tables_md_path, 'r', encoding='utf-8') as f:
        tables_content = f.read()

    # 分离三个表格
    tables_parts = tables_content.split('## Table ')
    table1_md = '## Table ' + tables_parts[1].split('\n## Table ')[0] if len(tables_parts) > 1 else ''
    table2_md = '## Table ' + tables_parts[2].split('\n## Table ')[0] if len(tables_parts) > 2 else ''
    table3_md = '## Table ' + tables_parts[3] if len(tables_parts) > 3 else ''

    print(f"表格分离: Table1({len(table1_md)} chars), Table2({len(table2_md)} chars), Table3({len(table3_md)} chars)")

    # 1. 插入Table 1
    print(f"\n🔧 插入Table 1...")
    # 找到"summarized in **Table 1**"的位置
    table1_ref_pattern = r'(Demographic and clinical characteristics are summarized in \*\*Table 1\*\*\.)'
    match = re.search(table1_ref_pattern, paper_content)
    if match:
        table1_ref = match.group(1)
        # 在引用后插入表格
        insertion_point = match.end()
        new_content = paper_content[:insertion_point] + '\n\n' + table1_md + '\n\n' + paper_content[insertion_point:]
        paper_content = new_content
        print("✅ Table 1 插入成功")
    else:
        print("⚠️ 未找到Table 1引用")

    # 2. 插入Table 2
    print(f"🔧 插入Table 2...")
    # 找到Table 2引用的位置 (在"3.2 Macular Structural Changes"部分后)
    # 先找到"3.2 Macular Structural Changes"章节
    section_pattern = r'(## 3\.2 Macular Structural Changes[\s\S]*?)(?=\n## |\n# |$)'
    match = re.search(section_pattern, paper_content)
    if match:
        section_text = match.group(1)
        # 在章节末尾插入Table 2
        section_end = match.end()
        # 检查是否已经有表格
        if '## Table 2.' not in section_text:
            new_content = paper_content[:section_end] + '\n\n' + table2_md + '\n\n' + paper_content[section_end:]
            paper_content = new_content
            print("✅ Table 2 插入成功")
        else:
            print("⚠️ Table 2 已存在")
    else:
        print("⚠️ 未找到3.2章节")

    # 3. 插入Table 3
    print(f"🔧 插入Table 3...")
    # 找到"3.3 Optic Disc Structural Changes"章节
    section_pattern = r'(## 3\.3 Optic Disc Structural Changes[\s\S]*?)(?=\n## |\n# |$)'
    match = re.search(section_pattern, paper_content)
    if match:
        section_text = match.group(1)
        # 在章节开始后插入Table 3
        # 找到章节的第一段结束位置
        first_para_end = section_text.find('\n\n')
        if first_para_end > 0:
            section_start = match.start()
            insertion_point = section_start + first_para_end
            # 检查是否已经有表格
            if '## Table 3.' not in section_text:
                new_content = paper_content[:insertion_point] + '\n\n' + table3_md + '\n\n' + paper_content[insertion_point:]
                paper_content = new_content
                print("✅ Table 3 插入成功")
            else:
                print("⚠️ Table 3 已存在")
        else:
            print("⚠️ 无法找到段落结束")
    else:
        print("⚠️ 未找到3.3章节")

    # 4. 更新图片引用部分（可选）
    print(f"🔧 更新图片引用部分...")
    # 找到表格图片引用的部分
    table_images_section = """| **Table 1** | Baseline Characteristics | `Table1_Baseline_Characteristics_ThreeLine_Final.png` |
    | **Table 2** | Macular Layer Analysis (5 layers) | `Table2_Macular_Layers_ThreeLine.png` |
    | **Table 3** | Optic Disc Parameters | `Table3_Optic_Disc_ThreeLine.png` |"""

    if table_images_section in paper_content:
        # 替换为说明文本
        updated_section = """| **Table 1** | Baseline Characteristics | See section 3.1 for detailed table |
    | **Table 2** | Macular Layer Analysis (5 layers) | See section 3.2 for detailed table |
    | **Table 3** | Optic Disc Parameters | See section 3.3 for detailed table |"""
        paper_content = paper_content.replace(table_images_section, updated_section)
        print("✅ 图片引用部分更新成功")

    # 5. 保存更新后的论文
    print(f"\n💾 保存更新后的论文: {output_path}")
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(paper_content)

    # 6. 生成Word文档
    print(f"🔧 生成Word文档...")
    word_output = output_path.replace('.md', '.docx')
    try:
        import subprocess
        result = subprocess.run(['pandoc', output_path, '-o', word_output], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✅ Word文档: {word_output}")
        else:
            print(f"⚠️ Word文档生成失败: {result.stderr}")
    except Exception as e:
        print(f"⚠️ 无法生成Word文档: {e}")

    # 7. 复制到Windows目录
    print(f"🔧 复制到Windows目录...")
    windows_path = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/01_Manuscript/OCT_MDD_Manuscript_463eyes_ML_with_Tables_20260315.md"
    try:
        with open(windows_path, 'w', encoding='utf-8') as f:
            f.write(paper_content)
        print(f"✅ Windows目录: {windows_path}")
    except Exception as e:
        print(f"⚠️ 复制失败: {e}")

    print(f"\n📊 更新完成:")
    print(f"  论文大小: {len(paper_content):,} 字符")
    print(f"  输出文件: {output_path}")
    print(f"  Windows文件: {windows_path}")

    print(f"\n" + "=" * 80)
    print("🎉 表格数据插入完成!")
    print("论文现在包含实际的Table 1, 2, 3数据")
    print("=" * 80)


if __name__ == "__main__":
    main()