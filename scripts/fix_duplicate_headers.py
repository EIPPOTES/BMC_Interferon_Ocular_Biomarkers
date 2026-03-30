#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *

def main():
    """主函数，包装原有执行代码"""
    """
    修复重复的章节标题
    """

    paper_path = '/root/.openclaw/workspace/revised_paper/manuscript_final_integrated_updated.md'

    with open(paper_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # 查找重复的3.10标题
    lines = content.split('\n')
    fixed_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        # 检查是否出现连续的相同章节标题
        if line.startswith('## 3.10') and i+1 < len(lines) and lines[i+1].startswith('## 3.10'):
            # 跳过第一个重复标题
            fixed_lines.append(line)
            # 跳过下一个重复标题
            i += 2
            # 添加后续内容
            while i < len(lines) and not lines[i].startswith('## 3.'):
                fixed_lines.append(lines[i])
                i += 1
            continue
        fixed_lines.append(line)
        i += 1

    fixed_content = '\n'.join(fixed_lines)

    # 保存修复后的文件
    fixed_path = '/root/.openclaw/workspace/revised_paper/manuscript_final_integrated_updated_fixed.md'
    with open(fixed_path, 'w', encoding='utf-8') as f:
        f.write(fixed_content)

    print(f"修复后的文件已保存: {fixed_path}")

    # 检查修复结果
    import re
    headers = re.findall(r'## 3\.\d+\s+.+', fixed_content)
    print("\n章节标题检查:")
    for h in headers[:15]:
        print(f"  {h}")


if __name__ == "__main__":
    main()