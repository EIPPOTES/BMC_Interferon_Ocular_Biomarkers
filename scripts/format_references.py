#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *

def main():
    """主函数，包装原有执行代码"""
    """
    快速参考文献格式化脚本，尝试转换为Vancouver风格。
    """

    import re
    import sys

    def format_reference_vancouver(line):
        """尝试将单条参考文献格式化为Vancouver风格"""
        # 跳过非参考文献行
        if not re.match(r'^\d+\.\s+', line):
            return line

        # 移除编号
        ref = re.sub(r'^\d+\.\s+', '', line).strip()

        # 常见模式匹配
        # 模式1: 作者. 标题. 期刊. 年份;卷(期):页码.
        # 模式2: 作者. 标题. 出版地: 出版社; 年份.
        # 模式3: 作者. 标题. 期刊. 年份;卷(Suppl X):页码.

        # 这里我们只做一些基本的清理和标准化
        # 1. 确保卷期格式一致: 卷(期)
        ref = re.sub(r'(\d+)\((\d+)\)', r'\1(\2)', ref)  # 确保没有空格

        # 2. 确保页码格式一致: 页码-页码
        ref = re.sub(r'(\d+)-(\d+)', r'\1-\2', ref)

        # 3. 确保年份后使用分号
        ref = re.sub(r'\. (\d{4})', r'. \1;', ref)

        # 4. 移除多余空格
        ref = re.sub(r'\s+', ' ', ref)

        # 返回格式化后的行（暂时保持原样，稍后添加更多逻辑）
        return ref

    def main():
        if len(sys.argv) < 2:
            print("使用方法: python format_references.py <输入文件>")
            return

        input_file = sys.argv[1]

        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        in_references = False
        formatted_lines = []

        for line in lines:
            stripped = line.strip()

            # 检测参考文献部分开始
            if stripped.startswith('# References'):
                in_references = True
                formatted_lines.append(line)
                continue

            # 检测参考文献部分结束
            if in_references and stripped.startswith('---') and 'Word count' in line:
                in_references = False
                formatted_lines.append(line)
                continue

            if in_references:
                # 格式化参考文献行
                formatted = format_reference_vancouver(line.rstrip('\n'))
                formatted_lines.append(formatted + '\n')
            else:
                formatted_lines.append(line)

        # 写回文件
        with open(input_file, 'w', encoding='utf-8') as f:
            f.writelines(formatted_lines)

        print(f"已处理文件: {input_file}")

    if __name__ == '__main__':
        main()


if __name__ == "__main__":
    main()