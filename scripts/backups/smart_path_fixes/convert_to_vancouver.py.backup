#!/usr/bin/env python3

def main():
    """主函数，包装原有执行代码"""
    """
    快速将参考文献转换为Vancouver风格。
    Vancouver格式：作者. 题目. 期刊. 年份;卷(期):页码.
    """

    import re
    import sys

    def parse_reference(ref_text):
        """解析参考文献条目，提取结构信息"""
        # 移除编号
        ref = re.sub(r'^\d+\.\s+', '', ref_text).strip()

        # 初始化结果字典
        result = {
            'raw': ref_text,
            'original': ref,
            'authors': '',
            'title': '',
            'journal': '',
            'year': '',
            'volume': '',
            'issue': '',
            'pages': '',
            'type': 'journal'  # 默认期刊文章
        }

        # 尝试检测不同类型
        # 1. 期刊文章：通常有卷(期):页码格式
        journal_match = re.search(r'\.\s+([^.]+)\.\s+(\d{4});(\d+)\((\d+)\):([\d\-]+)\.$', ref)
        if journal_match:
            result['journal'] = journal_match.group(1)
            result['year'] = journal_match.group(2)
            result['volume'] = journal_match.group(3)
            result['issue'] = journal_match.group(4)
            result['pages'] = journal_match.group(5)
            result['type'] = 'journal'
            # 提取作者和标题（假设在期刊名称之前）
            parts = ref.split('.' + result['journal'])[0].split('.')
            if len(parts) >= 2:
                result['authors'] = parts[0].strip()
                result['title'] = '.'.join(parts[1:]).strip() + '.'
            return result

        # 2. 带有分号的期刊文章变体
        journal_match2 = re.search(r'\.\s+([^.]+)\.\s+(\d{4});(\d+):([\d\-]+)\.$', ref)
        if journal_match2:
            result['journal'] = journal_match2.group(1)
            result['year'] = journal_match2.group(2)
            result['volume'] = journal_match2.group(3)
            result['pages'] = journal_match2.group(4)
            result['type'] = 'journal'
            parts = ref.split('.' + result['journal'])[0].split('.')
            if len(parts) >= 2:
                result['authors'] = parts[0].strip()
                result['title'] = '.'.join(parts[1:]).strip() + '.'
            return result

        # 3. 书籍或报告：通常有出版地和出版社
        book_match = re.search(r'\.\s+([^:]+):\s+([^;]+);\s+(\d{4})\.$', ref)
        if book_match:
            result['publisher'] = book_match.group(2)
            result['year'] = book_match.group(3)
            result['type'] = 'book'
            parts = ref.split(':' + result['publisher'])[0].split('.')
            if len(parts) >= 2:
                result['authors'] = parts[0].strip()
                result['title'] = '.'.join(parts[1:]).strip() + '.'
            return result

        # 4. 会议论文或其他格式
        # 暂时保持原样
        return result

    def format_vancouver(ref_info):
        """格式化为Vancouver风格"""
        if ref_info['type'] == 'journal' and ref_info['volume'] and ref_info['pages']:
            # 完整期刊文章
            issue_part = f"({ref_info['issue']})" if ref_info['issue'] else ''
            return f"{ref_info['authors']}. {ref_info['title']} {ref_info['journal']}. {ref_info['year']};{ref_info['volume']}{issue_part}:{ref_info['pages']}."
        elif ref_info['type'] == 'book':
            # 书籍
            return f"{ref_info['authors']}. {ref_info['title']} {ref_info.get('publisher', '')}; {ref_info['year']}."
        else:
            # 无法解析，返回原始（稍作清理）
            cleaned = re.sub(r'\s+', ' ', ref_info['original'])
            cleaned = cleaned.replace(';;', ';')
            return cleaned

    def main():
        if len(sys.argv) < 2:
            print("使用方法: python convert_to_vancouver.py <输入文件>")
            return

        input_file = sys.argv[1]

        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        in_references = False
        formatted_lines = []
        ref_count = 0

        for line in lines:
            stripped = line.strip()

            # 检测参考文献部分开始
            if stripped.startswith('# References'):
                in_references = True
                formatted_lines.append(line)
                # 移除格式说明
                if '*Note:' in line:
                    formatted_lines.append('\n')
                continue

            # 跳过格式说明行
            if in_references and stripped.startswith('*Note:'):
                continue

            # 检测参考文献部分结束
            if in_references and stripped.startswith('---') and 'Word count' in line:
                in_references = False
                formatted_lines.append(line)
                continue

            if in_references and stripped and not stripped.startswith('*'):
                # 是参考文献行
                ref_count += 1
                ref_info = parse_reference(stripped)
                formatted = format_vancouver(ref_info)
                formatted_lines.append(f"{ref_count}. {formatted}\n")
            else:
                formatted_lines.append(line)

        # 写回文件
        with open(input_file, 'w', encoding='utf-8') as f:
            f.writelines(formatted_lines)

        print(f"已处理 {ref_count} 条参考文献，转换为Vancouver风格")
        print(f"文件已保存: {input_file}")

    if __name__ == '__main__':
        main()


if __name__ == "__main__":
    main()