#!/usr/bin/env python3
"""
修复create_figures.py，添加if __name__ == "__main__":保护
"""

import re

def fix_create_figures():
    with open('create_figures.py', 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 检查是否已有main guard
    if 'if __name__ == "__main__":' in content:
        print("已有main guard，无需修复")
        return False
    
    # 分割内容：import部分和执行部分
    lines = content.split('\n')
    
    # 找到import语句的结束位置
    import_end = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped and not stripped.startswith('#') and not stripped.startswith('import') and not stripped.startswith('from'):
            # 找到第一个非import/非注释行
            import_end = i
            break
    
    # 如果没找到，可能在文件末尾
    if import_end == 0:
        import_end = len(lines)
    
    # 创建新内容
    new_lines = []
    
    # 添加import部分
    new_lines.extend(lines[:import_end])
    
    # 添加main函数定义
    new_lines.append('')
    new_lines.append('def main():')
    
    # 添加执行代码（缩进）
    for line in lines[import_end:]:
        if line.strip():
            new_lines.append(f'    {line}')
        else:
            new_lines.append('')
    
    # 添加main guard
    new_lines.append('')
    new_lines.append('')
    new_lines.append('if __name__ == "__main__":')
    new_lines.append('    main()')
    
    # 写入新文件
    new_content = '\n'.join(new_lines)
    
    with open('create_figures.py', 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    print("成功添加main guard保护")
    return True

if __name__ == "__main__":
    fix_create_figures()