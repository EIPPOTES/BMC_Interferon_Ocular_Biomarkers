#!/usr/bin/env python3
"""
脚本质量检查工具
检查Python脚本的质量问题：
1. 是否缺少if __name__ == "__main__":保护
2. 是否有硬编码路径
3. 其他常见问题
"""

import os
import sys
import re
from pathlib import Path

# 设置工作目录
current_dir = Path(__file__).parent
scripts_to_check = list(current_dir.glob("*.py"))

def check_main_guard(script_path):
    """检查脚本是否有if __name__ == "__main__":保护"""
    with open(script_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    has_main_guard = 'if __name__ == "__main__":' in content
    
    # 检查是否有直接执行的代码（不在函数或类中）
    lines = content.split('\n')
    has_top_level_code = False
    
    # 简单检查：是否有不在函数/类中的代码行（排除import和注释）
    in_function = False
    in_class = False
    for line in lines:
        line = line.strip()
        
        # 跳过空行、注释、import
        if not line or line.startswith('#') or line.startswith('import') or line.startswith('from'):
            continue
        
        # 检查是否进入/退出函数或类
        if line.startswith('def ') or line.startswith('class '):
            if line.endswith(':'):
                if line.startswith('def '):
                    in_function = True
                elif line.startswith('class '):
                    in_class = True
        elif line == 'if __name__ == "__main__":':
            break  # 遇到main guard，后面的代码不算顶层代码
        
        # 如果不在函数或类中，且不是import/注释，且不是if __name__...
        elif not in_function and not in_class and '=' in line or '(' in line:
            has_top_level_code = True
            break
    
    return has_main_guard, has_top_level_code

def check_hardcoded_paths(script_path):
    """检查是否有硬编码路径"""
    with open(script_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # 常见的硬编码路径模式
    patterns = [
        r'"/root/\.openclaw/workspace/',
        r"'/root/\.openclaw/workspace/",
        r'C:\\',  # Windows绝对路径
        r'D:\\',
        r'/mnt/c/Users/',  # WSL中的Windows路径
        r'\./data\.xlsx',  # 相对路径但应该使用配置
        r'data\.xlsx',     # 直接引用数据文件
    ]
    
    hardcoded_paths = []
    for pattern in patterns:
        matches = re.findall(pattern, content)
        if matches:
            hardcoded_paths.extend(matches)
    
    return list(set(hardcoded_paths))

def main():
    print("=== Python脚本质量检查 ===\n")
    
    issues_found = False
    scripts_without_main_guard = []
    scripts_with_top_level_code = []
    scripts_with_hardcoded_paths = []
    
    for script_path in scripts_to_check:
        if script_path.name == __file__:  # 跳过自己
            continue
            
        print(f"检查: {script_path.name}")
        
        # 检查main guard
        has_main_guard, has_top_level_code = check_main_guard(script_path)
        
        if not has_main_guard:
            scripts_without_main_guard.append(script_path.name)
            print("  ⚠️  缺少 if __name__ == '__main__': 保护")
        
        if has_top_level_code:
            scripts_with_top_level_code.append(script_path.name)
            print("  ⚠️  有顶层执行代码（可能导入时自动执行）")
        
        # 检查硬编码路径
        hardcoded_paths = check_hardcoded_paths(script_path)
        if hardcoded_paths:
            scripts_with_hardcoded_paths.append(script_path.name)
            print(f"  ⚠️  硬编码路径: {', '.join(hardcoded_paths[:3])}")
            if len(hardcoded_paths) > 3:
                print(f"    ... 还有 {len(hardcoded_paths)-3} 个")
        
        if has_main_guard and not has_top_level_code and not hardcoded_paths:
            print("  ✅ 质量良好")
        
        print()
    
    # 输出汇总报告
    print("=== 质量检查汇总 ===")
    print(f"检查脚本总数: {len(scripts_to_check)-1}")  # 排除自己
    
    if scripts_without_main_guard:
        print(f"\n❌ 缺少if __name__ == '__main__':保护的脚本 ({len(scripts_without_main_guard)}个):")
        for script in sorted(scripts_without_main_guard)[:10]:
            print(f"  - {script}")
        if len(scripts_without_main_guard) > 10:
            print(f"  ... 还有 {len(scripts_without_main_guard)-10} 个")
        issues_found = True
    
    if scripts_with_top_level_code:
        print(f"\n⚠️  有顶层执行代码的脚本 ({len(scripts_with_top_level_code)}个):")
        for script in sorted(scripts_with_top_level_code)[:10]:
            print(f"  - {script}")
        if len(scripts_with_top_level_code) > 10:
            print(f"  ... 还有 {len(scripts_with_top_level_code)-10} 个")
        issues_found = True
    
    if scripts_with_hardcoded_paths:
        print(f"\n⚠️  有硬编码路径的脚本 ({len(scripts_with_hardcoded_paths)}个):")
        for script in sorted(scripts_with_hardcoded_paths)[:10]:
            print(f"  - {script}")
        if len(scripts_with_hardcoded_paths) > 10:
            print(f"  ... 还有 {len(scripts_with_hardcoded_paths)-10} 个")
        issues_found = True
    
    if not issues_found:
        print("\n✅ 所有脚本质量良好！")
    else:
        print(f"\n📊 问题统计:")
        print(f"  - 缺少main guard: {len(scripts_without_main_guard)}")
        print(f"  - 有顶层执行代码: {len(scripts_with_top_level_code)}")
        print(f"  - 有硬编码路径: {len(scripts_with_hardcoded_paths)}")
    
    return 0 if not issues_found else 1

if __name__ == "__main__":
    sys.exit(main())