#!/usr/bin/env python3
"""
安全修复脚本：仅为Python脚本添加if __name__ == "__main__":保护
不修改其他内容，确保脚本功能不受影响
"""

import os
import sys
import shutil
from pathlib import Path

def backup_script(script_path):
    """备份脚本文件"""
    backup_dir = Path(__file__).parent / "backups" / "main_guard_fixes"
    backup_dir.mkdir(parents=True, exist_ok=True)
    
    backup_path = backup_dir / f"{script_path.name}.backup"
    shutil.copy2(script_path, backup_path)
    return backup_path

def needs_main_guard(content):
    """检查脚本是否需要添加main guard"""
    # 如果已经有main guard，不需要修复
    if 'if __name__ == "__main__":' in content:
        return False
    
    # 检查是否有顶层执行代码（不在函数中）
    lines = content.split('\n')
    
    # 简单启发式：检查是否有直接调用的函数或代码
    # 跳过import、注释、空行、函数/类定义
    in_function = False
    in_class = False
    
    for line in lines:
        stripped = line.strip()
        
        # 跳过空行、注释、import
        if not stripped or stripped.startswith('#') or stripped.startswith('import') or stripped.startswith('from'):
            continue
        
        # 检查函数或类定义
        if stripped.startswith('def ') or stripped.startswith('class '):
            if stripped.endswith(':'):
                if stripped.startswith('def '):
                    in_function = True
                elif stripped.startswith('class '):
                    in_class = True
            continue
        
        # 检查是否退出函数/类（通过缩进判断）
        if not line.startswith(' ') and not line.startswith('\t'):
            in_function = False
            in_class = False
        
        # 如果不在函数/类中，且看起来是执行代码
        if not in_function and not in_class:
            # 检查是否有函数调用、赋值、表达式
            if '(' in stripped and ')' in stripped:  # 函数调用
                return True
            if '=' in stripped and not '==' in stripped and not '!=' in stripped:  # 赋值
                return True
            if stripped.startswith('print(') or stripped.startswith('print '):  # print语句
                return True
    
    return False

def add_main_guard_safely(content):
    """安全地添加main guard"""
    lines = content.split('\n')
    
    # 找到import语句的结束位置
    import_end = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped and not stripped.startswith('#') and not stripped.startswith('import') and not stripped.startswith('from'):
            import_end = i
            break
    
    if import_end == 0:
        # 没有import语句，找到第一个非注释行
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped and not stripped.startswith('#'):
                import_end = i
                break
    
    # 如果整个文件都是注释或空行，在末尾添加
    if import_end == 0 and all(not line.strip() or line.strip().startswith('#') for line in lines):
        import_end = len(lines)
    
    # 创建新内容
    new_lines = []
    
    # 添加import部分
    new_lines.extend(lines[:import_end])
    
    # 添加main函数定义
    new_lines.append('')
    new_lines.append('def main():')
    new_lines.append('    """主函数，包装原有执行代码"""')
    
    # 添加原有代码（缩进）
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
    
    return '\n'.join(new_lines)

def fix_script(script_path):
    """修复单个脚本"""
    print(f"处理: {script_path.name}")
    
    # 读取内容
    try:
        with open(script_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError:
        try:
            with open(script_path, 'r', encoding='gbk') as f:
                content = f.read()
        except UnicodeDecodeError:
            with open(script_path, 'r', encoding='latin-1') as f:
                content = f.read()
    
    # 检查是否需要修复
    if not needs_main_guard(content):
        print(f"  ⏭️  已有main guard或无需修复")
        return False
    
    # 备份
    backup_path = backup_script(script_path)
    
    # 修复
    new_content = add_main_guard_safely(content)
    
    # 写入
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    print(f"  ✅ 添加main guard保护")
    print(f"  备份: {backup_path}")
    return True

def main():
    print("=== 安全修复：添加if __name__ == '__main__':保护 ===\n")
    
    # 获取所有Python脚本（排除自身）
    scripts_dir = Path(__file__).parent
    all_scripts = list(scripts_dir.glob("*.py"))
    scripts_to_fix = [s for s in all_scripts if s.name != __file__]
    
    print(f"找到 {len(scripts_to_fix)} 个Python脚本")
    
    # 询问修复范围
    print("\n选择修复范围:")
    print("1. 测试修复（前5个脚本）")
    print("2. 全部修复")
    print("3. 仅检查，不修复")
    
    choice = input("\n请选择 (1-3, 默认1): ").strip()
    
    if choice == '2':
        scripts_to_process = scripts_to_fix
        print("修复全部脚本")
    elif choice == '3':
        # 仅检查
        print("仅检查模式...")
        for script in scripts_to_fix[:10]:  # 只检查前10个
            with open(script, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            needs = needs_main_guard(content)
            has_guard = 'if __name__ == "__main__":' in content
            status = "需要修复" if needs else ("已有保护" if has_guard else "无需保护")
            print(f"{script.name}: {status}")
        return 0
    else:  # 默认或选择1
        scripts_to_process = scripts_to_fix[:5]
        print(f"测试修复: 处理前5个脚本")
    
    # 确认
    confirm = input("\n确认开始修复？(y/N): ").strip().lower()
    if confirm != 'y':
        print("取消修复")
        return 0
    
    # 执行修复
    fixed_count = 0
    for script_path in scripts_to_process:
        try:
            if fix_script(script_path):
                fixed_count += 1
        except Exception as e:
            print(f"  ❌ 修复失败: {e}")
    
    # 汇总
    print(f"\n=== 修复完成 ===")
    print(f"处理脚本数: {len(scripts_to_process)}")
    print(f"成功修复: {fixed_count}")
    print(f"无需修复: {len(scripts_to_process) - fixed_count}")
    
    if fixed_count > 0:
        print(f"\n备份文件保存在: {scripts_dir / 'backups' / 'main_guard_fixes'}")
        print("建议测试修复后的脚本是否正常工作")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())