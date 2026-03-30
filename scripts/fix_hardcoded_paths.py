#!/usr/bin/env python3
"""
硬编码路径修复工具
批量替换Python脚本中的硬编码路径为配置引用
"""

import os
import sys
import re
import shutil
from pathlib import Path

# 设置工作目录
current_dir = Path(__file__).parent
workspace_dir = current_dir.parent

# 导入路径配置
sys.path.insert(0, str(workspace_dir))
try:
    from configs.paths_config import (
        MAIN_DATA_FILE, FIGURES_OUTPUT_DIR,
        TABLES_OUTPUT_DIR, REPORTS_OUTPUT_DIR,
        RAW_DATA_DIR, VISUALIZATIONS_DIR
    )
    CONFIG_AVAILABLE = True
except ImportError as e:
    print(f"❌ 无法导入路径配置: {e}")
    print("请先运行 configs/paths_config.py")
    CONFIG_AVAILABLE = False
    sys.exit(1)

def backup_script(script_path):
    """备份脚本文件"""
    backup_dir = current_dir / "backups" / "hardcoded_path_fixes"
    backup_dir.mkdir(parents=True, exist_ok=True)
    
    backup_path = backup_dir / f"{script_path.name}.backup"
    shutil.copy2(script_path, backup_path)
    return backup_path

def read_script_content(script_path):
    """读取脚本内容，处理编码问题"""
    encodings = ['utf-8', 'gbk', 'latin-1']
    
    for encoding in encodings:
        try:
            with open(script_path, 'r', encoding=encoding) as f:
                return f.read()
        except UnicodeDecodeError:
            continue
    
    # 如果所有编码都失败，使用二进制读取
    with open(script_path, 'rb') as f:
        return f.read().decode('utf-8', errors='ignore')

def write_script_content(script_path, content):
    """写入脚本内容"""
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(content)

def add_config_import(content):
    """在脚本开头添加配置导入"""
    # 检查是否已经导入了配置
    config_imports = [
        'from configs.paths_config import',
        'import configs.paths_config',
        'sys.path.insert(0, str(Path(__file__).parent.parent))'
    ]
    
    for import_stmt in config_imports:
        if import_stmt in content:
            return content, False  # 已存在导入
    
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
    
    # 创建新内容
    new_lines = []
    
    # 添加import部分
    new_lines.extend(lines[:import_end])
    
    # 添加配置导入
    new_lines.append('# 路径配置')
    new_lines.append('from pathlib import Path')
    new_lines.append('import sys')
    new_lines.append('sys.path.insert(0, str(Path(__file__).parent.parent))')
    new_lines.append('from configs.paths_config import *')
    
    # 添加剩余内容
    new_lines.extend(lines[import_end:])
    
    new_content = '\n'.join(new_lines)
    return new_content, True

def replace_hardcoded_paths(content):
    """替换硬编码路径为配置引用"""
    if not CONFIG_AVAILABLE:
        return content, []
    
    replacements = []
    
    # 定义替换规则
    replace_rules = [
        # 主数据文件路径
        (r'"/root/\.openclaw/workspace/data\.xlsx"', f'str({MAIN_DATA_FILE})'),
        (r"'/root/\.openclaw/workspace/data\.xlsx'", f'str({MAIN_DATA_FILE})'),
        (r'"/root/\.openclaw/workspace/data/raw/data\.xlsx"', f'str({MAIN_DATA_FILE})'),
        (r"'/root/\.openclaw/workspace/data/raw/data\.xlsx'", f'str({MAIN_DATA_FILE})'),
        (r'"\./data\.xlsx"', f'str({MAIN_DATA_FILE})'),
        (r"'\./data\.xlsx'", f'str({MAIN_DATA_FILE})'),
        (r'data\.xlsx', f'str({MAIN_DATA_FILE})'),
        
        # 图表输出路径
        (r'"/mnt/c/Users/CUI/Desktop/论文及图表"', f'str({FIGURES_OUTPUT_DIR})'),
        (r"'/mnt/c/Users/CUI/Desktop/论文及图表'", f'str({FIGURES_OUTPUT_DIR})'),
        (r'"/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/"', f'str({TABLES_OUTPUT_DIR})'),
        (r"'/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/'", f'str({TABLES_OUTPUT_DIR})'),
        (r'"/mnt/c/Users/CUI/Desktop/final/02_Figures/"', f'str({FIGURES_OUTPUT_DIR})'),
        (r"'/mnt/c/Users/CUI/Desktop/final/02_Figures/'", f'str({FIGURES_OUTPUT_DIR})'),
        
        # 通用Windows路径模式
        (r'"/mnt/c/Users/CUI/Desktop/[^"]*"', f'str({FIGURES_OUTPUT_DIR})'),
        (r"'/mnt/c/Users/CUI/Desktop/[^']*'", f'str({FIGURES_OUTPUT_DIR})'),
        
        # 工作空间根目录
        (r'"/root/\.openclaw/workspace/"', 'str(BASE_DIR)'),
        (r"'/root/\.openclaw/workspace/'", 'str(BASE_DIR)'),
    ]
    
    new_content = content
    
    for pattern, replacement in replace_rules:
        # 使用正则表达式查找匹配
        matches = re.findall(pattern, new_content)
        if matches:
            # 记录替换
            replacements.extend([(pattern, replacement, match) for match in matches])
            # 执行替换
            new_content = re.sub(pattern, replacement, new_content)
    
    return new_content, replacements

def fix_script(script_path):
    """修复单个脚本的硬编码路径"""
    print(f"处理: {script_path.name}")
    
    # 备份
    backup_path = backup_script(script_path)
    
    # 读取内容
    content = read_script_content(script_path)
    original_content = content
    
    # 应用修复
    fixes_applied = []
    
    # 1. 添加配置导入
    content, import_added = add_config_import(content)
    if import_added:
        fixes_applied.append('添加配置导入')
    
    # 2. 替换硬编码路径
    content, path_replacements = replace_hardcoded_paths(content)
    if path_replacements:
        fixes_applied.append(f'替换{len(path_replacements)}个硬编码路径')
    
    # 写入修复后的内容
    if fixes_applied:
        write_script_content(script_path, content)
        print(f"  ✅ 应用修复: {', '.join(fixes_applied)}")
        if path_replacements:
            for i, (pattern, replacement, match) in enumerate(path_replacements[:3]):
                print(f"     {i+1}. {match[:50]}... → {replacement}")
            if len(path_replacements) > 3:
                print(f"     ... 还有 {len(path_replacements)-3} 个替换")
        print(f"  备份: {backup_path}")
        return True, len(path_replacements)
    else:
        print(f"  ⏭️  无需修复")
        # 删除备份文件，因为没有修改
        backup_path.unlink(missing_ok=True)
        return False, 0

def main():
    print("=== 硬编码路径批量修复工具 ===\n")
    
    if not CONFIG_AVAILABLE:
        print("❌ 路径配置不可用，无法继续")
        return 1
    
    print("路径配置状态:")
    print(f"  主数据文件: {MAIN_DATA_FILE}")
    print(f"  图表输出目录: {FIGURES_OUTPUT_DIR}")
    print(f"  表格输出目录: {TABLES_OUTPUT_DIR}")
    print(f"  报告输出目录: {REPORTS_OUTPUT_DIR}")
    print()
    
    # 获取所有Python脚本（排除自身和检查脚本）
    all_scripts = list(current_dir.glob("*.py"))
    scripts_to_fix = [
        s for s in all_scripts 
        if s.name not in [__file__, 'check_script_quality.py', 'fix_main_guard_only.py', 'fix_script_quality.py']
    ]
    
    print(f"找到 {len(scripts_to_fix)} 个需要检查的脚本")
    
    # 询问修复范围
    print("\n选择修复范围:")
    print("1. 测试修复（前5个脚本）")
    print("2. 修复前50个脚本")
    print("3. 修复所有脚本")
    print("4. 仅检查，不修复")
    
    choice = input("\n请选择 (1-4, 默认1): ").strip()
    
    if choice == '2':
        scripts_to_process = scripts_to_fix[:50]
        print(f"修复前50个脚本")
    elif choice == '3':
        scripts_to_process = scripts_to_fix
        print(f"修复所有脚本")
    elif choice == '4':
        # 仅检查模式
        print("仅检查模式...")
        print("扫描常见的硬编码路径模式:")
        
        path_patterns = [
            r'/root/\.openclaw/workspace/',
            r'/mnt/c/Users/CUI/Desktop/',
            r'\./data\.xlsx',
            r'data\.xlsx'
        ]
        
        problematic_scripts = []
        for script in scripts_to_fix[:20]:  # 只检查前20个
            content = read_script_content(script)
            has_hardcoded = False
            for pattern in path_patterns:
                if re.search(pattern, content):
                    has_hardcoded = True
                    break
            
            if has_hardcoded:
                problematic_scripts.append(script.name)
        
        print(f"发现 {len(problematic_scripts)} 个脚本有硬编码路径:")
        for script in problematic_scripts[:20]:
            print(f"  - {script}")
        if len(problematic_scripts) > 20:
            print(f"  ... 还有 {len(problematic_scripts)-20} 个")
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
    total_replacements = 0
    for script_path in scripts_to_process:
        try:
            fixed, replacements = fix_script(script_path)
            if fixed:
                fixed_count += 1
                total_replacements += replacements
        except Exception as e:
            print(f"  ❌ 修复失败: {e}")
    
    # 汇总报告
    print(f"\n=== 修复完成 ===")
    print(f"处理脚本总数: {len(scripts_to_process)}")
    print(f"成功修复: {fixed_count}")
    print(f"总替换数: {total_replacements}")
    print(f"无需修复: {len(scripts_to_process) - fixed_count}")
    
    if fixed_count > 0:
        print(f"\n备份文件保存在: {current_dir / 'backups' / 'hardcoded_path_fixes'}")
        print("建议测试修复后的脚本是否正常工作")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())