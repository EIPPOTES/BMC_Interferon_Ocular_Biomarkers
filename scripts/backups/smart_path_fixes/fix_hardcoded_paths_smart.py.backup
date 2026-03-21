#!/usr/bin/env python3
"""
智能硬编码路径修复工具
避免替换文档字符串和注释中的路径
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
        RAW_DATA_DIR, VISUALIZATIONS_DIR,
        BASE_DIR
    )
    CONFIG_AVAILABLE = True
except ImportError as e:
    print(f"❌ 无法导入路径配置: {e}")
    print("请先运行 configs.paths_config.py")
    CONFIG_AVAILABLE = False
    sys.exit(1)

def is_in_string_or_comment(content, position):
    """检查给定位置是否在字符串或注释中"""
    # 简化检查：查找前面的引号或注释符号
    lines = content[:position].split('\n')
    current_line = lines[-1]
    
    # 检查是否在注释中
    if '#' in current_line:
        comment_pos = current_line.find('#')
        char_pos = len(current_line) - len(lines[-1])
        if position - (position - len(current_line) + char_pos) > comment_pos:
            return True
    
    # 检查是否在字符串中（简化版）
    # 统计当前行中未匹配的引号
    single_quotes = current_line.count("'") - current_line.count("\\'")
    double_quotes = current_line.count('"') - current_line.count('\\"')
    
    if single_quotes % 2 == 1 or double_quotes % 2 == 1:
        return True
    
    return False

def smart_replace_hardcoded_paths(content):
    """智能替换硬编码路径，避免文档字符串和注释"""
    if not CONFIG_AVAILABLE:
        return content, []
    
    replacements = []
    
    # 定义替换规则（更精确的匹配）
    replace_rules = [
        # 精确的主数据文件路径（带引号）
        (r'(["\'])/root/\.openclaw/workspace/data\.xlsx(["\'"])', 
         lambda m: f'{m.group(1)}str({MAIN_DATA_FILE}){m.group(2)}'),
        
        (r'(["\'])\./data\.xlsx(["\'])', 
         lambda m: f'{m.group(1)}str({MAIN_DATA_FILE}){m.group(2)}'),
        
        (r'(["\'])data\.xlsx(["\'])', 
         lambda m: f'{m.group(1)}str({MAIN_DATA_FILE}){m.group(2)}'),
        
        # 精确的Windows路径
        (r'(["\'])/mnt/c/Users/CUI/Desktop/论文及图表(["\'])', 
         lambda m: f'{m.group(1)}str({FIGURES_OUTPUT_DIR}){m.group(2)}'),
        
        (r'(["\'])/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/(["\'])', 
         lambda m: f'{m.group(1)}str({TABLES_OUTPUT_DIR}){m.group(2)}'),
        
        (r'(["\'])/mnt/c/Users/CUI/Desktop/final/02_Figures/(["\'])', 
         lambda m: f'{m.group(1)}str({FIGURES_OUTPUT_DIR}){m.group(2)}'),
    ]
    
    new_content = content
    
    for pattern, replacement_func in replace_rules:
        # 查找所有匹配
        matches = list(re.finditer(pattern, new_content))
        for match in reversed(matches):  # 从后往前替换，避免位置偏移
            start, end = match.start(), match.end()
            
            # 检查是否在字符串或注释中（应该是在字符串中，否则不匹配）
            # 这里我们信任正则表达式只匹配引号内的内容
            old_text = match.group(0)
            new_text = replacement_func(match)
            
            if old_text != new_text:
                replacements.append((old_text, new_text))
                new_content = new_content[:start] + new_text + new_content[end:]
    
    return new_content, replacements

def add_config_import_smart(content):
    """智能添加配置导入，避免重复"""
    # 检查是否已经导入了配置
    config_imports = [
        'from configs.paths_config import',
        'import configs.paths_config',
    ]
    
    for import_stmt in config_imports:
        if import_stmt in content:
            return content, False  # 已存在导入
    
    lines = content.split('\n')
    
    # 找到文件开头的注释和shebang
    insert_index = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('#!'):
            insert_index = i + 1
        elif stripped and not stripped.startswith('#'):
            # 第一个非注释行
            break
    
    # 创建新内容
    new_lines = []
    
    # 添加shebang和初始注释
    new_lines.extend(lines[:insert_index])
    
    # 添加配置导入
    if insert_index > 0:
        new_lines.append('')
    new_lines.append('# 路径配置')
    new_lines.append('from pathlib import Path')
    new_lines.append('import sys')
    new_lines.append('sys.path.insert(0, str(Path(__file__).parent.parent))')
    new_lines.append('from configs.paths_config import *')
    
    # 添加剩余内容
    new_lines.extend(lines[insert_index:])
    
    new_content = '\n'.join(new_lines)
    return new_content, True

def fix_script_smart(script_path, dry_run=False):
    """智能修复单个脚本"""
    print(f"处理: {script_path.name}")
    
    # 备份
    if not dry_run:
        backup_dir = current_dir / "backups" / "smart_path_fixes"
        backup_dir.mkdir(parents=True, exist_ok=True)
        backup_path = backup_dir / f"{script_path.name}.backup"
        shutil.copy2(script_path, backup_path)
    
    # 读取内容
    with open(script_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    original_content = content
    
    # 应用修复
    fixes_applied = []
    
    # 1. 添加配置导入
    content, import_added = add_config_import_smart(content)
    if import_added:
        fixes_applied.append('添加配置导入')
    
    # 2. 智能替换硬编码路径
    content, path_replacements = smart_replace_hardcoded_paths(content)
    if path_replacements:
        fixes_applied.append(f'替换{len(path_replacements)}个硬编码路径')
    
    # 检查是否有需要手动处理的路径模式
    manual_check_patterns = [
        r'/root/\.openclaw/workspace/[^"\']+',  # 未加引号的路径
        r'/mnt/c/Users/[^"\']+',  # 未加引号的Windows路径
        r'data\.xlsx(?!["\'])',  # 未加引号的data.xlsx
    ]
    
    manual_checks = []
    for pattern in manual_check_patterns:
        if re.search(pattern, content):
            manual_checks.append(pattern)
    
    if manual_checks and not dry_run:
        print(f"  ⚠️  发现需要手动检查的路径模式")
    
    # 写入修复后的内容
    if fixes_applied and not dry_run:
        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(content)
        
        print(f"  ✅ 应用修复: {', '.join(fixes_applied)}")
        if path_replacements:
            for i, (old_text, new_text) in enumerate(path_replacements[:2]):
                old_display = old_text if len(old_text) < 50 else old_text[:47] + '...'
                new_display = new_text if len(new_text) < 50 else new_text[:47] + '...'
                print(f"     {i+1}. {old_display} → {new_display}")
            if len(path_replacements) > 2:
                print(f"     ... 还有 {len(path_replacements)-2} 个替换")
        if not dry_run:
            print(f"  备份: {backup_path}")
        return True, len(path_replacements), manual_checks
    elif dry_run:
        if fixes_applied:
            print(f"  🔍 模拟修复: {', '.join(fixes_applied)}")
            if path_replacements:
                print(f"     将替换 {len(path_replacements)} 个路径")
        else:
            print(f"  ⏭️  无需修复")
        return False, 0, manual_checks
    else:
        print(f"  ⏭️  无需修复")
        return False, 0, manual_checks

def main():
    print("=== 智能硬编码路径修复工具 ===\n")
    
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
        if s.name not in [__file__, 'check_script_quality.py', 
                         'fix_main_guard_only.py', 'fix_script_quality.py',
                         'fix_hardcoded_paths.py', 'fix_create_figures.py']
    ]
    
    print(f"找到 {len(scripts_to_fix)} 个需要检查的脚本")
    
    # 询问修复模式
    print("\n选择修复模式:")
    print("1. 模拟运行（不修改文件）")
    print("2. 修复关键脚本（前10个）")
    print("3. 修复所有脚本")
    print("4. 仅检查analyze_data.py和create_figures.py")
    
    choice = input("\n请选择 (1-4, 默认1): ").strip()
    
    if choice == '2':
        scripts_to_process = scripts_to_fix[:10]
        dry_run = False
        print(f"修复关键脚本（前10个）")
    elif choice == '3':
        scripts_to_process = scripts_to_fix
        dry_run = False
        print(f"修复所有脚本")
    elif choice == '4':
        scripts_to_process = [s for s in scripts_to_fix if s.name in ['analyze_data.py', 'create_figures.py']]
        dry_run = False
        print(f"修复核心脚本: {[s.name for s in scripts_to_process]}")
    else:  # 默认或选择1
        scripts_to_process = scripts_to_fix[:5]
        dry_run = True
        print(f"模拟运行: 检查前5个脚本")
    
    if dry_run:
        print("⚠️ 模拟运行模式：不会实际修改文件")
    
    # 确认
    if not dry_run:
        confirm = input("\n确认开始修复？(y/N): ").strip().lower()
        if confirm != 'y':
            print("取消修复")
            return 0
    
    # 执行修复
    fixed_count = 0
    total_replacements = 0
    total_manual_checks = []
    
    for script_path in scripts_to_process:
        try:
            fixed, replacements, manual_checks = fix_script_smart(script_path, dry_run)
            if fixed or (dry_run and replacements > 0):
                fixed_count += 1
                total_replacements += replacements
            if manual_checks:
                total_manual_checks.append((script_path.name, manual_checks))
        except Exception as e:
            print(f"  ❌ 处理失败: {e}")
    
    # 汇总报告
    print(f"\n=== {'模拟' if dry_run else '修复'}完成 ===")
    print(f"处理脚本总数: {len(scripts_to_process)}")
    print(f"{'发现可修复' if dry_run else '成功修复'}: {fixed_count}")
    print(f"总替换数: {total_replacements}")
    
    if total_manual_checks:
        print(f"\n⚠️ 需要手动检查的脚本 ({len(total_manual_checks)}个):")
        for script_name, checks in total_manual_checks[:5]:
            print(f"  - {script_name}: {len(checks)}个可疑模式")
    
    if not dry_run and fixed_count > 0:
        print(f"\n备份文件保存在: {current_dir / 'backups' / 'smart_path_fixes'}")
        print("建议测试修复后的脚本是否正常工作")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())