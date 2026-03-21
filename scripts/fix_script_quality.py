#!/usr/bin/env python3
"""
脚本质量修复工具
修复Python脚本的常见质量问题：
1. 添加if __name__ == "__main__":保护
2. 替换硬编码路径为配置引用
3. 重构顶层执行代码到main()函数
"""

import os
import sys
import re
import shutil
from pathlib import Path

# 设置工作目录
current_dir = Path(__file__).parent
configs_dir = current_dir.parent / "configs"

# 导入路径配置
sys.path.insert(0, str(configs_dir.parent))
try:
    from configs.paths_config import (
        MAIN_DATA_FILE, FIGURES_OUTPUT_DIR, 
        TABLES_OUTPUT_DIR, REPORTS_OUTPUT_DIR,
        ensure_directories
    )
    CONFIG_AVAILABLE = True
except ImportError:
    print("⚠️ 无法导入路径配置，请先运行 configs/paths_config.py")
    CONFIG_AVAILABLE = False

def backup_script(script_path):
    """备份脚本文件"""
    backup_dir = current_dir / "backups" / "quality_fixes"
    backup_dir.mkdir(parents=True, exist_ok=True)
    
    backup_path = backup_dir / f"{script_path.name}.backup"
    shutil.copy2(script_path, backup_path)
    return backup_path

def read_script_content(script_path):
    """读取脚本内容，处理编码问题"""
    try:
        with open(script_path, 'r', encoding='utf-8') as f:
            return f.read()
    except UnicodeDecodeError:
        try:
            with open(script_path, 'r', encoding='gbk') as f:
                return f.read()
        except UnicodeDecodeError:
            with open(script_path, 'r', encoding='latin-1') as f:
                return f.read()

def write_script_content(script_path, content):
    """写入脚本内容"""
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(content)

def fix_main_guard(script_path, content):
    """为脚本添加if __name__ == "__main__":保护"""
    if 'if __name__ == "__main__":' in content:
        return content, False  # 已存在，不需要修复
    
    # 分析脚本结构，找到顶层代码
    lines = content.split('\n')
    
    # 找到import语句的结束位置
    import_end = 0
    for i, line in enumerate(lines):
        line_strip = line.strip()
        if line_strip and not line_strip.startswith('#') and not line_strip.startswith('import') and not line_strip.startswith('from'):
            import_end = i
            break
    
    # 如果import_end为0，说明没有import语句
    if import_end == 0:
        import_end = len([l for l in lines if l.strip().startswith(('import', 'from'))])
    
    # 将剩余代码包装到main()函数中
    new_lines = lines[:import_end]
    
    # 添加main函数定义
    new_lines.append('')
    new_lines.append('def main():')
    
    # 添加缩进的代码
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
    
    new_content = '\n'.join(new_lines)
    return new_content, True

def fix_hardcoded_paths(script_path, content):
    """替换硬编码路径为配置引用"""
    if not CONFIG_AVAILABLE:
        return content, False
    
    replacements = []
    
    # 替换数据文件路径
    data_patterns = [
        (r'"/root/\.openclaw/workspace/data\.xlsx"', f'str({MAIN_DATA_FILE})'),
        (r"'/root/\.openclaw/workspace/data\.xlsx'", f'str({MAIN_DATA_FILE})'),
        (r'"\./data\.xlsx"', f'str({MAIN_DATA_FILE})'),
        (r"'\./data\.xlsx'", f'str({MAIN_DATA_FILE})'),
        (r'data\.xlsx', f'str({MAIN_DATA_FILE})'),
    ]
    
    # 替换输出路径
    output_patterns = [
        (r'"/mnt/c/Users/CUI/Desktop/论文及图表"', f'str({FIGURES_OUTPUT_DIR})'),
        (r"'/mnt/c/Users/CUI/Desktop/论文及图表'", f'str({FIGURES_OUTPUT_DIR})'),
        (r'"/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/"', f'str({TABLES_OUTPUT_DIR})'),
        (r"'/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/'", f'str({TABLES_OUTPUT_DIR})'),
    ]
    
    new_content = content
    
    for pattern, replacement in data_patterns + output_patterns:
        matches = re.findall(pattern, new_content)
        if matches:
            new_content = re.sub(pattern, replacement, new_content)
            replacements.extend(matches)
    
    return new_content, len(replacements) > 0

def add_config_import(script_path, content):
    """在脚本开头添加配置导入"""
    if not CONFIG_AVAILABLE:
        return content, False
    
    # 检查是否已经导入了配置
    config_import = 'from configs.paths_config import'
    if config_import in content:
        return content, False
    
    # 找到import语句的位置
    lines = content.split('\n')
    import_insert_index = 0
    
    for i, line in enumerate(lines):
        if line.strip().startswith(('import', 'from')):
            import_insert_index = i + 1
        elif line.strip() and not line.strip().startswith('#'):
            # 遇到非import/非注释行，停止查找
            break
    
    # 插入配置导入
    new_lines = lines[:import_insert_index]
    new_lines.append('# 路径配置')
    new_lines.append('from pathlib import Path')
    new_lines.append('import sys')
    new_lines.append('sys.path.insert(0, str(Path(__file__).parent.parent))')
    new_lines.append('from configs.paths_config import *')
    new_lines.extend(lines[import_insert_index:])
    
    new_content = '\n'.join(new_lines)
    return new_content, True

def fix_script(script_path, fix_types=None):
    """修复单个脚本的质量问题"""
    if fix_types is None:
        fix_types = ['main_guard', 'paths', 'import']
    
    print(f"修复: {script_path.name}")
    
    # 备份原始脚本
    backup_path = backup_script(script_path)
    
    # 读取内容
    content = read_script_content(script_path)
    original_content = content
    
    # 应用修复
    fixes_applied = []
    
    if 'import' in fix_types:
        content, fixed = add_config_import(script_path, content)
        if fixed:
            fixes_applied.append('添加配置导入')
    
    if 'paths' in fix_types:
        content, fixed = fix_hardcoded_paths(script_path, content)
        if fixed:
            fixes_applied.append('替换硬编码路径')
    
    if 'main_guard' in fix_types:
        content, fixed = fix_main_guard(script_path, content)
        if fixed:
            fixes_applied.append('添加main guard')
    
    # 写入修复后的内容
    if fixes_applied:
        write_script_content(script_path, content)
        print(f"  ✅ 应用修复: {', '.join(fixes_applied)}")
        print(f"  备份: {backup_path}")
        return True
    else:
        print(f"  ⏭️  无需修复")
        # 删除备份文件，因为没有修改
        backup_path.unlink(missing_ok=True)
        return False

def main():
    print("=== Python脚本质量修复工具 ===\n")
    
    if not CONFIG_AVAILABLE:
        print("⚠️ 警告: 路径配置不可用，路径修复功能将受限")
        print("请先运行 configs/paths_config.py 确保配置可用")
        print("继续仅修复main guard问题...\n")
    
    # 获取所有Python脚本（排除自身和检查脚本）
    all_scripts = list(current_dir.glob("*.py"))
    scripts_to_fix = [
        s for s in all_scripts 
        if s.name not in [__file__, 'check_script_quality.py']
    ]
    
    print(f"找到 {len(scripts_to_fix)} 个需要检查的脚本")
    
    # 询问修复类型
    print("\n选择修复类型:")
    print("1. 仅添加main guard保护")
    print("2. 仅修复硬编码路径")
    print("3. 仅添加配置导入")
    print("4. 全面修复（推荐）")
    print("5. 测试修复（仅处理前5个脚本）")
    
    choice = input("\n请选择 (1-5, 默认4): ").strip()
    
    if choice == '1':
        fix_types = ['main_guard']
        scripts_to_fix = scripts_to_fix[:]  # 所有脚本
    elif choice == '2':
        fix_types = ['paths']
        scripts_to_fix = scripts_to_fix[:]  # 所有脚本
    elif choice == '3':
        fix_types = ['import']
        scripts_to_fix = scripts_to_fix[:]  # 所有脚本
    elif choice == '5':
        fix_types = ['main_guard', 'paths', 'import']
        scripts_to_fix = scripts_to_fix[:5]  # 仅前5个
        print(f"测试模式: 仅处理前5个脚本")
    else:  # 默认或选择4
        fix_types = ['main_guard', 'paths', 'import']
        scripts_to_fix = scripts_to_fix[:]  # 所有脚本
    
    print(f"修复类型: {fix_types}")
    print(f"处理脚本数: {len(scripts_to_fix)}")
    
    # 确认
    confirm = input("\n确认开始修复？(y/N): ").strip().lower()
    if confirm != 'y':
        print("取消修复")
        return 0
    
    # 执行修复
    fixed_count = 0
    for script_path in scripts_to_fix:
        try:
            if fix_script(script_path, fix_types):
                fixed_count += 1
        except Exception as e:
            print(f"  ❌ 修复失败: {e}")
    
    # 汇总报告
    print(f"\n=== 修复完成 ===")
    print(f"处理脚本总数: {len(scripts_to_fix)}")
    print(f"成功修复: {fixed_count}")
    print(f"无需修复: {len(scripts_to_fix) - fixed_count}")
    
    if fixed_count > 0:
        print(f"\n备份文件保存在: {current_dir / 'backups' / 'quality_fixes'}")
        print("建议运行 check_script_quality.py 验证修复效果")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())