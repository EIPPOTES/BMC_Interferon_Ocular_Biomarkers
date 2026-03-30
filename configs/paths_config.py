"""
路径配置文件 - 统一管理所有文件路径
避免硬编码路径，提高可维护性和可移植性
"""

import os
from pathlib import Path

# 基础路径
BASE_DIR = Path(__file__).parent.parent  # workspace根目录

# 数据目录
DATA_DIR = BASE_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
VISUALIZATIONS_DIR = DATA_DIR / "visualizations"

# 脚本目录
SCRIPTS_DIR = BASE_DIR / "scripts"

# 报告目录
REPORTS_DIR = BASE_DIR / "reports"

# 配置目录
CONFIGS_DIR = BASE_DIR / "configs"

# 输出目录（根据环境变量可配置）
# 默认使用Windows桌面路径（WSL兼容）
DEFAULT_OUTPUT_DIR = Path("/mnt/c/Users/CUI/Desktop/论文及图表")
FIGURES_OUTPUT_DIR = Path(os.environ.get("FIGURES_OUTPUT_DIR", DEFAULT_OUTPUT_DIR))
TABLES_OUTPUT_DIR = Path(os.environ.get("TABLES_OUTPUT_DIR", DEFAULT_OUTPUT_DIR))
REPORTS_OUTPUT_DIR = Path(os.environ.get("REPORTS_OUTPUT_DIR", DEFAULT_OUTPUT_DIR))

# 关键文件路径
MAIN_DATA_FILE = RAW_DATA_DIR / "data.xlsx"  # 主数据文件

# 外部数据源目录（用于需要访问原始数据文件的脚本）
EXTERNAL_DATA_DIR = Path("/mnt/c/Users/CUI/Desktop/最终修改")

# 确保所有目录存在
def ensure_directories():
    """确保所有必要的目录都存在"""
    directories = [
        RAW_DATA_DIR,
        PROCESSED_DATA_DIR,
        VISUALIZATIONS_DIR,
        REPORTS_DIR,
        CONFIGS_DIR,
        FIGURES_OUTPUT_DIR,
        TABLES_OUTPUT_DIR,
        REPORTS_OUTPUT_DIR,
    ]
    
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)
        print(f"✓ 目录就绪: {directory}")
    
    return True

# 路径验证
def validate_paths():
    """验证关键路径是否存在"""
    errors = []
    
    # 检查数据文件
    if not MAIN_DATA_FILE.exists():
        errors.append(f"数据文件不存在: {MAIN_DATA_FILE}")
    
    # 检查目录
    if not SCRIPTS_DIR.exists():
        errors.append(f"脚本目录不存在: {SCRIPTS_DIR}")
    
    if not REPORTS_DIR.exists():
        errors.append(f"报告目录不存在: {REPORTS_DIR}")
    
    if not CONFIGS_DIR.exists():
        errors.append(f"配置目录不存在: {CONFIGS_DIR}")
    
    if errors:
        print("⚠️ 路径验证发现问题:")
        for error in errors:
            print(f"  - {error}")
        return False
    
    print("✅ 所有关键路径验证通过")
    return True

# 工具函数
def get_relative_path(absolute_path):
    """获取相对于工作空间的相对路径"""
    try:
        return Path(absolute_path).relative_to(BASE_DIR)
    except ValueError:
        return absolute_path

def resolve_path(path):
    """解析路径，支持相对路径和绝对路径"""
    path_obj = Path(path)
    if path_obj.is_absolute():
        return path_obj
    else:
        return BASE_DIR / path

if __name__ == "__main__":
    # 测试路径配置
    print("=== 路径配置测试 ===")
    ensure_directories()
    validate_paths()
    print(f"工作空间根目录: {BASE_DIR}")
    print(f"主数据文件: {MAIN_DATA_FILE}")
    print(f"图表输出目录: {FIGURES_OUTPUT_DIR}")
    print("=== 测试完成 ===")