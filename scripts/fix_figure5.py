#!/usr/bin/env python3
"""
修复generate_figures_english.py中的Figure 5，将硬编码数据改为动态计算
"""

import re

def add_effect_size_function(file_content):
    """在合适位置添加效应量计算函数"""
    
    # 寻找读取数据后的位置（在"Data loaded"行之后）
    lines = file_content.split('\n')
    insert_pos = -1
    
    for i, line in enumerate(lines):
        if 'df = pd.read_excel(data_file)' in line:
            # 在读取数据行之后插入函数
            insert_pos = i + 2  # 跳过空行
            break
    
    if insert_pos == -1:
        # 如果没有找到，在import部分之后插入
        for i, line in enumerate(lines):
            if 'warnings.filterwarnings' in line:
                insert_pos = i + 2
                break
    
    if insert_pos == -1:
        insert_pos = 20  # 默认位置
    
    # 要插入的函数代码
    effect_size_function = '''
def calculate_effect_sizes(df):
    """计算选定OCT参数的Cohen's d效应量"""
    from scipy import stats
    import numpy as np
    
    def cohens_d(group1, group2):
        """计算Cohen's d效应量"""
        n1, n2 = len(group1), len(group2)
        mean1, mean2 = np.mean(group1), np.mean(group2)
        var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
        
        # 合并标准差
        pooled_sd = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
        
        if pooled_sd == 0:
            return 0
        return (mean1 - mean2) / pooled_sd
    
    # 分组
    mdd_mask = df['分组'] == '抑郁症'
    control_mask = df['分组'] == '健康对照'
    
    if not control_mask.any():
        print("错误: 未找到对照组")
        return []
    
    mdd_group = df[mdd_mask]
    control_group = df[control_mask]
    
    print(f"  MDD组样本量: {len(mdd_group)}")
    print(f"  对照组样本量: {len(control_group)}")
    
    # 参数映射 (英文标签: 数据列名)
    param_mapping = [
        ('Mean Macular Thickness', 'Retina_平均厚度'),
        ('Total Macular Volume', 'Retina_总体积'),
        ('Outer Temporal', 'Retina_外环颞侧'),
        ('Inner Temporal', 'Retina_内环颞侧'),
        ('Superior RNFL', 'RNFL_上方'),
        ('Cup Area', 'Cup Area'),
        ('Rim Volume', 'Rim Volume'),
        ('C/D Area Ratio', 'C/D Area Ratio')
    ]
    
    results = []
    
    for eng_label, col_name in param_mapping:
        if col_name not in df.columns:
            print(f"警告: 列 '{col_name}' 不存在，跳过 {eng_label}")
            continue
        
        # 获取数据，移除NaN
        mdd_data = mdd_group[col_name].dropna()
        control_data = control_group[col_name].dropna()
        
        if len(mdd_data) < 5 or len(control_data) < 5:
            print(f"警告: {eng_label} 数据不足")
            continue
        
        # 计算效应量
        d = cohens_d(mdd_data, control_data)
        
        # 计算Mann-Whitney U检验p值
        try:
            u_stat, p_value = stats.mannwhitneyu(mdd_data, control_data, alternative='two-sided')
        except:
            # 如果Mann-Whitney失败，使用t检验
            t_stat, p_value = stats.ttest_ind(mdd_data, control_data, equal_var=False, nan_policy='omit')
        
        results.append((eng_label, d, p_value))
        
        print(f"    {eng_label}: d = {d:.3f}, p = {p_value:.6f}")
    
    return results
'''
    
    # 在指定位置插入函数
    lines.insert(insert_pos, effect_size_function)
    
    return '\n'.join(lines)

def update_figure5_section(file_content):
    """更新Figure 5部分，使用动态计算的数据"""
    
    # 查找Figure 5部分
    start_marker = '# ==================== FIGURE 5: Forest Plot ===================='
    end_marker = '# ==================== FIGURE 6: Subgroup Analysis ===================='
    
    start_idx = file_content.find(start_marker)
    if start_idx == -1:
        print("错误: 未找到Figure 5部分")
        return file_content
    
    end_idx = file_content.find(end_marker, start_idx)
    if end_idx == -1:
        print("错误: 未找到Figure 6部分")
        return file_content
    
    # 原始的Figure 5部分
    original_section = file_content[start_idx:end_idx]
    
    # 新的Figure 5部分
    new_section = '''# ==================== FIGURE 5: Forest Plot ====================
print("\\n【Figure 5: Forest Plot】")
print("  计算效应量...")

# 动态计算效应量
forest_data = calculate_effect_sizes(df)

if not forest_data:
    print("错误: 无法计算效应量，使用默认数据")
    # 备用数据（原始硬编码值）
    forest_data = [
        ('Mean Macular Thickness', -0.415, 0.000003),
        ('Total Macular Volume', -0.416, 0.000003),
        ('Outer Temporal', -0.497, 0.000003),
        ('Inner Temporal', -0.375, 0.000032),
        ('Superior RNFL', -0.311, 0.002229),
        ('Cup Area', 0.224, 0.022329),
        ('Rim Volume', -0.303, 0.010735),
        ('C/D Area Ratio', 0.246, 0.021236),
    ]

fig, ax = plt.subplots(figsize=(10, 8))

labels = [d[0] for d in forest_data]
effects = [d[1] for d in forest_data]
p_values = [d[2] for d in forest_data]

colors = ['red' if e < 0 else 'blue' for e in effects]

# 绘制森林图
y_pos = np.arange(len(labels))
bars = ax.barh(y_pos, effects, color=colors, alpha=0.7, edgecolor='black')

# 添加参考线
ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
ax.axvline(x=-0.2, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
ax.axvline(x=0.2, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

# 添加P值标记
for i, (effect, p) in enumerate(zip(effects, p_values)):
    if p < 0.001:
        sig = '***'
    elif p < 0.01:
        sig = '**'
    elif p < 0.05:
        sig = '*'
    else:
        sig = ''
    ax.text(effect + 0.02 if effect > 0 else effect - 0.02, i, 
            f'{effect:.2f} {sig}', va='center', fontsize=9, ha='left' if effect > 0 else 'right')

ax.set_yticks(y_pos)
ax.set_yticklabels(labels)
ax.set_xlabel("Cohen's d (Effect Size)", fontsize=12)
ax.set_title('Figure 5. Effect Sizes of Retinal Changes in MDD Patients', fontsize=14, fontweight='bold')
ax.grid(True, axis='x', alpha=0.3)

# 添加图例
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='red', alpha=0.7, label='Reduced in MDD'),
                   Patch(facecolor='blue', alpha=0.7, label='Increased in MDD')]
ax.legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
plt.savefig(f'{figures_dir}/Figure5_Forest_Plot.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"✓ Figure 5 saved")
'''
    
    # 替换原始部分
    updated_content = file_content[:start_idx] + new_section + file_content[end_idx:]
    
    return updated_content

def main():
    input_file = '/root/.openclaw/workspace/generate_figures_english.py'
    backup_file = '/root/.openclaw/workspace/generate_figures_english.py.backup'
    
    try:
        # 读取原始文件
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        print(f"读取文件: {input_file}")
        
        # 创建备份
        import shutil
        shutil.copy2(input_file, backup_file)
        print(f"创建备份: {backup_file}")
        
        # 第一步：添加效应量计算函数
        print("添加效应量计算函数...")
        content = add_effect_size_function(content)
        
        # 第二步：更新Figure 5部分
        print("更新Figure 5部分...")
        content = update_figure5_section(content)
        
        # 写入更新后的文件
        output_file = '/root/.openclaw/workspace/generate_figures_english_fixed.py'
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(content)
        
        print(f"修复完成! 输出文件: {output_file}")
        print(f"原始文件已备份到: {backup_file}")
        
        # 测试新文件
        print("\n测试新文件...")
        print("=" * 60)
        
        # 导入并运行测试
        import subprocess
        result = subprocess.run(['python3', output_file], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("测试成功! 新文件运行正常。")
            # 显示Figure 5相关的输出
            for line in result.stdout.split('\\n'):
                if 'Figure 5' in line or '效应量' in line or 'd =' in line:
                    print(line)
        else:
            print(f"测试失败! 错误信息:")
            print(result.stderr)
            
    except Exception as e:
        print(f"修复过程中出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()