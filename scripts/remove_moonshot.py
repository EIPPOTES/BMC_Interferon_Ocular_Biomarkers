#!/usr/bin/env python3
# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
直接从openclaw.json中移除moonshot-v1-128k模型
"""

import json
import sys
import os

def remove_moonshot_model():
    config_file = "/root/.openclaw/openclaw.json"
    
    # 备份原文件
    backup_file = config_file + ".backup"
    os.system(f"cp {config_file} {backup_file}")
    print(f"创建备份: {backup_file}")
    
    # 读取配置文件
    with open(config_file, 'r', encoding='utf-8') as f:
        # 使用json5解析，因为文件使用js5风格（有注释）
        try:
            import json5
            config = json5.load(f)
        except ImportError:
            # 如果没有json5，尝试用json解析（移除注释）
            content = f.read()
            # 简单移除单行注释
            lines = []
            for line in content.split('\n'):
                if '//' in line:
                    line = line.split('//')[0]
                lines.append(line)
            content = '\n'.join(lines)
            config = json.loads(content)
    
    print(f"原始配置已读取")
    
    # 修改配置
    # 1. 从models.providers.moonshot.models中移除moonshot-v1-128k
    moonshot_models = config.get('models', {}).get('providers', {}).get('moonshot', {}).get('models', [])
    
    if moonshot_models:
        new_models = [m for m in moonshot_models if m.get('id') != 'moonshot-v1-128k']
        print(f"过滤模型: {len(moonshot_models)} -> {len(new_models)}")
        
        if len(new_models) == 0:
            print("警告: moonshot provider将没有模型")
        elif len(new_models) == 1 and new_models[0].get('id') == 'moonshot-v1-32k':
            print("只保留moonshot-v1-32k模型")
        
        config['models']['providers']['moonshot']['models'] = new_models
    
    # 2. 确保models.mode是'replace'
    if config.get('models', {}).get('mode') != 'replace':
        config['models']['mode'] = 'replace'
        print("设置models.mode = 'replace'")
    
    # 写入文件
    with open(config_file, 'w', encoding='utf-8') as f:
        # 使用json5格式写入
        try:
            import json5
            json5.dump(config, f, indent=2, quote_keys=True, trailing_commas=False)
        except ImportError:
            json.dump(config, f, indent=2)
    
    print(f"配置已更新: {config_file}")
    
    # 验证更改
    with open(config_file, 'r', encoding='utf-8') as f:
        content = f.read()
        if 'moonshot-v1-128k' in content:
            print("警告: 文件中可能仍然包含moonshot-v1-128k")
        else:
            print("✅ moonshot-v1-128k已从文件中移除")
    
    return True

if __name__ == "__main__":
    try:
        remove_moonshot_model()
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)