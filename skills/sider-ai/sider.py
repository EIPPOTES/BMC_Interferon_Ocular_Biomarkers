#!/usr/bin/env python3
"""
sider-ai Python Wrapper
通过 sider.ai API 调用 GPT、Claude、Gemini、DeepSeek、o1 等全部模型
"""

import os
import json
import sys
from pathlib import Path

# 配置路径
CONFIG_DIR = Path.home() / ".openclaw" / "workspace" / "skills" / "sider-ai"
CONFIG_FILE = CONFIG_DIR / "config.json"

# 默认模型
DEFAULT_MODEL = "gpt-4o-mini"

# 支持的模型列表
MODELS = {
    # OpenAI
    "gpt-4o": "OpenAI GPT-4o",
    "gpt-4o-mini": "OpenAI GPT-4o Mini",
    "gpt-4-turbo": "OpenAI GPT-4 Turbo",
    "gpt-3.5-turbo": "OpenAI GPT-3.5 Turbo",
    "o1": "OpenAI o1 推理模型",
    "o1-mini": "OpenAI o1-mini",
    "o1-preview": "OpenAI o1-preview",
    
    # Anthropic
    "claude-3.5-sonnet": "Anthropic Claude 3.5 Sonnet",
    "claude-3.5-haiku": "Anthropic Claude 3.5 Haiku",
    "claude-3-opus": "Anthropic Claude 3 Opus",
    "claude-3-sonnet": "Anthropic Claude 3 Sonnet",
    "claude-3-haiku": "Anthropic Claude 3 Haiku",
    
    # Google
    "gemini-2.0-flash": "Google Gemini 2.0 Flash",
    "gemini-1.5-pro": "Google Gemini 1.5 Pro",
    "gemini-1.5-flash": "Google Gemini 1.5 Flash",
    
    # DeepSeek
    "deepseek-chat": "DeepSeek Chat",
    "deepseek-coder": "DeepSeek Coder",
    
    # Meta
    "llama-3.1-70b": "Meta Llama 3.1 70B",
    "llama-3.1-8b": "Meta Llama 3.1 8B",
}

def load_config():
    """加载配置"""
    if CONFIG_FILE.exists():
        with open(CONFIG_FILE, 'r') as f:
            return json.load(f)
    return {"token": None, "cookie": None}

def save_config(token=None, cookie=None):
    """保存配置"""
    config = load_config()
    if token:
        config["token"] = token
    if cookie:
        config["cookie"] = cookie
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_FILE, 'w') as f:
        json.dump(config, f, indent=2)
    return True

def get_session():
    """获取 sider session"""
    from sider_ai_api import Session
    
    config = load_config()
    token = config.get("token")
    cookie = config.get("cookie")
    
    if not token:
        raise ValueError("请先配置 sider token！使用 set_sider_token(token='xxx', cookie='xxx')")
    
    return Session(token=token, cookie=cookie)

def set_sider_token(token: str, cookie: str = None):
    """
    设置 sider token 和 cookie
    
    Args:
        token: sider.ai 的认证 token
        cookie: 可选的 cookie 字符串
    """
    save_config(token, cookie)
    print(f"✅ Token 已保存到 {CONFIG_FILE}")
    
    # 验证 token
    try:
        session = get_session()
        print(f"✅ Token 验证成功！剩余次数: {session.remain}/{session.total}")
    except Exception as e:
        print(f"⚠️ Token 验证失败: {e}")

def sider_chat(prompt: str, model: str = DEFAULT_MODEL) -> dict:
    """
    与 AI 对话
    
    Args:
        prompt: 用户提示词
        model: 使用的模型 (默认 gpt-4o-mini)
    
    Returns:
        dict: {
            "success": True/False,
            "response": "AI回复内容",
            "model": "使用的模型",
            "remain": 剩余次数,
            "total": 总次数
        }
    """
    try:
        session = get_session()
        
        # 流式获取响应
        response_parts = []
        for part in session.chat(prompt, model):
            response_parts.append(part)
        
        full_response = "".join(response_parts)
        
        return {
            "success": True,
            "response": full_response,
            "model": model,
            "remain": session.remain,
            "total": session.total
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "model": model
        }

def sider_ocr(image_path: str, model: str = "gemini-2.0-flash") -> dict:
    """
    对图像进行 OCR 识别
    
    Args:
        image_path: 图像文件路径
        model: 使用的 OCR 模型 (默认 gemini-2.0-flash)
    
    Returns:
        dict: {
            "success": True/False,
            "response": "识别结果",
            "model": "使用的模型"
        }
    """
    try:
        session = get_session()
        
        # OCR 识别
        response_parts = []
        for part in session.ocr(image_path, model):
            response_parts.append(part)
        
        full_response = "".join(response_parts)
        
        return {
            "success": True,
            "response": full_response,
            "model": model
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "model": model
        }

def list_models():
    """列出所有支持的模型"""
    print("\n📋 支持的模型列表:\n")
    print(f"{'模型ID':<25} {'供应商':<40}")
    print("-" * 65)
    for model_id, desc in MODELS.items():
        print(f"{model_id:<25} {desc:<40}")
    print()

def get_status():
    """获取 API 状态"""
    try:
        session = get_session()
        return {
            "success": True,
            "remain": session.remain,
            "total": session.total,
            "status": "active"
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "status": "inactive"
        }

# 命令行接口
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="sider-ai 工具")
    parser.add_argument("command", choices=["chat", "ocr", "status", "models", "config"], help="命令")
    parser.add_argument("--prompt", "-p", help="对话内容")
    parser.add_argument("--model", "-m", default=DEFAULT_MODEL, help="模型名称")
    parser.add_argument("--image", "-i", help="图像路径 (OCR用)")
    parser.add_argument("--token", "-t", help="设置 token")
    parser.add_argument("--cookie", "-c", help="设置 cookie")
    
    args = parser.parse_args()
    
    if args.command == "models":
        list_models()
        
    elif args.command == "config":
        if args.token:
            set_sider_token(args.token, args.cookie)
        else:
            config = load_config()
            if config.get("token"):
                print(f"✅ Token 已配置: {config['token'][:20]}...")
                print(f"✅ Cookie 已配置: {'是' if config.get('cookie') else '否'}")
            else:
                print("❌ 未配置 token")
                
    elif args.command == "status":
        status = get_status()
        if status["success"]:
            print(f"✅ API 状态: {status['status']}")
            print(f"📊 剩余次数: {status['remain']}/{status['total']}")
        else:
            print(f"❌ 错误: {status['error']}")
            
    elif args.command == "chat":
        if not args.prompt:
            print("❌ 请输入 --prompt 或 -p")
            sys.exit(1)
        result = sider_chat(args.prompt, args.model)
        if result["success"]:
            print(f"\n🤖 {result['model']}:\n")
            print(result["response"])
            print(f"\n📊 剩余次数: {result['remain']}/{result['total']}")
        else:
            print(f"❌ 错误: {result['error']}")
            
    elif args.command == "ocr":
        if not args.image:
            print("❌ 请输入 --image 或 -i")
            sys.exit(1)
        result = sider_ocr(args.image, args.model)
        if result["success"]:
            print(f"\n🔤 OCR 结果 ({result['model']}):\n")
            print(result["response"])
        else:
            print(f"❌ 错误: {result['error']}")