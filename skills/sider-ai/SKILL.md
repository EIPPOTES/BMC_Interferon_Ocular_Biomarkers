# sider-ai Skill

通过 sider.ai API 调用 GPT、Claude、Gemini、DeepSeek、o1 等全部模型。

## 功能

- **chat()** - 与 AI 对话（流式输出）
- **ocr()** - 图像文字识别
- **多模型支持** - 全部 sider.ai 支持的模型

## 依赖

```bash
pip install sider-ai-api
```

## 模型列表

| 模型ID | 供应商 |
|--------|--------|
| `gpt-4o-mini` | OpenAI |
| `gpt-4o` | OpenAI |
| `gpt-4-turbo` | OpenAI |
| `gpt-3.5-turbo` | OpenAI |
| `claude-3.5-sonnet` | Anthropic |
| `claude-3.5-haiku` | Anthropic |
| `claude-3-opus` | Anthropic |
| `gemini-2.0-flash` | Google |
| `gemini-1.5-pro` | Google |
| `deepseek-chat` | DeepSeek |
| `deepseek-coder` | DeepSeek |
| `o1` | OpenAI推理 |
| `o1-mini` | OpenAI推理 |
| `llama-3.1-70b` | Meta |

## 使用方式

### 1. 对话调用

```python
import json
import sys
sys.path.insert(0, '/root/.openclaw/workspace/skills/sider-ai')

from sider import sider_chat

# 调用 GPT-4o
result = sider_chat("解释什么是OCT血管成像", model="gpt-4o-mini")
print(result)

# 调用 Claude
result = sider_chat("解释抑郁症与视网膜变薄的关联", model="claude-3.5-sonnet")
print(result)

# 调用 Gemini
result = sider_chat("解释青光眼的发病机制", model="gemini-2.0-flash")
print(result)

# 调用 o1 推理模型
result = sider_chat("分析这个研究设计的统计学效力", model="o1")
print(result)
```

### 2. OCR 识别

```python
from sider import sider_ocr

# 图像文字识别
result = sider_ocr("/path/to/image.jpg", model="gemini-2.0-flash")
print(result)
```

## 配置 Token

首次使用需配置 token：

```python
from sider import set_sider_token

# 设置 token 和 cookie
set_sider_token(token="your_token_here", cookie="your_cookie_here")

# 保存到配置文件
save_config()
```

Token 获取方式：
1. 登录 sider.ai
2. 浏览器 F12 → Application → Cookies → sider.ai
3. 或访问 `edge://settings/cookies/detail?site=sider.ai`

## 返回格式

```python
{
    "success": True,
    "response": "AI回复内容",
    "model": "gpt-4o-mini",
    "remain": 999,
    "total": 1000
}
```

## 错误处理

| 错误码 | 说明 |
|--------|------|
| 401 | Token 无效，请检查配置 |
| 403 | Cookie 无效或过期 |
| 429 | 请求过于频繁，请稍后重试 |
| 500 | 服务端错误，请重试 |

---

*Created: 2026-04-01*
*Version: 1.0*