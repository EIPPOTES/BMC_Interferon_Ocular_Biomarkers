#!/usr/bin/env python3
"""
用户偏好自动学习系统 - 从历史会话中提取用户习惯
"""

import json
import os
import re
from collections import Counter
from datetime import datetime, timedelta, timezone
from pathlib import Path

SESSIONS_BASE = Path.home() / ".openclaw" / "agents"
MEMORY_DIR = Path.home() / ".openclaw" / "workspace" / "memory"


def get_recent_sessions(days=7):
    """获取最近N天的session文件"""
    cutoff = (datetime.now() - timedelta(days=days)).timestamp()
    results = []
    
    if not SESSIONS_BASE.is_dir():
        return results
    
    for agent_dir in SESSIONS_BASE.iterdir():
        if not agent_dir.is_dir():
            continue
        sessions_dir = agent_dir / "sessions"
        if not sessions_dir.is_dir():
            continue
            
        for session_file in sessions_dir.glob("*.jsonl"):
            if session_file.stat().st_mtime >= cutoff:
                results.append(session_file)
    
    return results


def extract_text_from_content(content) -> str:
    """提取消息内容"""
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts = []
        for item in content:
            if isinstance(item, dict) and item.get("type") == "text":
                parts.append(item.get("text", ""))
        return " ".join(parts)
    return ""


def extract_preferences(messages):
    """从消息中提取用户偏好"""
    # 收集用户消息
    user_messages = []
    assistant_messages = []
    
    for msg in messages:
        if msg["role"] == "user":
            user_messages.append(msg["text"])
        elif msg["role"] == "assistant":
            assistant_messages.append(msg["text"])
    
    preferences = {
        "communication_style": [],
        "research_interests": [],
        "time_patterns": [],
        "preferred_tools": [],
        "language": "zh-CN",  # 默认中文
    }
    
    # 分析沟通风格
    long_responses = sum(1 for m in user_messages if len(m) > 200)
    short_responses = sum(1 for m in user_messages if len(m) < 50)
    total = len(user_messages) or 1
    
    if long_responses / total > 0.5:
        preferences["communication_style"].append("详细")
    if short_responses / total > 0.5:
        preferences["communication_style"].append("简洁")
    
    # 分析研究兴趣关键词
    research_keywords = [
        "眼科", "OCT", "视网膜", "抑郁症", "MDD", "精神",
        "论文", "投稿", "统计", "分析", "数据",
        "青光眼", "AMD", "糖尿病", "黄斑"
    ]
    
    all_text = " ".join(user_messages + assistant_messages)
    interests = [kw for kw in research_keywords if kw in all_text]
    preferences["research_interests"] = list(set(interests))[:10]
    
    # 分析时间模式
    timestamps = []
    for msg in messages:
        ts = msg.get("timestamp")
        if ts:
            try:
                if isinstance(ts, (int, float)):
                    ts = ts / 1000 if ts > 1e12 else ts
                    dt = datetime.fromtimestamp(ts, tz=timezone.utc)
                    timestamps.append(dt.hour)
            except:
                pass
    
    if timestamps:
        hour_counts = Counter(timestamps)
        peak_hours = hour_counts.most_common(3)
        preferences["time_patterns"] = [f"{h}:00" for h, _ in peak_hours]
    
    # 分析常用工具
    tool_keywords = [
        ("python", "Python"),
        ("R ", "R语言"),
        ("stat", "统计分析"),
        ("web", "网络搜索"),
        ("pubmed", "文献检索"),
        ("cron", "定时任务"),
    ]
    
    tools = []
    for kw, name in tool_keywords:
        if kw.lower() in all_text.lower():
            tools.append(name)
    preferences["preferred_tools"] = list(set(tools))[:5]
    
    return preferences


def generate_preference_report(preferences):
    """生成偏好报告"""
    report = f"""# 用户偏好报告 - {datetime.now().strftime('%Y-%m-%d')}

## 沟通风格
{", ".join(preferences["communication_style"]) or "待分析"}

## 研究兴趣
{", ".join(preferences["research_interests"]) or "待分析"}

## 时间模式 (活跃时段)
{", ".join(preferences["time_patterns"]) or "待分析"}

## 常用工具
{", ".join(preferences["preferred_tools"]) or "待分析"}

## 语言
{preferences["language"]}

---
*基于最近7天会话分析，自动生成于 {datetime.now().strftime('%H:%M:%S')}*
"""
    return report


def main():
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] 开始分析用户偏好...")
    
    sessions = get_recent_sessions(days=7)
    print(f"找到 {len(sessions)} 个最近会话")
    
    all_messages = []
    for session_file in sessions:
        try:
            with open(session_file, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    if rec.get("type") != "message":
                        continue
                    msg = rec.get("message", {})
                    role = msg.get("role", "")
                    text = extract_text_from_content(msg.get("content", ""))
                    if text:
                        ts = rec.get("timestamp") or msg.get("timestamp")
                        all_messages.append({
                            "role": role,
                            "text": text,
                            "timestamp": ts
                        })
        except Exception as e:
            print(f"Error reading {session_file}: {e}")
    
    print(f"分析 {len(all_messages)} 条消息...")
    
    preferences = extract_preferences(all_messages)
    
    # 保存偏好报告
    output_file = MEMORY_DIR / "user-preferences.md"
    report = generate_preference_report(preferences)
    
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(report)
    
    print(f"✅ 偏好分析已写入: {output_file}")
    print(f"📊 研究兴趣: {', '.join(preferences['research_interests'])}")
    print(f"🕐 活跃时段: {', '.join(preferences['time_patterns'])}")


if __name__ == "__main__":
    main()