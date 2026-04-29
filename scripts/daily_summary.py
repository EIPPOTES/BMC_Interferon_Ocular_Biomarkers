#!/usr/bin/env python3
"""
每日会话自动摘要 - 自动提取关键信息并写入memory
"""

import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

SESSIONS_BASE = Path.home() / ".openclaw" / "agents"
MEMORY_DIR = Path.home() / ".openclaw" / "workspace" / "memory"


def get_today_sessions():
    """获取今天所有agent的session文件"""
    today = datetime.now().strftime("%Y-%m-%d")
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
            # 检查文件修改时间
            mtime = datetime.fromtimestamp(session_file.stat().st_mtime)
            if mtime.strftime("%Y-%m-%d") == today:
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


def parse_session(session_file):
    """解析session文件，提取关键信息"""
    messages = []
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
                    messages.append({
                        "role": role,
                        "text": text[:500],  # 限制长度
                        "timestamp": ts
                    })
    except Exception as e:
        print(f"Error reading {session_file}: {e}")
    
    return messages


def extract_key_info(messages):
    """从消息中提取关键信息"""
    tasks = []
    decisions = []
    results = []
    
    user_keywords = ["需要", "帮我", "请", "要", "做", "完成", "处理", "开始", "执行"]
    decision_keywords = ["决定", "选择", "采用", "使用", "确定", "批准"]
    result_keywords = ["完成", "成功", "已", "好了", "可以", "结束"]
    
    for msg in messages:
        text = msg["text"]
        # 清理元数据
        lines = text.split("\n")
        clean_lines = []
        for line in lines:
            line = line.strip()
            # 移除系统元数据
            if any(line.startswith(k) for k in ["```", "Conversation", "Sender", "{"]) and not line.startswith("```json"):
                continue
            if line and len(line) > 2:
                clean_lines.append(line)
        
        clean_text = " ".join(clean_lines[:3])  # 取前3行
        
        if not clean_text or len(clean_text) < 5:
            continue
            
        # 检测任务
        if any(kw in clean_text for kw in user_keywords):
            tasks.append(clean_text[:150])
        # 检测决策
        if any(kw in clean_text for kw in decision_keywords):
            decisions.append(clean_text[:150])
        # 检测结果
        if any(kw in clean_text for kw in result_keywords):
            results.append(clean_text[:150])
    
    # 去重
    tasks = list(dict.fromkeys(tasks))[:5]
    decisions = list(dict.fromkeys(decisions))[:3]
    results = list(dict.fromkeys(results))[:5]
    
    return tasks, decisions, results


def generate_summary():
    """生成每日摘要"""
    sessions = get_today_sessions()
    
    if not sessions:
        return None
    
    all_tasks = []
    all_decisions = []
    all_results = []
    agent_stats = {}
    
    for session_file in sessions:
        agent_name = session_file.parent.parent.name
        messages = parse_session(session_file)
        
        if agent_name not in agent_stats:
            agent_stats[agent_name] = 0
        agent_stats[agent_name] += len(messages)
        
        tasks, decisions, results = extract_key_info(messages)
        all_tasks.extend(tasks)
        all_decisions.extend(decisions)
        all_results.extend(results)
    
    # 去重
    all_tasks = list(dict.fromkeys(all_tasks))[:8]
    all_decisions = list(dict.fromkeys(all_decisions))[:5]
    all_results = list(dict.fromkeys(all_results))[:8]
    
    today = datetime.now().strftime("%Y-%m-%d")
    
    summary = f"""# 每日会话摘要 - {today}

## 概览
- 会话数: {len(sessions)}
- 总消息: {sum(agent_stats.values())}
- Agent: {", ".join(agent_stats.keys())}

## 关键任务
"""
    for i, task in enumerate(all_tasks, 1):
        summary += f"{i}. {task}\n"
    
    summary += "\n## 重要决策\n"
    for i, decision in enumerate(all_decisions, 1):
        summary += f"{i}. {decision}\n"
    
    summary += "\n## 完成结果\n"
    for i, result in enumerate(all_results, 1):
        summary += f"{i}. {result}\n"
    
    summary += f"\n---\n*自动生成于 {datetime.now().strftime('%H:%M:%S')}*"
    
    return summary


def main():
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] 开始每日会话摘要...")
    
    summary = generate_summary()
    
    if summary:
        today = datetime.now().strftime("%Y-%m-%d")
        output_file = MEMORY_DIR / f"{today}-auto-summary.md"
        
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(summary)
        
        print(f"✅ 摘要已写入: {output_file}")
    else:
        print("📝 今日无会话记录")


if __name__ == "__main__":
    main()