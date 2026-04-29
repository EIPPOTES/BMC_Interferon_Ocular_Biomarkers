#!/usr/bin/env python3
"""
项目状态感知系统 - 自动追踪研究项目进度
"""

import json
import os
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path

SESSIONS_BASE = Path.home() / ".openclaw" / "agents"
WORKSPACE = Path.home() / ".openclaw" / "workspace"
TASK_FILE = WORKSPACE / "TASK_TRACKER.md"


# 项目关键词映射
PROJECT_KEYWORDS = {
    "OCT-MDD论文": ["OCT-MDD", "视网膜", "抑郁症", "MDD", "论文投稿", "JAD", "Journal of Affective"],
    "年轻组论文": ["年轻组", "v35", "年轻患者"],
    "知识库构建": ["知识库", "眼科", "OpenAlex", "文献检索"],
    "系统优化": ["OpenClaw", "记忆系统", "自动化", "优化"],
}

# 状态关键词
STATUS_KEYWORDS = {
    "已完成": ["完成", "成功", "好了", "已提交", "已投稿", "done", "finished", "submitted"],
    "进行中": ["进行", "开始", "正在", "执行", "running", "working", "processing"],
    "待处理": ["待", "需要", "未完成", "pending", "todo", "还没"],
    "有问题": ["错误", "失败", "问题", "error", "fail", "issue"],
}


def get_recent_sessions(days=3):
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


def detect_projects(messages):
    """检测消息中的项目"""
    project_mentions = {name: [] for name in PROJECT_KEYWORDS}
    
    for msg in messages:
        text = msg.get("text", "")
        # 清理元数据
        lines = text.split("\n")
        clean_text = " ".join([l.strip() for l in lines[:5] if l.strip() and not l.startswith("```")])
        
        for proj_name, keywords in PROJECT_KEYWORDS.items():
            for kw in keywords:
                if kw in clean_text:
                    project_mentions[proj_name].append(clean_text[:200])
                    break
    
    return {k: v for k, v in project_mentions.items() if v}


def detect_status(text):
    """检测文本中的状态"""
    for status, keywords in STATUS_KEYWORDS.items():
        for kw in keywords:
            if kw in text:
                return status
    return "进行中"


def analyze_project_updates(messages):
    """分析项目状态更新"""
    project_updates = {}
    
    for msg in messages:
        text = msg.get("text", "")
        
        # 检测项目
        for proj_name, keywords in PROJECT_KEYWORDS.items():
            if any(kw in text for kw in keywords):
                status = detect_status(text)
                if proj_name not in project_updates:
                    project_updates[proj_name] = {
                        "status": status,
                        "mentions": [],
                        "latest": text[:300]
                    }
                project_updates[proj_name]["mentions"].append(text[:200])
                
                # 更新为更高优先级状态
                status_priority = {"已完成": 4, "进行中": 3, "待处理": 2, "有问题": 1}
                current_priority = status_priority.get(project_updates[proj_name]["status"], 0)
                new_priority = status_priority.get(status, 0)
                if new_priority >= current_priority:
                    project_updates[proj_name]["status"] = status
    
    return project_updates


def generate_project_report(project_updates):
    """生成项目状态报告"""
    if not project_updates:
        return None
    
    report = f"""## 🔄 近期项目状态 (自动检测) - {datetime.now().strftime('%Y-%m-%d')}

| 项目 | 检测状态 | 活动数 |
|------|----------|--------|
"""
    for proj, data in project_updates.items():
        count = len(data["mentions"])
        status = data["status"]
        report += f"| {proj} | {status} | {count} |\n"
    
    report += "\n### 详细活动\n"
    for proj, data in list(project_updates.items())[:5]:
        report += f"\n**{proj}** ({data['status']}):\n"
        for mention in data["mentions"][:2]:
            report += f"- {mention[:100]}...\n"
    
    return report


def update_task_tracker(project_updates):
    """更新TASK_TRACKER.md"""
    if not project_updates:
        return False
    
    report = generate_project_report(project_updates)
    if not report:
        return False
    
    today = datetime.now().strftime('%Y-%m-%d %H:%M')
    
    # 追加到现有文件
    try:
        existing = TASK_FILE.read_text(encoding="utf-8")
    except:
        existing = "# TASK_TRACKER.md\n"
    
    # 检查是否已有今日记录
    if f"### {datetime.now().strftime('%Y-%m-%d')}" in existing:
        # 更新现有记录
        lines = existing.split("\n")
        new_lines = []
        for i, line in enumerate(lines):
            new_lines.append(line)
            if f"### {datetime.now().strftime('%Y-%m-%d')}" in line:
                new_lines.append("")
                new_lines.append(report)
                break
        content = "\n".join(new_lines)
    else:
        # 追加新记录
        content = existing + f"\n\n---\n\n{report}\n"
    
    TASK_FILE.write_text(content, encoding="utf-8")
    return True


def main():
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] 开始项目状态分析...")
    
    sessions = get_recent_sessions(days=3)
    print(f"分析 {len(sessions)} 个会话...")
    
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
            print(f"Error: {e}")
    
    print(f"处理 {len(all_messages)} 条消息...")
    
    # 分析项目更新
    project_updates = analyze_project_updates(all_messages)
    
    if project_updates:
        # 更新TASK_TRACKER
        if update_task_tracker(project_updates):
            print(f"✅ 已更新TASK_TRACKER.md")
            for proj, data in project_updates.items():
                print(f"  📌 {proj}: {data['status']}")
        else:
            print("⚠️ 未能更新TASK_TRACKER")
    else:
        print("📝 未检测到项目状态变化")


if __name__ == "__main__":
    main()