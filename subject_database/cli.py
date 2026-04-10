#!/usr/bin/env python3
"""
受试者数据库命令行工具
简单易用的交互界面
"""

import os
import sys
import json
from pathlib import Path

# 添加当前目录到路径
sys.path.insert(0, str(Path(__file__).parent))

from database_manager import SubjectDatabase
from ocr_handler import OCRHandler

def print_header():
    print("\n" + "=" * 60)
    print("      👁️ 眼科研究受试者数据库管理系统")
    print("=" * 60)

def print_menu():
    print("""
┌─────────────────────────────────────────────────────────────┐
│ 菜单选项:                                                    │
│  [1] 📸 拍照识别 - 从图片提取受试者信息                     │
│  [2] 📋 手动录入 - 手动输入受试者信息                       │
│  [3] 🔍 查询    - 查看已有受试者                            │
│  [4] 📊 统计    - 查看统计信息                              │
│  [5] ✏️  修改   - 修改受试者信息                             │
│  [6] 🗑️  删除   - 删除受试者                                 │
│  [7] 💾 导出    - 导出为CSV文件                              │
│  [0] 🚪 退出    - 退出系统                                   │
└─────────────────────────────────────────────────────────────┘
""")

def option_1_ocr(db, handler):
    """拍照识别"""
    print("\n📸 拍照识别模式")
    print("-" * 40)
    path = input("请输入图片路径 (或拖入图片): ").strip().strip('"\'')
    
    if not os.path.exists(path):
        print("❌ 文件不存在！")
        return
    
    print("🔄 正在识别...")
    result = handler.process_image(path)
    
    if result["success"]:
        print(f"\n✅ 识别完成！")
        print(f"   方法: {result.get('ocr_method', 'N/A')}")
        print(f"   置信度: {result.get('confidence', 0):.0%}")
        
        data = result.get("structured_data", {})
        print("\n识别到的信息:")
        for k, v in data.items():
            print(f"   {k}: {v}")
        
        # 保存
        save = input("\n是否保存到数据库? (y/n): ").strip().lower()
        if save == 'y':
            # 分配ID
            max_id = 0
            for s in db.db["subjects"]:
                if s["subject_id"].startswith("S"):
                    try:
                        max_id = max(max_id, int(s["subject_id"][1:]))
                    except:
                        pass
            new_id = f"S{max_id + 1:04d}"
            
            # 添加照片
            import shutil
            ext = os.path.splitext(path)[1]
            new_photo = db.photos_dir / f"{new_id}{ext}"
            shutil.copy2(path, new_photo)
            
            # 保存到数据库
            db.add_subject({
                "subject_id": new_id,
                "status": "已录入",
                "original_id": data.get("original_id", ""),
                "age": data.get("age", 0),
                "gender": data.get("gender", ""),
                "group": data.get("group", ""),
                "phq9_score": data.get("phq9_score", 0),
                "gad7_score": data.get("gad7_score", 0),
                "axial_length_od": data.get("axial_length_od", 0),
                "axial_length_os": data.get("axial_length_os", 0),
                "enrollment_date": data.get("enrollment_date", ""),
                "photo_path": str(new_photo)
            }, None)
            
            # 保存OCR结果
            handler.save_ocr_result(new_id, result)
            print(f"✅ 已保存！受试者ID: {new_id}")
    else:
        print(f"❌ 识别失败: {result.get('error', '未知错误')}")

def option_2_manual(db):
    """手动录入"""
    print("\n📋 手动录入模式")
    print("-" * 40)
    
    data = {}
    data["age"] = int(input("年龄: ").strip())
    data["gender"] = input("性别 (男/女): ").strip()
    data["group"] = input("组别 (患者/对照): ").strip()
    
    opt = input("是否记录PHQ-9? (y/n): ").strip().lower()
    if opt == 'y':
        data["phq9_score"] = float(input("PHQ-9评分: ").strip())
    
    opt = input("是否记录GAD-7? (y/n): ").strip().lower()
    if opt == 'y':
        data["gad7_score"] = float(input("GAD-7评分: ").strip())
    
    opt = input("是否记录眼轴? (y/n): ").strip().lower()
    if opt == 'y':
        data["axial_length_od"] = float(input("眼轴(右): ").strip())
        data["axial_length_os"] = float(input("眼轴(左): ").strip())
    
    data["notes"] = input("备注: ").strip()
    data["status"] = "已录入"
    
    subject_id = db.add_subject(data)
    print(f"✅ 已保存！受试者ID: {subject_id}")

def option_3_query(db):
    """查询"""
    print("\n🔍 查询模式")
    print("-" * 40)
    print("1. 按ID查询")
    print("2. 按状态查询")
    print("3. 按组别查询")
    print("4. 显示所有")
    
    choice = input("选择: ").strip()
    
    if choice == "1":
        sid = input("受试者ID: ").strip()
        subjects = [s for s in db.db["subjects"] if s["subject_id"] == sid]
    elif choice == "2":
        status = input("状态 (待处理/已录入/已核实/已排除): ").strip()
        subjects = db.get_subjects(status=status)
    elif choice == "3":
        group = input("组别 (患者/对照): ").strip()
        subjects = db.get_subjects(group=group)
    else:
        subjects = db.get_subjects()
    
    print(f"\n找到 {len(subjects)} 条记录:")
    for s in subjects[:20]:
        print(f"  [{s['subject_id']}] {s.get('gender','')} {s.get('age','')}岁 {s.get('group','')} | {s.get('status','')}")
    
    if len(subjects) > 20:
        print(f"  ... 还有 {len(subjects)-20} 条")

def option_4_stats(db):
    """统计"""
    print("\n📊 统计信息")
    print("-" * 40)
    stats = db.get_stats()
    for k, v in stats.items():
        print(f"  {k}: {v}")

def main():
    db = SubjectDatabase()
    handler = OCRHandler()
    
    while True:
        print_header()
        print_menu()
        option_4_stats(db)
        
        choice = input("\n请选择操作: ").strip()
        
        if choice == "0":
            print("\n👋 再见！")
            break
        elif choice == "1":
            option_1_ocr(db, handler)
        elif choice == "2":
            option_2_manual(db)
        elif choice == "3":
            option_3_query(db)
        elif choice == "4":
            option_4_stats(db)
        elif choice == "5":
            sid = input("受试者ID: ").strip()
            field = input("要修改的字段: ").strip()
            value = input("新值: ").strip()
            if db.update_subject(sid, {field: value}):
                print("✅ 更新成功")
            else:
                print("❌ 更新失败")
        elif choice == "6":
            sid = input("要删除的受试者ID: ").strip()
            confirm = input(f"确认删除 {sid}? (y/n): ").strip().lower()
            if confirm == 'y' and db.delete_subject(sid):
                print("✅ 删除成功")
            else:
                print("❌ 删除取消")
        elif choice == "7":
            filename = input("导出文件名: ").strip()
            if db.export_csv(filename):
                print(f"✅ 已导出到 {filename}")
            else:
                print("❌ 导出失败或无数据")
        
        input("\n按回车继续...")

if __name__ == "__main__":
    main()
