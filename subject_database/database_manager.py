#!/usr/bin/env python3
"""
受试者信息数据库管理系统
功能：拍照识别、OCR提取、数据管理
"""

import os
import json
import shutil
from datetime import datetime
from pathlib import Path

DATABASE_FILE = "database.json"
PHOTOS_DIR = "photos"
OCR_DIR = "ocr_results"

class SubjectDatabase:
    def __init__(self, base_dir="/root/.openclaw/workspace/subject_database"):
        self.base_dir = Path(base_dir)
        self.db_file = self.base_dir / DATABASE_FILE
        self.photos_dir = self.base_dir / PHOTOS_DIR
        self.ocr_dir = self.base_dir / OCR_DIR
        
        # 确保目录存在
        self.photos_dir.mkdir(exist_ok=True)
        self.ocr_dir.mkdir(exist_ok=True)
        
        # 加载或初始化数据库
        if self.db_file.exists():
            with open(self.db_file, 'r', encoding='utf-8') as f:
                self.db = json.load(f)
        else:
            self._init_database()
    
    def _init_database(self):
        self.db = {
            "last_update": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "total_subjects": 0,
            "patients": 0,
            "controls": 0,
            "subjects": []
        }
        self._save()
    
    def _save(self):
        self.db["last_update"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(self.db_file, 'w', encoding='utf-8') as f:
            json.dump(self.db, f, ensure_ascii=False, indent=2)
    
    def add_subject(self, subject_data, photo_path=None):
        """添加新受试者"""
        # 生成ID
        max_id = 0
        for s in self.db["subjects"]:
            if s["subject_id"].startswith("S"):
                try:
                    max_id = max(max_id, int(s["subject_id"][1:]))
                except:
                    pass
        new_id = f"S{max_id + 1:04d}"
        
        subject = {
            "subject_id": new_id,
            "status": "待处理",
            "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "photo_path": None,
            "ocr_data": {}
        }
        
        # 合并数据
        subject.update(subject_data)
        
        # 复制照片
        if photo_path and os.path.exists(photo_path):
            ext = os.path.splitext(photo_path)[1]
            new_photo_name = f"{new_id}{ext}"
            new_photo_path = self.photos_dir / new_photo_name
            shutil.copy2(photo_path, new_photo_path)
            subject["photo_path"] = str(new_photo_path)
        
        # 更新统计
        self.db["subjects"].append(subject)
        self._update_stats()
        self._save()
        
        return new_id
    
    def _update_stats(self):
        self.db["total_subjects"] = len(self.db["subjects"])
        self.db["patients"] = sum(1 for s in self.db["subjects"] if s.get("group") == "患者")
        self.db["controls"] = sum(1 for s in self.db["subjects"] if s.get("group") == "对照")
    
    def get_subjects(self, status=None, group=None):
        """查询受试者"""
        results = self.db["subjects"]
        if status:
            results = [s for s in results if s.get("status") == status]
        if group:
            results = [s for s in results if s.get("group") == group]
        return results
    
    def update_subject(self, subject_id, data):
        """更新受试者信息"""
        for s in self.db["subjects"]:
            if s["subject_id"] == subject_id:
                s.update(data)
                self._save()
                return True
        return False
    
    def delete_subject(self, subject_id):
        """删除受试者"""
        for i, s in enumerate(self.db["subjects"]):
            if s["subject_id"] == subject_id:
                # 删除照片
                if s.get("photo_path") and os.path.exists(s["photo_path"]):
                    os.remove(s["photo_path"])
                del self.db["subjects"][i]
                self._update_stats()
                self._save()
                return True
        return False
    
    def save_ocr_result(self, subject_id, ocr_data):
        """保存OCR结果"""
        ocr_file = self.ocr_dir / f"{subject_id}_ocr.json"
        with open(ocr_file, 'w', encoding='utf-8') as f:
            json.dump(ocr_data, f, ensure_ascii=False, indent=2)
        
        # 更新数据库
        self.update_subject(subject_id, {
            "ocr_data": ocr_data,
            "status": "已录入"
        })
        return str(ocr_file)
    
    def get_stats(self):
        """获取统计信息"""
        return {
            "总数": self.db["total_subjects"],
            "患者": self.db["patients"],
            "对照": self.db["controls"],
            "待处理": len(self.get_subjects("待处理")),
            "已录入": len(self.get_subjects("已录入")),
            "已核实": len(self.get_subjects("已核实")),
            "已排除": len(self.get_subjects("已排除"))
        }
    
    def export_csv(self, filename):
        """导出CSV"""
        import csv
        if not self.db["subjects"]:
            return False
        
        fields = list(self.db["subjects"][0].keys())
        with open(filename, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(self.db["subjects"])
        return True


if __name__ == "__main__":
    db = SubjectDatabase()
    print("受试者数据库管理系统")
    print(f"当前统计: {db.get_stats()}")
