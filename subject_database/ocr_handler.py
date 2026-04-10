#!/usr/bin/env python3
"""
受试者信息OCR识别模块
支持从图片中提取文本和结构化信息
"""

import os
import json
import re
from datetime import datetime
from pathlib import Path

class OCRHandler:
    def __init__(self, base_dir="/root/.openclaw/workspace/subject_database"):
        self.base_dir = Path(base_dir)
        self.photos_dir = self.base_dir / "photos"
        self.ocr_dir = self.base_dir / "ocr_results"
    
    def process_image(self, image_path, use_external=True):
        """
        处理图片，提取文本信息
        use_external: 是否尝试使用外部OCR服务
        """
        from PIL import Image
        
        # 检查图片
        if not os.path.exists(image_path):
            return {"success": False, "error": "图片不存在"}
        
        # 加载图片
        img = Image.open(image_path)
        
        result = {
            "success": True,
            "image_path": image_path,
            "image_size": img.size,
            "processed_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "raw_text": "",
            "structured_data": {},
            "confidence": 0
        }
        
        # 方法1: 使用EasyOCR (如果可用)
        try:
            import easyocr
            reader = easyocr.Reader(['ch_sim', 'en'])
            raw_text = reader.readtext(image_path, detail=0)
            result["raw_text"] = "\n".join(raw_text)
            result["ocr_method"] = "EasyOCR"
            result["confidence"] = 0.8
        except ImportError:
            result["ocr_method"] = "EasyOCR not available"
        
        # 方法2: 尝试pytesseract
        if not result["raw_text"]:
            try:
                import pytesseract
                result["raw_text"] = pytesseract.image_to_string(image_path, lang='chi_sim+eng')
                result["ocr_method"] = "Tesseract"
                result["confidence"] = 0.7
            except ImportError:
                result["ocr_method"] = "Tesseract not available"
        
        # 解析结构化数据
        result["structured_data"] = self._parse_text(result["raw_text"])
        
        return result
    
    def _parse_text(self, text):
        """从文本中解析结构化信息"""
        data = {}
        
        # 提取年龄
        age_match = re.search(r'年龄[：:]\s*(\d+)', text)
        if not age_match:
            age_match = re.search(r'(\d{1,2})\s*岁', text)
        if age_match:
            data["age"] = int(age_match.group(1))
        
        # 提取性别
        if "男" in text and "女" not in text:
            data["gender"] = "男"
        elif "女" in text:
            data["gender"] = "女"
        gender_match = re.search(r'性别[：:]\s*([男女])', text)
        if gender_match:
            data["gender"] = gender_match.group(1)
        
        # 提取ID
        id_match = re.search(r'编号[：:]\s*([A-Z0-9]+)', text)
        if not id_match:
            id_match = re.search(r'ID[：:]\s*([A-Z0-9]+)', text)
        if id_match:
            data["original_id"] = id_match.group(1)
        
        # 提取PHQ-9
        phq9_match = re.search(r'PHQ-9[：:]\s*(\d+)', text)
        if phq9_match:
            data["phq9_score"] = int(phq9_match.group(1))
        
        # 提取GAD-7
        gad7_match = re.search(r'GAD-7[：:]\s*(\d+)', text)
        if gad7_match:
            data["gad7_score"] = int(gad7_match.group(1))
        
        # 提取组别
        if "患者" in text or "病例" in text:
            data["group"] = "患者"
        elif "对照" in text or "正常" in text or "健康" in text:
            data["group"] = "对照"
        
        # 提取眼轴
        axial_od = re.search(r'眼轴.*?右.*?(\d+\.?\d*)', text)
        axial_os = re.search(r'眼轴.*?左.*?(\d+\.?\d*)', text)
        if axial_od:
            data["axial_length_od"] = float(axial_od.group(1))
        if axial_os:
            data["axial_length_os"] = float(axial_os.group(1))
        
        # 提取日期
        date_match = re.search(r'(\d{4})[年/-](\d{1,2})[月/-](\d{1,2})', text)
        if date_match:
            data["enrollment_date"] = f"{date_match.group(1)}-{int(date_match.group(2)):02d}-{int(date_match.group(3)):02d}"
        
        return data
    
    def process_folder(self, folder_path):
        """批量处理文件夹中的图片"""
        results = []
        folder = Path(folder_path)
        
        for img_file in folder.glob("*.[jp][pn][g]"):
            result = self.process_image(str(img_file))
            result["filename"] = img_file.name
            results.append(result)
        
        return results
    
    def generate_report(self, ocr_result):
        """生成OCR识别报告"""
        data = ocr_result.get("structured_data", {})
        
        report = f"""
╔══════════════════════════════════════════════════════════════╗
║                    OCR识别结果报告                            ║
╠══════════════════════════════════════════════════════════════╣
║ 图像: {ocr_result.get('image_path', 'N/A'):<45} ║
║ 方法: {ocr_result.get('ocr_method', 'N/A'):<45} ║
║ 置信度: {ocr_result.get('confidence', 0):.0%}{'':>38} ║
╠══════════════════════════════════════════════════════════════╣
║ 识别到的信息:                                                 ║
║   - 年龄: {data.get('age', '未识别'):<40} ║
║   - 性别: {data.get('gender', '未识别'):<40} ║
║   - 组别: {data.get('group', '未识别'):<40} ║
║   - PHQ-9: {data.get('phq9_score', '未识别'):<40} ║
║   - GAD-7: {data.get('gad7_score', '未识别'):<40} ║
║   - 眼轴(右): {data.get('axial_length_od', '未识别'):<37} ║
║   - 眼轴(左): {data.get('axial_length_os', '37} ║
╚══════════════════════════════════════════════════════════════╝
        """
        return report


def main():
    import sys
    
    if len(sys.argv) < 2:
        print("""
╔══════════════════════════════════════════════════════════════╗
║           受试者信息OCR识别工具 v1.0                          ║
╠══════════════════════════════════════════════════════════════╣
║ 使用方法:                                                     ║
║   python3 ocr_handler.py <图片路径>                          ║
║   python3 ocr_handler.py --folder <文件夹>                   ║
║   python3 ocr_handler.py --interactive                       ║
╚══════════════════════════════════════════════════════════════╝
        """)
        return
    
    handler = OCRHandler()
    
    if sys.argv[1] == "--folder":
        folder = sys.argv[2] if len(sys.argv) > 2 else "."
        results = handler.process_folder(folder)
        print(f"处理了 {len(results)} 张图片")
        
    elif sys.argv[1] == "--interactive":
        print("请输入图片路径:")
        path = input().strip()
        result = handler.process_image(path)
        print(handler.generate_report(result))
        
    else:
        path = sys.argv[1]
        result = handler.process_image(path)
        print(handler.generate_report(result))


if __name__ == "__main__":
    main()
