from PIL import Image, ImageStat
import os
import numpy as np


def main():
    """主函数，包装原有执行代码"""
    def analyze_image(img_path):
        """分析图片，检测可能的乱码区域"""
        results = {
            'filename': os.path.basename(img_path),
            'size': None,
            'mode': None,
            'has_transparency': False,
            'color_distribution': {},
            'suspicious_areas': [],
            'text_like_regions': 0
        }

        try:
            img = Image.open(img_path)
            results['size'] = img.size
            results['mode'] = img.mode

            # 转换为RGB分析
            if img.mode in ('RGBA', 'P'):
                if img.mode == 'RGBA':
                    # 检查透明度
                    alpha = img.split()[-1]
                    results['has_transparency'] = any(p < 255 for p in alpha.getdata())
                img_rgb = img.convert('RGB')
            else:
                img_rgb = img.convert('RGB') if img.mode != 'RGB' else img

            # 分析颜色分布
            stat = ImageStat.Stat(img_rgb)
            results['color_distribution'] = {
                'mean_r': stat.mean[0],
                'mean_g': stat.mean[1],
                'mean_b': stat.mean[2],
                'stddev': stat.stddev
            }

            # 转换为灰度分析对比度
            img_gray = img_rgb.convert('L')
            gray_stat = ImageStat.Stat(img_gray)

            # 检测高对比度区域（可能是文字）
            # 如果标准差很低，可能是纯色块（乱码方块）
            if gray_stat.stddev[0] < 10:
                results['suspicious_areas'].append(f"低对比度区域 (std={gray_stat.stddev[0]:.1f}) - 可能是纯色块")

            # 分析像素分布
            pixels = np.array(img_gray)
            white_ratio = np.sum(pixels > 240) / pixels.size
            black_ratio = np.sum(pixels < 20) / pixels.size
            gray_ratio = np.sum((pixels >= 20) & (pixels <= 240)) / pixels.size

            results['pixel_distribution'] = {
                'white_ratio': white_ratio,
                'black_ratio': black_ratio,
                'gray_ratio': gray_ratio
            }

            # 检测可能的乱码特征
            # 乱码方块通常有特定的颜色模式
            if gray_ratio > 0.3 and gray_stat.stddev[0] > 30:
                results['text_like_regions'] = int(gray_ratio * 100)

            # 如果黑色像素比例异常高，可能有黑色方块（乱码）
            if black_ratio > 0.05:
                results['suspicious_areas'].append(f"高黑色像素比例 ({black_ratio*100:.1f}%) - 可能有黑色方块")

            return results

        except Exception as e:
            return {'filename': os.path.basename(img_path), 'error': str(e)}

    def print_analysis(results):
        """打印分析结果"""
        print(f"\n{'='*60}")
        print(f"文件: {results['filename']}")
        print(f"{'='*60}")

        if 'error' in results:
            print(f"错误: {results['error']}")
            return

        print(f"尺寸: {results['size'][0]} x {results['size'][1]} 像素")
        print(f"模式: {results['mode']}")
        print(f"透明度: {'是' if results['has_transparency'] else '否'}")

        print(f"\n颜色分布:")
        cd = results['color_distribution']
        print(f"  平均RGB: ({cd['mean_r']:.1f}, {cd['mean_g']:.1f}, {cd['mean_b']:.1f})")
        print(f"  标准差: {cd['stddev']}")

        print(f"\n像素分布:")
        pd = results['pixel_distribution']
        print(f"  白色: {pd['white_ratio']*100:.1f}%")
        print(f"  黑色: {pd['black_ratio']*100:.1f}%")
        print(f"  灰色(文字): {pd['gray_ratio']*100:.1f}%")

        print(f"\n文字区域评估: {results['text_like_regions']}%")

        if results['suspicious_areas']:
            print(f"\n⚠️ 可疑区域:")
            for area in results['suspicious_areas']:
                print(f"  - {area}")
        else:
            print(f"\n✅ 未检测到明显异常")

    # 分析所有表格
    table_dir = "/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/"
    print("开始分析表格图片...")
    print(f"目录: {table_dir}")

    for filename in sorted(os.listdir(table_dir)):
        if filename.endswith('.png'):
            filepath = os.path.join(table_dir, filename)
            results = analyze_image(filepath)
            print_analysis(results)

    print(f"\n{'='*60}")
    print("分析完成!")
    print(f"{'='*60}")



if __name__ == "__main__":
    main()