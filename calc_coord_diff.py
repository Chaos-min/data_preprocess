#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
脚本功能：
    读取两期 TXT 文件（每行格式：ID X Y Z），计算第二期相对于第一期的坐标变化量
    （dX = X2 - X1, dY = Y2 - Y1, dZ = Z2 - Z1），结果保留正负号，单位为米。

使用方法：
    python calc_coord_diff.py 第一期文件.txt 第二期文件.txt [输出文件.txt]
    若未指定输出文件，结果将打印到屏幕。

文件格式示例：
    P001 123.456789 234.567890 345.678901
    P002 124.456789 235.567890 346.678901
    ...
"""

import sys
import os

def read_coords_file(filepath):
    """
    读取坐标文件，返回字典 {ID: (x, y, z)}
    """
    coords = {}
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:          # 跳过空行
                    continue
                parts = line.split()
                if len(parts) != 4:
                    print(f"警告：文件 {filepath} 第 {line_num} 行格式不正确，跳过：{line}")
                    continue
                try:
                    id_str = parts[0]
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    coords[id_str] = (x, y, z)
                except ValueError as e:
                    print(f"警告：文件 {filepath} 第 {line_num} 行数据转换失败，跳过：{line}，错误：{e}")
                    continue
    except FileNotFoundError:
        print(f"错误：文件 {filepath} 不存在。")
        sys.exit(1)
    except Exception as e:
        print(f"读取文件 {filepath} 时发生未知错误：{e}")
        sys.exit(1)
    return coords

def calculate_differences(coords1, coords2):
    """
    计算两个字典的坐标差（第二期 - 第一期），
    返回字典 {ID: (dX, dY, dZ)} 以及缺失 ID 列表。
    """
    diffs = {}
    missing_in_first = []
    missing_in_second = []

    # 获取所有 ID 的并集
    all_ids = set(coords1.keys()) | set(coords2.keys())

    for id_str in all_ids:
        if id_str not in coords1:
            missing_in_first.append(id_str)
            continue
        if id_str not in coords2:
            missing_in_second.append(id_str)
            continue
        x1, y1, z1 = coords1[id_str]
        x2, y2, z2 = coords2[id_str]
        diffs[id_str] = (x2 - x1, y2 - y1, z2 - z1)

    return diffs, missing_in_first, missing_in_second

def main():
    # 检查命令行参数
    if len(sys.argv) < 3:
        print("用法：python calc_coord_diff.py <第一期文件> <第二期文件> [输出文件]")
        print("示例：python calc_coord_diff.py phase1.txt phase2.txt result.txt")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    out_file = sys.argv[3] if len(sys.argv) >= 4 else None

    # 读取两期数据
    print(f"正在读取 {file1} ...")
    coords1 = read_coords_file(file1)
    print(f"成功读取 {len(coords1)} 个点。")

    print(f"正在读取 {file2} ...")
    coords2 = read_coords_file(file2)
    print(f"成功读取 {len(coords2)} 个点。")

    # 计算差值
    diffs, missing1, missing2 = calculate_differences(coords1, coords2)

    # 输出警告信息（缺失 ID）
    if missing1:
        print(f"警告：以下 ID 在第二期文件中存在，但第一期文件中缺失：{', '.join(missing1)}")
    if missing2:
        print(f"警告：以下 ID 在第一期文件中存在，但第二期文件中缺失：{', '.join(missing2)}")

    # 准备输出内容（按 ID 排序）
    if diffs:
        # 按 ID 字母/数字顺序排序
        sorted_ids = sorted(diffs.keys())
        lines = []
        # 表头
        header = "ID          dX (m)       dY (m)       dZ (m)"
        lines.append(header)
        lines.append("-" * 50)
        for id_str in sorted_ids:
            dx, dy, dz = diffs[id_str]
            # 保留 6 位小数，对齐输出（可根据需要调整格式）
            lines.append(f"{id_str:12s} {dx:12.6f} {dy:12.6f} {dz:12.6f}")
    else:
        lines = ["没有可计算的公共点。"]

    # 输出结果
    if out_file:
        try:
            with open(out_file, 'w', encoding='utf-8') as f:
                f.write("\n".join(lines))
            print(f"结果已保存到 {out_file}")
        except Exception as e:
            print(f"写入输出文件失败：{e}")
            sys.exit(1)
    else:
        print("\n计算结果（第二期 - 第一期）：")
        print("\n".join(lines))

if __name__ == "__main__":
    main()