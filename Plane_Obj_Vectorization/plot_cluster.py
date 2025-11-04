# -*- coding：utf-8 -*-
# @Author: Peitong Li
# @Email: 2311797@tongji.edu.cn
# @Time: 2025/08/26
# @Description: 根据main.cpp输出的csv数据绘制图形生成的结果

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def generate_n_colors(n):
    """生成 n 种不同的颜色"""
    cmap = plt.cm.get_cmap('tab20b', n)  # 使用 'tab20b' 色彩映射，可以根据需要选择其他色彩映射
    colors = [cmap(i) for i in range(n)]
    return colors


def generate_40_colors():
    """生成40种具有区分性的颜色"""
    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
        '#ad494a', '#8c6d31', '#843c39', '#e66101', '#55a868', '#0072b2', '#009e73', '#f0e442', '#cc79a7', '#0072b2',
        '#d55e00', '#009e73', '#f0e442', '#cc79a7', '#f0e442', '#0072b2', '#d55e00', '#009e73', '#f0e442', '#cc79a7',
        '#f0e442', '#0072b2', '#d55e00', '#009e73', '#f0e442', '#cc79a7', '#f0e442', '#0072b2', '#d55e00', '#009e73'
    ]
    return colors


def plot_line(ax, csv_file, colors, type=0):
    """
    在指定的子图轴对象上绘制一条线，数据从 CSV 文件中读取。
    
    参数:
        ax: 子图的轴对象 (matplotlib.axes.Axes)
        csv_file: CSV 文件路径
    """
    # 读取 CSV 文件
    # 读取数据，指定列名
    df = pd.read_csv(csv_file, header=None, names=['line_id', 'sx', 'sy', 'ex', 'ey', 'cluster'])
    class_ids = []
    for i, line_id in enumerate(df['line_id'].unique()):
        line_data = df[df['line_id'] == line_id]
        for _, row in line_data.iterrows():
            # sx, sy, ex, ey, cluster = row[['sx', 'sy', 'ex', 'ey', 'cluster']]
            # 转换坐标系, 让图形变得更合理
            sy, sx, ey, ex, cluster = row[['sx', 'sy', 'ex', 'ey', 'cluster']]
            sy = -sy
            ey = -ey
            cluster = int(cluster)
            if type == 0:
                ax.plot([sx, ex], [sy, ey], marker='.', color=colors[i % len(colors)], label=f'Line {line_id}')
            if type == 1:
                if cluster not in class_ids:
                    class_ids.append(cluster)
                    ax.plot([sx, ex], [sy, ey], marker='.', color=colors[cluster % len(colors)], label=f'class {str(cluster)}')
                else:
                    ax.plot([sx, ex], [sy, ey], marker='.', color=colors[cluster % len(colors)])
            if type == 2:
                ax.plot([sx, ex], [sy, ey], marker='.', color=colors[cluster % len(colors)], label=f'Line {line_id}')
            if type == 3:
                ax.plot([sx, ex], [sy, ey], color=colors[cluster % len(colors)], label=f'Line {line_id}')
    
    # 在子图中绘制线并在线的两端绘制点
    ax.set_title(f"Line from {csv_file.split('/')[-1]}")  # 设置子图标题
    # 添加图例并将其放置在图像外部
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1, fontsize='small')
    ax.axis('equal')  # 保持图像比例



def main():
    csv_path = "../info/res/loop/"
    csv_filenames = os.listdir(csv_path)
    for csv_filename in csv_filenames:
        if csv_filename != "test.csv": continue
        colors = generate_40_colors()
        # 创建一个图形和四个子图
        fig, axs = plt.subplots(1, 4, figsize=(12.8, 6.4))

        # 遍历每个 CSV 文件并绘制线
        plot_line(axs[0], "../info/csv/" + csv_filename, colors)
        plot_line(axs[1], "../info/res/cluster/" + csv_filename, colors, type=1)
        plot_line(axs[2], "../info/res/erased/" + csv_filename, colors, type=2)
        plot_line(axs[3], "../info/res/loop/" + csv_filename, colors, type=3)

        
        # 调整布局以防止图例遮挡图像
        plt.tight_layout()
        # 显示图形
        plt.savefig(f'../info/res/{csv_filename.split(".")[0]}.png', dpi=300, bbox_inches='tight')
        print("save fig: ", f'../info/res/{csv_filename.split(".")[0]}.png')


if __name__ == '__main__':
    main()

