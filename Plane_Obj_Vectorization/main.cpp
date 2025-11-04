// -*- coding：utf-8 -*-
// @Author: Peitong Li
// @Email: 2311797@tongji.edu.cn
// @Time: 2025/08/26
// @Description: 矢量平面对象生成核心代码 -- 主程序


#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "ldbscan.h"
#include "Plane_Obj_Vectorization.h"
#include "simple_data.h"
#include <filesystem>
#include <matplotlibcpp.h>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;
namespace plt = matplotlibcpp;

std::string filename;
std::string savePath = "../info/res/";

// 生成 n 种色差较大的颜色（返回 HEX 格式字符串，如 "#FF0000"）
std::vector<std::string> generateDistinctColors(int n)
{
    std::vector<std::string> colors;
    if (n <= 0)
        return colors;

    // HSV 参数：固定饱和度(S)和亮度(V)，色相(H)均匀分布
    const float s = 0.8f; // 饱和度 (0-1)
    const float v = 0.9f; // 亮度 (0-1)

    for (int i = 0; i < n; ++i)
    {
        float h = static_cast<float>(i) * (360.0f / n); // 色相均匀分布 (0-360°)

        // HSV 转 RGB
        float c = v * s;
        float x = c * (1 - std::abs(std::fmod(h / 60.0f, 2) - 1));
        float m = v - c;

        float r, g, b;
        if (h < 60)
        {
            r = c;
            g = x;
            b = 0;
        }
        else if (h < 120)
        {
            r = x;
            g = c;
            b = 0;
        }
        else if (h < 180)
        {
            r = 0;
            g = c;
            b = x;
        }
        else if (h < 240)
        {
            r = 0;
            g = x;
            b = c;
        }
        else if (h < 300)
        {
            r = x;
            g = 0;
            b = c;
        }
        else
        {
            r = c;
            g = 0;
            b = x;
        }

        // RGB 转 HEX
        char hex[8];
        snprintf(hex, sizeof(hex), "#%02X%02X%02X", static_cast<int>((r + m) * 255), static_cast<int>((g + m) * 255), static_cast<int>((b + m) * 255));
        colors.push_back(hex);
    }

    return colors;
}

// 获取某文件夹下的所有文件（返回文件名列表）
std::vector<std::string> getFilesInDirectory(const std::string &folderPath)
{
    std::vector<std::string> files;

    try
    {
        for (const auto &entry : fs::directory_iterator(folderPath))
        {
            if (entry.is_regular_file())
            { // 只获取文件，不包含子目录
                files.push_back(entry.path().filename().string());
            }
        }
    }
    catch (const fs::filesystem_error &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return files;
}

void saveToCSV(const std::vector<Simple_Line> &lines, const std::string &filename)
{
    std::ofstream file(filename); // 打开文件
    if (!file.is_open())
    {
        std::cerr << "无法打开文件：" << filename << std::endl;
        return;
    }

    // 遍历 vector 并写入每一行数据
    for (const auto &line : lines)
    {
        file << line.id << "," << line.sx << "," << line.sy << "," << line.ex << "," << line.ey << "," << line.cluster << "\n";
    }

    file.close(); // 关闭文件
    std::cout << "数据已成功保存到文件：" << filename << std::endl;
}

void drawFourConnectedLines(const std::vector<Simple_Point> &lines, const std::string &line_style = "--", const std::string &line_label = "")
{
    if (lines.size() != 4)
    {
        throw std::invalid_argument("lines must contain exactly 4 points!");
    }

    // 定义四条线的连接顺序：0-1, 1-2, 2-3, 3-0
    const std::vector<std::pair<int, int>> connections = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    for (size_t i = 0; i < connections.size(); ++i)
    {
        const auto &[start_idx, end_idx] = connections[i];
        const auto &p1 = lines[start_idx];
        const auto &p2 = lines[end_idx];

        // 提取当前线段的两个点
        std::vector<double> x = {p1.x, p2.x};
        std::vector<double> y = {p1.y, p2.y};

        // 绘制线段（使用相同的线型和标签）
        plt::plot(x, y,
                  {
                      {"linestyle", line_style}, {"c", "black"}, {"label", (i == 0) ? line_label : ""} // 只在第一条线段显示标签
                  });
    }
}

void save_lines(std::vector<Simple_Line> lines, std::vector<Simple_Line> clustered_lines, Plane_Obj_Vectorization loop_ob)
{
    saveToCSV(lines, savePath + "cluster/" + filename + ".csv");
    saveToCSV(clustered_lines, savePath + "erased/" + filename + ".csv");
    std::vector<Simple_Line> lines_looped;
    for (int i = 0; i < clustered_lines.size(); i++)
    {
        if (loop_ob._m_line_coord_points.count(clustered_lines[i].id))
        {
            std::vector<double> coord_points = loop_ob._m_line_coord_points[clustered_lines[i].id];
            Simple_Line new_line = clustered_lines[i];
            new_line.sx = coord_points[0];
            new_line.sy = coord_points[1];
            new_line.ex = coord_points[2];
            new_line.ey = coord_points[3];
            lines_looped.push_back(new_line);
        }
    }
    saveToCSV(lines_looped, savePath + "loop/" + filename + ".csv");

    std::vector<std::string> custom_colors = generateDistinctColors(clustered_lines.size());
    plt::figure_size(1280, 640);
    plt::subplot(1, 1, 1);
    for (int i = 0; i < lines.size(); i++)
    {
        std::vector<double> new_x = {lines[i].sx, lines[i].ex};
        std::vector<double> new_y = {lines[i].sy, lines[i].ey};
        int cluster_num = lines[i].cluster;
        plt::plot(new_x, new_y, {{"c", custom_colors[cluster_num - 1]}, {"marker", ""}, {"label", "line " + std::to_string(lines[i].id)}});

        drawFourConnectedLines(lines[i].bbox_lines);
    }
    plt::title(filename + " box");
    plt::axis("equal");
    plt::legend();
    plt::save(savePath + filename + "_box.png");
    std::cout << "save " << savePath + filename + "_box.png" << std::endl;
}

void plot_lines(std::vector<Simple_Line> lines, std::vector<Simple_Line> clustered_lines, Plane_Obj_Vectorization loop_ob)
{
    std::vector<std::string> custom_colors = generateDistinctColors(clustered_lines.size());

    // save_res
    plt::figure_size(1280, 640);
    plt::subplot(1, 4, 1);
    for (int i = 0; i < lines.size(); i++)
    {
        std::vector<double> new_x = {lines[i].sx, lines[i].ex};
        std::vector<double> new_y = {lines[i].sy, lines[i].ey};
        plt::plot(new_x, new_y, {{"c", custom_colors[i % custom_colors.size()]}, {"marker", "."}, {"label", "Line " + std::to_string(lines[i].id)}});
    }
    plt::title(filename + "orignal");
    plt::axis("equal");
    plt::tight_layout(); // 自动调整布局
    plt::margins(0.1, 0.1);
    plt::legend();

    plt::subplot(1, 4, 2);
    std::map<int, int> cluster_map;
    for (int i = 0; i < lines.size(); i++)
    {
        std::vector<double> new_x = {lines[i].sx, lines[i].ex};
        std::vector<double> new_y = {lines[i].sy, lines[i].ey};
        int cluster_num = lines[i].cluster;
        if (!cluster_map.count(cluster_num))
        {
            cluster_map[cluster_num] = 1;
            plt::plot(new_x, new_y, {{"c", custom_colors[cluster_num - 1]}, {"marker", "."}, {"label", "class " + std::to_string(cluster_num)}});
        }
        else
        {
            plt::plot(new_x, new_y, {{"c", custom_colors[cluster_num - 1]}, {"marker", "."}});
        }
    }
    plt::title(filename + "clustered");
    plt::axis("equal");
    plt::legend();

    plt::subplot(1, 4, 3);
    cluster_map.clear();
    for (int i = 0; i < clustered_lines.size(); i++)
    {
        std::vector<double> new_x = {clustered_lines[i].sx, clustered_lines[i].ex};
        std::vector<double> new_y = {clustered_lines[i].sy, clustered_lines[i].ey};
        plt::plot(new_x, new_y, {{"c", custom_colors[clustered_lines[i].cluster - 1]}, {"marker", "."}, {"label", "line " + std::to_string(clustered_lines[i].id)}});
    }
    plt::title(filename + " erased");
    plt::axis("equal");
    plt::legend();

    plt::subplot(1, 4, 4);
    cluster_map.clear();
    for (int i = 0; i < clustered_lines.size(); i++)
    {
        if (loop_ob._m_line_coord_points.count(clustered_lines[i].id))
        {
            std::vector<double> coord_points = loop_ob._m_line_coord_points[clustered_lines[i].id];
            std::vector<double> new_x = {coord_points[0], coord_points[2]};
            std::vector<double> new_y = {coord_points[1], coord_points[3]};
            plt::plot(new_x, new_y, {{"c", custom_colors[clustered_lines[i].cluster - 1]}, {"label", "line " + std::to_string(clustered_lines[i].id)}});
        }
    }
    plt::title(filename + " loop");
    plt::axis("equal");
    plt::legend();

    plt::save(savePath + filename + ".png");
    std::cout << "save " << savePath + filename + ".png" << std::endl;

    plt::figure_size(1280, 640);
    plt::subplot(1, 1, 1);
    for (int i = 0; i < clustered_lines.size(); i++)
    {
        std::vector<double> new_x = {clustered_lines[i].sx, clustered_lines[i].ex};
        std::vector<double> new_y = {clustered_lines[i].sy, clustered_lines[i].ey};
        int cluster_num = clustered_lines[i].cluster;
        plt::plot(new_x, new_y, {{"c", custom_colors[cluster_num - 1]}, {"marker", ""}, {"label", "line " + std::to_string(clustered_lines[i].id)}});

        drawFourConnectedLines(clustered_lines[i].bbox_lines);
    }
    plt::title(filename + " loop");
    plt::axis("equal");
    plt::legend();
    // plt::save(savePath + filename + "_cluster.png");
    std::cout << "save " << savePath + filename + "_cluster.png" << std::endl;

    plt::figure_size(1280, 640);
    plt::subplot(1, 1, 1);
    for (int i = 0; i < clustered_lines.size(); i++)
    {
        if (loop_ob._m_line_coord_points.count(clustered_lines[i].id))
        {
            std::vector<double> coord_points = loop_ob._m_line_coord_points[clustered_lines[i].id];
            std::vector<double> new_x = {coord_points[0], coord_points[2]};
            std::vector<double> new_y = {coord_points[1], coord_points[3]};
            int cluster_num = clustered_lines[i].cluster;
            plt::plot(new_x, new_y, {{"c", custom_colors[cluster_num - 1]}, {"marker", ""}, {"label", "line " + std::to_string(clustered_lines[i].id)}});

            drawFourConnectedLines(clustered_lines[i].bbox_lines);
        }
    }
    plt::title(filename + " loop");
    plt::axis("equal");
    plt::legend();

    // plt::save(savePath + filename + "_big.png");
    std::cout << "save " << savePath + filename + "_big.png" << std::endl;
    plt::show();
}

int main()
{
    std::string folderPath = "../info/csv/";
    std::vector<std::string> files = getFilesInDirectory(folderPath);
    for (int i = 0; i < files.size(); i++)
    {
        fs::path p(files[i]);
        filename = p.stem().string();
        if (filename != "test")
            continue;
        std::ifstream csv_file("../info/csv/" + filename + ".csv");

        if (!csv_file.is_open())
        {
            std::cerr << "Error: Could not open file 'lines_info.csv'!" << std::endl;
            return -1;
        }

        std::string line;
        Line_Creater line_creater;

        // 1、读取数据
        while (std::getline(csv_file, line))
        {
            std::istringstream iss(line);
            std::string token;
            std::vector<double> tmp_line;

            // 按逗号分割每行数据
            while (std::getline(iss, token, ','))
            {
                tmp_line.push_back(std::stod(token)); // 将字符串转换为 double
            }

            // 确保数据完整（至少 6 列：id, sx, sy, ex, ey, cluster）
            if (tmp_line.size() >= 5)
            {
                int id = static_cast<int>(tmp_line[0]);
                double sx = tmp_line[1];
                double sy = tmp_line[2];
                double ex = tmp_line[3];
                double ey = tmp_line[4];
                std::vector<double> line{sx, sy, ex, ey};

                // 创建 Simple_Line 对象并添加到 vector
                if (!line_creater.lines.count(id))
                {
                    line_creater.lines[id] = line;
                }
            }
        }

        // 2、聚类
        LDBSCAN ldbscan = LDBSCAN();
        std::vector<Simple_Line> lines = line_creater.gen_lines(ldbscan._box_r_scale_ratio);

        ldbscan.set_minimum_length(line_creater.minimum_length); // 动态设置 DBSCAN 算法的半径
        ldbscan.ldbscan(lines);

        std::vector<int> _line_labels_after_ldbscan;
        std::vector<int> _longest_line_index;
        for (size_t i = 0; i < lines.size(); i++)
        {
            _line_labels_after_ldbscan.push_back(lines[i].cluster);
        }
        std::vector<Simple_Line> clustered_lines = ldbscan.findLongestSegments(lines, _line_labels_after_ldbscan, _longest_line_index);

        // 3、矢量化
        Plane_Obj_Vectorization loop_obj = Plane_Obj_Vectorization();
        loop_obj.loop_arrow(clustered_lines);
        if (!loop_obj._loop_state)
        {
            return false;
        }
        std::vector<int> _better_line_index;
        std::vector<Simple_Line> best_lines = ldbscan.findBestSegments(lines, _better_line_index, loop_obj);
        loop_obj.loop_arrow(best_lines);

        plot_lines(lines, clustered_lines, loop_obj);
        save_lines(lines, best_lines, loop_obj);
    }
    system("python ../plot_cluster.py");
    return 0;
}