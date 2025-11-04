// -*- coding：utf-8 -*-
// @Author: Peitong Li
// @Email: 2311797@tongji.edu.cn
// @Time: 2025/08/26
// @Description: 模拟基于线的路径规划, 寻找最短路径


#include "matplotlibcpp.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <omp.h>
#include <queue>
#include <unordered_map>
#include <vector>

using uint64 = unsigned long long;
namespace plt = matplotlibcpp;

// 基本数据结构
struct Simple_Point
{
    double x;
    double y;
    // 计算两点之间的欧氏距离
    double distance(const Simple_Point &other) const { return std::sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y)); }
};

struct Simple_Line
{
    int id;
    double line_length;
    Simple_Point start;
    Simple_Point end;
    Simple_Line() : id(0), line_length(0.0) {}
    Simple_Line(int _id, const Simple_Point &_start, const Simple_Point &_end) : id(_id), start(_start), end(_end) { line_length = start.distance(end); }
};

// 状态数据结构
struct State
{
    double cost; // 累计成本
    int current; // 当前节点索引
    uint64 mask; // 访问掩码（64位）

    bool operator<(const State &other) const
    {
        return cost > other.cost; // 最小堆
    }
};

// 辅助函数
double euclidean_distance(const Simple_Point &a, const Simple_Point &b) { return std::hypot(a.x - b.x, a.y - b.y); }

// 预计算距离矩阵
std::vector<std::vector<double>> precompute_distances(const std::vector<Simple_Point> &points)
{
    const int n = points.size();
    std::vector<std::vector<double>> dist_matrix(n, std::vector<double>(n));

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            dist_matrix[i][j] = euclidean_distance(points[i], points[j]);
        }
    }

    return dist_matrix;
}

// 贪心算法计算初始解
double calculate_greedy_solution(const std::vector<std::vector<double>> &dist_matrix, int n, const std::unordered_map<int, int> &pair_map)
{
    std::vector<bool> visited(n, false);
    visited[0] = true;
    int current = 0;
    double total_cost = 0.0;
    for (int i = 1; i < n; ++i)
    {
        int next = -1;
        double min_dist = std::numeric_limits<double>::max();

        // 优先考虑配对点
        if (pair_map.count(current) && !visited[pair_map.at(current)])
        {
            next = pair_map.at(current);
        }
        else
        {
            for (int j = 0; j < n; ++j)
            {
                if (!visited[j] && dist_matrix[current][j] < min_dist)
                {
                    min_dist = dist_matrix[current][j];
                    next = j;
                }
            }
        }

        if (next != -1)
        {
            visited[next] = true;
            total_cost += dist_matrix[current][next];
            current = next;
        }
    }

    // 回到起点
    total_cost += dist_matrix[current][0];

    return total_cost;
}

// 数据处理函数
std::tuple<std::vector<Simple_Point>, std::unordered_map<int, int>, std::vector<int>> process_lines(const std::vector<Simple_Line> &lines)
{
    std::vector<Simple_Point> points;
    std::unordered_map<int, int> pair_map;
    std::vector<int> edge_ids;

    int point_id = 0;
    for (const auto &line : lines)
    {
        points.push_back(line.start);
        edge_ids.push_back(line.id);
        points.push_back(line.end);
        edge_ids.push_back(line.id);

        pair_map[point_id] = point_id + 1;
        pair_map[point_id + 1] = point_id;
        point_id += 2;
    }
    return {points, pair_map, edge_ids};
}

// 贪心优化算法
std::pair<std::vector<int>, double> find_approximate_path(const std::vector<Simple_Point> &points, const std::unordered_map<int, int> &pair_map)
{
    const int n = points.size();

    // 预分配内存
    std::vector<int> path;
    path.reserve(n + 1);
    path.push_back(0);

    // 使用位图代替布尔数组，提高缓存命中率
    std::vector<uint64> visited_bits((n + 63) / 64, 0);
    visited_bits[0] |= 1ULL; // 标记起点已访问

    // 预计算距离矩阵 - 使用OpenMP加速
    std::vector<std::vector<double>> dist_matrix(n, std::vector<double>(n, 0.0));

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double dist = euclidean_distance(points[i], points[j]);
            dist_matrix[i][j] = dist;
            dist_matrix[j][i] = dist; // 对称矩阵
        }
    }

    // 贪心构造初始解
    while (path.size() < n)
    {
        int last = path.back();
        int next_best = -1;
        double min_dist = std::numeric_limits<double>::max();

        // 优先考虑配对点
        if (pair_map.count(last))
        {
            int pair_node = pair_map.at(last);
            if (!(visited_bits[pair_node / 64] & (1ULL << (pair_node % 64))))
            {
                next_best = pair_node;
            }
        }

        // 如果没有可用的配对点，寻找最近点
        if (next_best == -1)
        {
            // 使用局部变量缓存，减少访问开销
            const auto &distances = dist_matrix[last];

            for (int i = 0; i < n; ++i)
            {
                if (!(visited_bits[i / 64] & (1ULL << (i % 64))) && distances[i] < min_dist)
                {
                    min_dist = distances[i];
                    next_best = i;
                }
            }
        }

        if (next_best != -1)
        {
            path.push_back(next_best);
            visited_bits[next_best / 64] |= (1ULL << (next_best % 64));
        }
    }

    // 回到起点
    path.push_back(0);

    // 注: 本来一遍贪心算法就可以完成闭合的任务，但由于一些线的长度或缺失线的问题，导致贪心算法失效，最终的距离总和不是最优解，需要进行优化

    // 局部优化：2-opt 交换
    bool improved = true;
    double total_cost = 0.0;

    // 计算初始路径长度
    for (size_t i = 0; i < path.size() - 1; ++i)
    {
        total_cost += dist_matrix[path[i]][path[i + 1]];
    }

    // 限制迭代次数，确保时间在1-10ms范围内
    int max_iterations = 3;
    int iteration = 0;

    while (improved && iteration < max_iterations)
    {
        improved = false;
        iteration++;

        for (size_t i = 1; i < path.size() - 2; ++i)
        {
            for (size_t j = i + 1; j < path.size() - 1; ++j)
            {
                // 检查是否可以通过交换两条边来减少总距离
                int a = path[i - 1], b = path[i];
                int c = path[j], d = path[j + 1];

                // 考虑配对点约束: 不能交换配对点, 即线必须是闭合的
                bool can_swap = true;
                if (pair_map.count(a) && pair_map.at(a) == b)
                    can_swap = false;
                if (pair_map.count(c) && pair_map.at(c) == d)
                    can_swap = false;

                if (can_swap)
                {
                    double current_dist = dist_matrix[a][b] + dist_matrix[c][d];
                    double new_dist = dist_matrix[a][c] + dist_matrix[b][d];

                    if (new_dist < current_dist)
                    {
                        // 交换path中的两个片段路径
                        std::reverse(path.begin() + i, path.begin() + j + 1);
                        total_cost = total_cost - current_dist + new_dist;
                        improved = true;
                        break;
                    }
                }
            }
            if (improved)
                break;
        }
    }

    return {path, total_cost};
}

std::pair<std::vector<int>, double> AStar(const std::vector<Simple_Point> &points, const std::unordered_map<int, int> &pair_map)
{
    const int n = points.size();
    if (n > 64)
        throw std::runtime_error("Exceeding maximum 64 points support");

    std::vector<std::vector<double>> dist_matrix(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            dist_matrix[i][j] = euclidean_distance(points[i], points[j]);
        }
    }

    std::priority_queue<State> pq;
    pq.push({0.0, 0, 1ull << 0});

    struct StateKey
    {
        int current;
        uint64 mask;
        bool operator==(const StateKey &other) const { return current == other.current && mask == other.mask; }
    };

    struct Hash
    {
        size_t operator()(const StateKey &k) const { return std::hash<int>()(k.current) ^ (std::hash<uint64>()(k.mask) << 1); }
    };

    struct StateInfo
    {
        double cost;
        int predecessor;
    };

    std::unordered_map<StateKey, StateInfo, Hash> state_map;
    state_map[{0, 1ull << 0}] = {0.0, -1};

    double min_cost = std::numeric_limits<double>::max();
    int best_end_node = -1;
    uint64 best_end_mask = 0;

    while (!pq.empty())
    {
        State current = pq.top();
        pq.pop();

        StateKey key{current.current, current.mask};
        if (current.cost > state_map[key].cost)
            continue;

        if (current.mask == (1ull << n) - 1)
        {
            double final_cost = current.cost + dist_matrix[current.current][0];
            if (final_cost < min_cost)
            {
                min_cost = final_cost;
                best_end_node = current.current;
                best_end_mask = current.mask;
            }
            continue;
        }

        if (pair_map.count(current.current))
        {
            int pair_node = pair_map.at(current.current);
            if (!(current.mask & (1ull << pair_node)))
            {
                State new_state;
                new_state.current = pair_node;
                new_state.mask = current.mask | (1ull << pair_node);
                new_state.cost = current.cost + dist_matrix[current.current][pair_node];

                double estimated_total = new_state.cost + dist_matrix[pair_node][0];
                if (estimated_total >= min_cost)
                    continue;

                StateKey new_key{new_state.current, new_state.mask};
                if (!state_map.count(new_key) || new_state.cost < state_map[new_key].cost)
                {
                    state_map[new_key] = {new_state.cost, current.current};
                    pq.push(new_state);
                }
                continue;
            }
        }

        for (int next = 0; next < n; ++next)
        {
            if (!(current.mask & (1ull << next)))
            {
                State new_state;
                new_state.current = next;
                new_state.mask = current.mask | (1ull << next);
                new_state.cost = current.cost + dist_matrix[current.current][next];

                double estimated_total = new_state.cost + dist_matrix[next][0];
                if (estimated_total >= min_cost)
                    continue;

                StateKey new_key{new_state.current, new_state.mask};
                if (!state_map.count(new_key) || new_state.cost < state_map[new_key].cost)
                {
                    state_map[new_key] = {new_state.cost, current.current};
                    pq.push(new_state);
                }
            }
        }
    }

    std::vector<int> best_path;
    if (best_end_node != -1)
    {
        best_path.push_back(0);
        int current_node = best_end_node;
        uint64 current_mask = best_end_mask;

        while (current_node != 0)
        {
            best_path.push_back(current_node);
            StateKey key{current_node, current_mask};
            auto it = state_map.find(key);
            if (it == state_map.end())
                break;
            int prev_node = it->second.predecessor;
            current_mask ^= (1ull << current_node);
            current_node = prev_node;
        }
        best_path.push_back(0);
        std::reverse(best_path.begin(), best_path.end());
    }

    return {best_path, min_cost};
}

// 调用算法
std::pair<std::vector<int>, double> find_optimal_path(const std::vector<Simple_Point> &points, const std::unordered_map<int, int> &pair_map)
{
    const int n = points.size();
    if (n > 64)
        throw std::runtime_error("Exceeding maximum 64 points support");

    // 当点数超过阈值时，使用近似算法
    // if (n > 20) {
    //     return find_approximate_path(points, pair_map);
    // }

    // return find_approximate_path(points, pair_map);
    return AStar(points, pair_map);
}

void visualize(const std::vector<Simple_Point> &points, const std::vector<int> &path, const std::vector<int> &edge_ids)
{
    plt::figure_size(1200, 800);

    // 绘制所有原始边
    std::unordered_map<int, std::pair<Simple_Point, Simple_Point>> edge_map;
    for (size_t i = 0; i < points.size(); i += 2)
    {
        int edge_id = edge_ids[i];
        edge_map[edge_id] = {points[i], points[i + 1]};
        plt::plot({points[i].x, points[i + 1].x}, {points[i].y, points[i + 1].y}, "gray");
        // 标注边序号
        double mid_x = (points[i].x + points[i + 1].x) / 2;
        double mid_y = (points[i].y + points[i + 1].y) / 2;
        plt::text(mid_x, mid_y, "L" + std::to_string(edge_id));
    }

    // 绘制路径
    for (size_t i = 0; i < path.size() - 1; ++i)
    {
        int from = path[i];
        int to = path[i + 1];
        Simple_Point p1 = points[from];
        Simple_Point p2 = points[to];

        // 判断线段类型
        std::string style = (edge_ids[from] == edge_ids[to]) ? "r-" : "k--";
        plt::plot({p1.x, p2.x}, {p1.y, p2.y}, style);
        plt::pause(0.2); // 动态显示
    }

    // 标注点集
    std::vector<double> x, y;
    for (const auto &p : points)
    {
        x.push_back(p.x);
        y.push_back(p.y);
    }
    plt::scatter(x, y, 50, {{"color", "blue"}, {"edgecolors", "black"}});
    plt::scatter(std::vector<double>{points[0].x}, std::vector<double>{points[0].y}, 100, {{"color", "red"}, {"marker", "*"}});

    plt::title("Optimal Path Visualization");
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::grid(true);
    plt::save("../optimal_path.png");
    std::cout << "Optimal Path Visualization saved to 'Path_Planning_Base_Lines/optimal_path.png'\n";
    plt::show();
}

// 主函数
int main()
{
    // 直行+左转
    std::vector<Simple_Line> lines = {{721, {-0.706955, -0.149191}, {-0.655177, -0.00025847}}, {29, {-0.386368, -0.277053}, {1.35308, -0.21472}}, {30, {0.033409, 0.283683}, {0.338169, -0.0746077}}, {31, {-0.0459915, 0.537551}, {1.11383, 0.358801}}, {430, {0.709603, 0.297978}, {0.993771, -0.0768412}}, {33, {0.978334, -0.0793122}, {1.33925, -0.056238}}, {34, {1.38416, -0.0714839}, {1.39715, -0.214701}}, {61, {-0.488059, -0.126961}, {0.365489, -0.109888}}, {72, {0.70433, 0.317606}, {1.13036, 0.324319}}, {585, {-0.641193, 0.297982}, {0.0416619, 0.49791}}, {761, {-0.600559, 0.25589}, {0.034039, 0.271183}}, {583, {-1.9331, -0.223407}, {-0.687037, 0.0122987}}, {584, {-1.94803, -0.265655}, {-0.873097, -0.42813}}, {0, {-0.75, -0.4}, {-0.75, -0.27}}};

    // 直箭头
    // std::vector<Simple_Line> lines = {{59, {-1.91913, 0.216423}, {0.156583, -0.11701}}, {74, {-1.91849, 0.251759}, {-1.89785, 0.394598}}, {205, {0.209993, 0.203923}, {1.50843, -0.258727}}, {73, {-1.74428, 0.412557}, {-0.149193, 0.12554}}, {207, {0.0195754, 0.242293}, {0.0380241, 0.10464}}, {58, {0.0126466, -0.223691}, {1.27177, -0.240028}}, {206, {-0.0424564, -0.206575}, {-0.00309438, -0.0605256}}};

    // 数据处理
    auto [points, pair_map, edge_ids] = process_lines(lines);

    auto start_time = std::chrono::high_resolution_clock::now();
    // 路径计算
    auto [best_path, min_cost] = find_optimal_path(points, pair_map);
    // 算法结束后
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "算法耗时: " << duration.count() << " 毫秒\n";
    // 结果输出
    std::cout << "Optimal Path: ";
    for (int idx : best_path)
    {
        std::cout << "P" << idx << " ";
    }
    std::cout << "\nMin Distance: " << min_cost << std::endl;

    std::vector<int> looped_arrow;
    std::unordered_map<int, std::pair<int, int>> looped_arrow_order;
    for (size_t i = 0; i < best_path.size() - 1; i += 2)
    {
        int from = best_path[i];
        looped_arrow.push_back(edge_ids[from]);
    }
    std::cout << "Arrow Line Order: ";
    for (size_t i = 0; i < looped_arrow.size(); i++)
    {
        std::cout << "L" << looped_arrow[i] << "--";
        if (i == 0)
        {
            looped_arrow_order[looped_arrow[i]] = {looped_arrow[looped_arrow.size() - 1], looped_arrow[i + 1]};
        }
        else if (i == looped_arrow.size() - 1)
        {
            looped_arrow_order[looped_arrow[i]] = {looped_arrow[i - 1], looped_arrow[0]};
        }
        else
        {
            looped_arrow_order[looped_arrow[i]] = {looped_arrow[i - 1], looped_arrow[i + 1]};
        }
    }
    std::cout << "L" << looped_arrow[0] << "\n";

    // 可视化
    visualize(points, best_path, edge_ids);

    return 0;
}