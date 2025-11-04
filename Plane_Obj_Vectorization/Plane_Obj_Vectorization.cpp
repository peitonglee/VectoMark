// -*- coding：utf-8 -*-
// @Author: Peitong Li
// @Email: 2311797@tongji.edu.cn
// @Time: 2025/08/26
// @Description: 量平面对象生成核心代码 -- 矢量对象生成算法


#include "Plane_Obj_Vectorization.h"

// 距离
double euclidean_distance(const Simple_Point &a, const Simple_Point &b) { return std::hypot(a.x - b.x, a.y - b.y); }

double calculateThreePointArea(const Simple_Point &A, const Simple_Point &B, const Simple_Point &C) {
    // 向量叉积法计算面积
    double area = 0.5 * std::abs((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x));
    return area;
}

// 两直线交点, abc, def是两条直线的方程参数
std::pair<double, double> _solveLinearEquations(double a, double b, double c, double d, double e, double f) {
    double determinant = a * d - b * c;

    if (determinant != 0) {
        double x = (e * d - b * f) / determinant;
        double y = (a * f - e * c) / determinant;
        return std::make_pair(x, y);
    } else {
        // Handle the case when the determinant is 0 (parallel lines or identical lines)
        // Return a special value to indicate no intersection
        return std::make_pair(NAN, NAN); // NAN represents "Not a Number"
    }
}

// 计算两条直线的交点
std::pair<double, double> _calculateIntersection(const Simple_Line &l1, const Simple_Line &l2) {
    double x1 = l1.sx, y1 = l1.sy, x2 = l1.ex, y2 = l1.ey;
    double x3 = l2.sx, y3 = l2.sy, x4 = l2.ex, y4 = l2.ey;

    double a = y2 - y1;
    double b = x1 - x2;
    double c = x1 * (y1 - y2) + y1 * (x2 - x1);
    double d = y4 - y3;
    double e = x3 - x4;
    double f = x3 * (y3 - y4) + y3 * (x4 - x3);

    return _solveLinearEquations(a, b, d, e, -c, -f);
}

double get_angle(double k1, double k2) {
    // 重合部分占比：重合/线长
    return atan(fabs((k1 - k2) / (1 + k1 * k2))) * 180 / M_PI;
}

// 数据处理
std::tuple<std::vector<Simple_Point>, std::unordered_map<int, int>, std::vector<int>> Plane_Obj_Vectorization::process_lines(const std::vector<Simple_Line> &lines) {
    std::vector<Simple_Point> points;
    std::unordered_map<int, int> pair_map;
    std::vector<int> edge_ids;

    int point_id = 0;
    double sum = 0.0;
    double min_l = 10.0;
    double max_l = 0.0;
    for (size_t i = 0; i < lines.size(); ++i) {
        points.push_back(lines[i].start);
        edge_ids.push_back(lines[i].id);
        points.push_back(lines[i].end);
        edge_ids.push_back(lines[i].id);

        pair_map[point_id] = point_id + 1;
        pair_map[point_id + 1] = point_id;
        point_id += 2;

        sum += lines[i].line_length;
        if (lines[i].line_length < min_l) {
            min_l = lines[i].line_length;
        }
        if (lines[i].line_length > max_l) {
            max_l = lines[i].line_length;
        }
    }
    _mean_length = sum / lines.size();
    _min_length = min_l;
    _max_length = max_l;
    return {points, pair_map, edge_ids};
}

void Plane_Obj_Vectorization::set_A_star() { _A_star_flag = true; }

// 算法调用: 点多点少的情况
std::pair<std::vector<int>, double> Plane_Obj_Vectorization::find_optimal_path(const std::vector<Simple_Point> &points, const std::unordered_map<int, int> &pair_map) {
    const int n = points.size();
    if (n > 64)
        throw std::runtime_error("Exceeding maximum 64 points support");

    // if (_A_star_flag) {
    //     return AStar(points, pair_map);
    // }
    return find_approximate_path(points, pair_map);
    // A* 算法比较准确，但当点数变多时，计算时间变长
}

// 核心算法: 贪心+局部优化
std::pair<std::vector<int>, double> Plane_Obj_Vectorization::find_approximate_path(const std::vector<Simple_Point> &points, const std::unordered_map<int, int> &pair_map) {
    const int n = points.size();

    // 预分配内存
    std::vector<int> path;
    path.reserve(n + 1);
    path.push_back(0);

    // 使用位图代替布尔数组，提高缓存命中率
    // 创建足够多的标记位
    std::vector<uint64> visited_bits((n + 63) / 64, 0);
    // 标记起点已访问
    visited_bits[0] |= 1ULL;

    // 预计算距离矩阵 - 使用OpenMP加速
    std::vector<std::vector<double>> dist_matrix(n, std::vector<double>(n, 0.0));

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        // 对称矩阵
        for (int j = i + 1; j < n; ++j) {
            double dist = euclidean_distance(points[i], points[j]);
            dist_matrix[i][j] = dist;
            dist_matrix[j][i] = dist;
        }
    }

    // 贪心构造初始解: 每一步都获取最近的点
    while (path.size() < n) {
        // 路径最后的节点序号
        int last = path.back();
        // 保存下一个最优节点
        int next_best = -1;
        // 用来保存最小距离
        double min_dist = std::numeric_limits<double>::max();

        // 优先考虑配对点: 线的另一个端点
        if (pair_map.count(last)) {
            // 获取到线段另一个端点节点
            int pair_node = pair_map.at(last);
            // 1. 确定在哪个子图中visited_bits[pair_node / 64]
            // 2. 创建掩码: 1ULL << (pair_node % 64)
            // 3. 按位与
            // 4. 取反表示未访问过
            if (!(visited_bits[pair_node / 64] & (1ULL << (pair_node % 64)))) {
                next_best = pair_node;
            }
        }

        // 如果没有可用的配对点，寻找最近点
        if (next_best == -1) {
            // 当前节点与其余所有节点的距离list, 使用局部变量缓存, 减少访问开销
            const auto &distances = dist_matrix[last];

            for (int i = 0; i < n; ++i) {
                // 若该节点未访问、距离是最近的, 则作为下个节点
                if (!(visited_bits[i / 64] & (1ULL << (i % 64))) && distances[i] < min_dist) {
                    min_dist = distances[i];
                    next_best = i;
                }
            }
        }
        // 如果找到下一个节点
        if (next_best != -1) {
            path.push_back(next_best);
            // 标记一下
            visited_bits[next_best / 64] |= (1ULL << (next_best % 64));
        }
    }

    // 回到起点
    path.push_back(0);

    // 局部优化：2-opt 交换
    bool improved = true;
    double total_cost = 0.0;

    // 计算初始路径长度
    for (size_t i = 0; i < path.size() - 1; ++i) {
        total_cost += dist_matrix[path[i]][path[i + 1]];
    }

    // 限制迭代次数，确保时间在1-10ms范围内
    int max_iterations = 3;
    int iteration = 0;

    while (improved && iteration < max_iterations) {
        improved = false;
        iteration++;

        for (size_t i = 1; i < path.size() - 2; ++i) {
            for (size_t j = i + 1; j < path.size() - 1; ++j) {
                // 检查是否可以通过交换两条边来减少总距离
                int a = path[i - 1], b = path[i];
                int c = path[j], d = path[j + 1];

                // 考虑配对点约束
                bool can_swap = true;
                if (pair_map.count(a) && pair_map.at(a) == b)
                    can_swap = false;
                if (pair_map.count(c) && pair_map.at(c) == d)
                    can_swap = false;

                if (can_swap) {
                    double current_dist = dist_matrix[a][b] + dist_matrix[c][d];
                    double new_dist = dist_matrix[a][c] + dist_matrix[b][d];

                    if (new_dist < current_dist) {
                        // 反转子路径
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

// 核心算法: A* 算法
std::pair<std::vector<int>, double> Plane_Obj_Vectorization::AStar(const std::vector<Simple_Point> &points, const std::unordered_map<int, int> &pair_map) {
    const int n = points.size();
    if (n > 64)
        throw std::runtime_error("Exceeding maximum 64 points support");

    // 预计算距离矩阵
    std::vector<std::vector<double>> dist_matrix(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dist_matrix[i][j] = euclidean_distance(points[i], points[j]);
        }
    }

    // 优先队列初始化: 会根据代价自动排序
    std::priority_queue<State> pq;
    pq.push({0.0, 0, 1ull << 0, {0}});

    std::unordered_map<StateKey, double, Hash> state_map;
    state_map[{0, 1ull << 0}] = 0.0; // 0号点

    std::vector<int> best_path;
    double min_cost = std::numeric_limits<double>::max();

    while (!pq.empty()) {
        State current = pq.top();
        pq.pop();

        // 状态剪枝检查
        StateKey key{current.current, current.mask};
        if (current.cost > state_map[key]) // 比较花销
            continue;

        // 终止条件处理: 1ull << n 表示左移n为, 相当与 1 * 10**n, 判断是否所有点都被访问过
        if (current.mask == (1ull << n) - 1) {
            double final_cost = current.cost + dist_matrix[current.current][0];
            if (final_cost < min_cost) {
                min_cost = final_cost;
                best_path = current.path;
                best_path.push_back(0);
            }
            continue;
        }

        // 处理配对约束
        if (pair_map.count(current.current)) {
            int pair_node = pair_map.at(current.current);
            if (!(current.mask & (1ull << pair_node))) {
                State new_state;
                new_state.current = pair_node;
                // 表示该位被访问过
                new_state.mask = current.mask | (1ull << pair_node);
                // 获取当前点到直线另一端的距离
                new_state.cost = current.cost + dist_matrix[current.current][pair_node];
                // 给当前这个点的序号加入路径list
                new_state.path = current.path;
                new_state.path.push_back(pair_node);

                StateKey new_key{new_state.current, new_state.mask};
                if (!state_map.count(new_key) || new_state.cost < state_map[new_key]) {
                    state_map[new_key] = new_state.cost;
                    pq.push(new_state);
                }
                continue;
            }
        }

        // 自由节点选择: 如当前为0-1, 则将0-1-2, 0-1-3, 0-1-4, 0-1-5, 0-1-6 加入优先队列
        for (int next = 0; next < n; ++next) {
            // 判断该点是否走过
            if (!(current.mask & (1ull << next))) {

                // 新建一个节点
                State new_state;
                // 节点序号
                new_state.current = next;
                // 标记已走过该点
                new_state.mask = current.mask | (1ull << next);
                // 累积花销
                new_state.cost = current.cost + dist_matrix[current.current][next];
                // 路径记录
                new_state.path = current.path;
                new_state.path.push_back(next);
                // 判断该状态之前是否达到过: 没有或比之前花销小
                StateKey new_key{new_state.current, new_state.mask};
                if (!state_map.count(new_key) || new_state.cost < state_map[new_key]) {
                    state_map[new_key] = new_state.cost;
                    pq.push(new_state);
                }
            }
        }
    }

    return {best_path, min_cost};
}

void Plane_Obj_Vectorization::update_endPoint(std::unordered_map<int, Simple_Line> lines, int cur_id, int target_id, Simple_Point &inter, int &mask, int type) {
    Simple_Point cur_line_joint_point = _joints[cur_id].second.first;
    Simple_Point target_line_joint_point = _joints[cur_id].second.second;
    // 获取非更新点, 用更新点和非更新点来代表线
    double cds = euclidean_distance(cur_line_joint_point, lines[cur_id].start);
    double cde = euclidean_distance(cur_line_joint_point, lines[cur_id].end);
    double tds = euclidean_distance(target_line_joint_point, lines[target_id].start);
    double tde = euclidean_distance(target_line_joint_point, lines[target_id].end);
    // 线的非连接点, 即非当前更新的端点
    Simple_Point c_not_joint_point = cds > cde ? lines[cur_id].start : lines[cur_id].end;
    Simple_Point t_not_joint_point = tds > tde ? lines[target_id].start : lines[target_id].end;

    std::pair<double, double> intersection = _calculateIntersection(lines[cur_id], lines[target_id]);

    if (std::isnan(intersection.first) || std::isnan(intersection.second))
        return; // 若两条线平行或重合，则跳过

    // 计算交点到当前直线的两端点距离
    Simple_Point inter_p{intersection.first, intersection.second};
    double inter_2_cj = euclidean_distance(inter_p, cur_line_joint_point);
    double inter_2_cnj = euclidean_distance(inter_p, c_not_joint_point);
    double inter_2_tj = euclidean_distance(inter_p, target_line_joint_point);
    double inter_2_tnj = euclidean_distance(inter_p, t_not_joint_point);
    double angle = get_angle(lines[cur_id].slope, lines[target_id].slope);

    double area1 = calculateThreePointArea(inter_p, cur_line_joint_point, target_line_joint_point);
    double area2 = calculateThreePointArea({0, 0}, {_min_length, 0}, {0, _max_length / 2});
    // 更新点距离两连接点距离
    if (area1 > area2) {
        return;
    }
    // 更新点在两线中间
    if (inter_2_cnj < lines[cur_id].line_length * 0.5 || inter_2_tnj < lines[target_id].line_length * 0.5) {
        return;
    }
    // 更新点在非更新点一侧
    if (inter_2_cj > inter_2_cnj || inter_2_tj > inter_2_tnj) {
        return;
    }
    if (type == 1) {
        // 更新距离不能大于所有线的长度均值, 或者角度够大
        if (angle > 20 || inter_2_cj < 2 * _mean_length) {
            inter = {intersection.first, intersection.second};
            // 交点是给被更新的线的连接点
            mask = cds < cde ? 1 : 2;
        }
    } else {
        // 更新距离不能大于所有线的长度均值, 或者角度够大
        if (angle > 20 || inter_2_tj < 2 * _mean_length) {
            inter = {intersection.first, intersection.second};
            // 交点是给被更新的线
            mask = tds < tde ? 1 : 2;
        }
    }
}

// 获取所有线段的交点和交点是由哪两个线索引对产生的
void Plane_Obj_Vectorization::get_intersections(std::unordered_map<int, Simple_Line> lines) {
    for (auto it = lines.begin(); it != lines.end(); ++it) {
        // 和当前线连接的下一根线id
        int next_line_id = _joints[it->first].first;
        int prev_line_id = _prv_id[it->first];

        Simple_Point next_new_point = Simple_Point(0, 0);
        Simple_Point prev_new_point = Simple_Point(0, 0);
        int mask_next = -1;
        int mask_prev = -1;

        update_endPoint(lines, it->first, next_line_id, next_new_point, mask_next, 1);
        update_endPoint(lines, prev_line_id, it->first, prev_new_point, mask_prev, 0);
        if (mask_next != -1 || mask_prev != -1) {
            _m_line_coord_points[lines[it->first].id] = {lines[it->first].start.x, lines[it->first].start.y, lines[it->first].end.x, lines[it->first].end.y};
        }
        if (next_new_point.x != 0 && next_new_point.y != 0.) {
            // 更新起点、终点不变
            if (mask_next == 1) {
                _m_line_coord_points[it->first][0] = next_new_point.x;
                _m_line_coord_points[it->first][1] = next_new_point.y;
            }
            // 更新终点、起点不变
            if (mask_next == 2) {
                _m_line_coord_points[it->first][2] = next_new_point.x;
                _m_line_coord_points[it->first][3] = next_new_point.y;
            }
        }

        if (prev_new_point.x != 0 && prev_new_point.y != 0.) {
            // 更新起点、终点不变
            if (mask_prev == 1) {
                _m_line_coord_points[it->first][0] = prev_new_point.x;
                _m_line_coord_points[it->first][1] = prev_new_point.y;
            }
            // 更新终点、起点不变
            if (mask_prev == 2) {
                _m_line_coord_points[it->first][2] = prev_new_point.x;
                _m_line_coord_points[it->first][3] = prev_new_point.y;
            }
        }
    }
}

// 入口
void Plane_Obj_Vectorization::loop_arrow(std::vector<Simple_Line> lines) {
    if (lines.size() < 4) {
        _loop_state = false;
        return;
    }
    _m_line_coord_points.clear();
    _parallel_line_pairs.clear();
    _vertical_line_pairs.clear();
    _figure_cycle_info.clear();
    _joints.clear();
    _prv_id.clear();
    _cls_id.clear();

    std::unordered_map<int, Simple_Line> lines_map;
    for (size_t i = 0; i < lines.size(); i++) {
        lines_map[lines[i].id] = lines[i];
        _cls_id[lines[i].cluster] = lines[i].id;
    }

    // 1. 数据处理
    auto [points, pair_map, edge_ids] = process_lines(lines);

    // 2. 最优路径
    auto start_time = std::chrono::high_resolution_clock::now();
    auto [best_path, min_cost] = find_optimal_path(points, pair_map);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "算法耗时: " << duration.count() << " ms\n";

    // 所有线的长度和
    double all_length = 0;
    for (auto it = lines_map.begin(); it != lines_map.end(); ++it) {
        all_length += it->second.line_length;
    }
    // if (min_cost > all_length * 1.5) {
    //     _loop_state = false;
    //     return;
    // }

    // 3. 输出结果
    // std::cout << "Optimal Path: ";
    // for (int idx : best_path) {
    //     std::cout << "P" << idx << " ";
    // }
    // std::cout << "\nMin Cost: " << min_cost << std::endl;

    // 记录线-线之间的连接点 <当前线id, <另一端线id, 交点连接>>
    // 4. 保存line的路径
    // std::cout << "path: ";
    for (size_t i = 0; i < best_path.size(); i++) {
        int from = best_path[i];
        if (((i % 2) - 1) == 0) {
            _joints[edge_ids[from]] = {edge_ids[best_path[i + 1]], {points[best_path[i]], points[best_path[i + 1]]}};
            // std::cout << "L" << edge_ids[from] << "->";
        }
        if (i % 2 == 0 && i != 0) {
            _prv_id[edge_ids[from]] = edge_ids[best_path[i - 1]];
        }
    }

    // std::cout << "L" << edge_ids[0] << std::endl;

    // 6. 获取线的更新点
    get_intersections(lines_map);

    // 获取平行线对
    for (size_t i = 0; i < lines.size() - 1; i++) {
        for (size_t j = i + 1; j < lines.size(); j++) {
            double angle = get_angle(lines[i].slope, lines[j].slope);
            if (angle < 2) {
                _parallel_line_pairs.push_back({lines[i].id, lines[j].id});
            }
            if (angle > 88) {
                _vertical_line_pairs.push_back({lines[i].id, lines[j].id});
            }
        }
    }

    // 获得完美闭环
    int line_num = 0;
    for (size_t i = 0; i < best_path.size() - 1; i++) {
        // 优化路径上的当前点
        Simple_Point c_p = points[i];
        // 对应的线id
        int c_line_id = edge_ids[best_path[i]];
        // 寻找和当前点对应的优化后的点
        std::vector<double> lineEndPoint;
        if (_m_line_coord_points.count(c_line_id)) {
            lineEndPoint = _m_line_coord_points[c_line_id];
        } else {
            lineEndPoint = {lines_map[c_line_id].start.x, lines_map[c_line_id].start.y, lines_map[c_line_id].end.x, lines_map[c_line_id].end.y};
        }
        Simple_Point start_point{lineEndPoint[0], lineEndPoint[1]};
        Simple_Point end_point{lineEndPoint[2], lineEndPoint[3]};
        double d1 = euclidean_distance(c_p, start_point);
        double d2 = euclidean_distance(c_p, end_point);
        Simple_Point start_end_flag = d1 < d2 ? start_point : end_point;

        if (((i % 2) - 1) == 0) {
            _figure_cycle_info[line_num].second.push_back(start_end_flag.x);
            _figure_cycle_info[line_num].second.push_back(start_end_flag.y);
            line_num++;
        }
        if (i % 2 == 0) {
            std::vector<double> p;
            p.push_back(start_end_flag.x);
            p.push_back(start_end_flag.y);
            _figure_cycle_info.push_back({c_line_id, p});
        }
    }
    _loop_state = true;
}

