// -*- coding：utf-8 -*-
// @Author: Peitong Li
// @Email: 2311797@tongji.edu.cn
// @Time: 2025/08/26
// @Description: 矢量平面对象生成核心代码 -- 线聚类算法


#include "ldbscan.h"
#include <cstdlib>

// 两点间欧式距离
double LDBSCAN::euclideanDistance(Simple_Point p1, Simple_Point p2) { return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2)); }

// 两线段的交点
Simple_Point LDBSCAN::findIntersection(Simple_Segment l1, Simple_Segment l2) {
    double a1 = l1.end.y - l1.start.y;
    double b1 = l1.start.x - l1.end.x;
    double c1 = a1 * l1.start.x + b1 * l1.start.y;

    double a2 = l2.end.y - l2.start.y;
    double b2 = l2.start.x - l2.end.x;
    double c2 = a2 * l2.start.x + b2 * l2.start.y;

    double determinant = a1 * b2 - a2 * b1;
    double x = (b2 * c1 - b1 * c2) / determinant;
    double y = (a1 * c2 - a2 * c1) / determinant;
    return {x, y};
}

// 计算二维空间中三个点构成的三角形面积
double calculateTriangleArea(const Simple_Point &A, const Simple_Point &B, const Simple_Point &C) {
    // 向量叉积法计算面积
    double area = 0.5 * std::abs((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x));
    return area;
}

// 计算点p是否在矩形rect内
bool LDBSCAN::isPointInsideRectangle(Simple_Point P, std::vector<Simple_Point> polygon, double area) {
    double a1 = calculateTriangleArea(P, polygon[0], polygon[1]);
    double a2 = calculateTriangleArea(P, polygon[1], polygon[2]);
    double a3 = calculateTriangleArea(P, polygon[2], polygon[3]);
    double a4 = calculateTriangleArea(P, polygon[3], polygon[0]);
    double sum_a = a1 + a2 + a3 + a4;
    if (std::abs(sum_a - area) <= 1e-4)
        return true;
    else
        return false;
}

// 判断两个线段是否相交的函数
bool LDBSCAN::doIntersect(Simple_Segment l1, Simple_Segment l2) {
    Simple_Point A = l1.start;
    Simple_Point B = l1.end;
    Simple_Point C = l2.start;
    Simple_Point D = l2.end;

    // 计算分母
    double denominator = (B.x - A.x) * (D.y - C.y) - (B.y - A.y) * (D.x - C.x);

    // 如果分母为0，说明两线段平行或共线
    if (denominator == 0) {
        return false; // 这里简化处理，不考虑共线情况
    }

    // 计算 t 和 s
    double t = ((C.x - A.x) * (D.y - C.y) - (C.y - A.y) * (D.x - C.x)) / denominator;
    double s = ((C.x - A.x) * (B.y - A.y) - (C.y - A.y) * (B.x - A.x)) / denominator;

    // 如果 t 和 s 都在 [0, 1] 之间，则线段相交
    return (t >= 0 && t <= 1 && s >= 0 && s <= 1);
}

// 计算线段l2在线段l1形成的矩形中的长度
double LDBSCAN::calc_overlap_area_dist(Simple_Line l1, Simple_Line l2) {
    double dist = 0.0;
    if (isPointInsideRectangle(l2.start, l1.bbox_lines, l1.bbox_area) && isPointInsideRectangle(l2.end, l1.bbox_lines, l1.bbox_area)) {
        // 线段在矩形内
        dist = l2.line_length;
    } else {
        // 线段不在矩形内，看有几个交点，最多两个交点，计算交线长度
        std::vector<Simple_Point> intersections;

        Simple_Segment s1 = {l1.bbox_lines[0], l1.bbox_lines[1]};
        Simple_Segment s2 = {l1.bbox_lines[1], l1.bbox_lines[2]};
        Simple_Segment s3 = {l1.bbox_lines[2], l1.bbox_lines[3]};
        Simple_Segment s4 = {l1.bbox_lines[3], l1.bbox_lines[0]};
        Simple_Segment s = {l2.start, l2.end};
        if (doIntersect(s1, s)) {
            Simple_Point i1 = findIntersection(s1, s);
            intersections.push_back(i1);
        }
        if (doIntersect(s2, s)) {
            Simple_Point i2 = findIntersection(s2, s);
            intersections.push_back(i2);
        }

        if (doIntersect(s3, s)) {
            Simple_Point i3 = findIntersection(s3, s);
            intersections.push_back(i3);
        }

        if (doIntersect(s4, s)) {
            Simple_Point i4 = findIntersection(s4, s);
            intersections.push_back(i4);
        }
        if (intersections.size() == 0) {
            dist = 0.0;
        } else if (intersections.size() == 1) {
            Simple_Point p = isPointInsideRectangle(l2.start, l1.bbox_lines, l1.bbox_area) ? l2.start : l2.end;
            dist = euclideanDistance(p, intersections[0]);
        } else {
            dist = euclideanDistance(intersections[0], intersections[1]);
        }
    }

    return dist;
}

// 重合部分占比: l2在l1形成的矩形中的长度/l1的长度
double LDBSCAN::calc_box_coverages(Simple_Line l1, Simple_Line l2) {
    // 重合部分占比：重合/线长
    double dist = calc_overlap_area_dist(l1, l2);
    double coverage_ratio1 = dist / l2.line_length;
    double coverage_ratio2 = dist / l1.line_length;
    return coverage_ratio1 > coverage_ratio2 ? coverage_ratio1 : coverage_ratio2;
}

// 计算两直线夹角
double LDBSCAN::get_angle(double k1, double k2) {
    const bool is_k1_inf = std::isinf(k1);
    const bool is_k2_inf = std::isinf(k2);

    // 情况1：两条直线均为垂直线（斜率均为inf）
    if (is_k1_inf && is_k2_inf) {
        return 0.0; // 平行直线，夹角0°
    }
    // 情况2：只有k1为垂直线（斜率inf）
    else if (is_k1_inf) {
        return 90.0 - std::atan(std::fabs(k2)) * 180.0 / M_PI;
    }
    // 情况3：只有k2为垂直线（斜率inf）
    else if (is_k2_inf) {
        return 90.0 - std::atan(std::fabs(k1)) * 180.0 / M_PI;
    }
    // 情况4：两条直线均为普通斜率（有限值）
    else {
        return std::atan(std::fabs((k1 - k2) / (1 + k1 * k2))) * 180.0 / M_PI;
    }
}

// 寻找线i的epsilon邻域
std::vector<int> LDBSCAN::regionQuery(const std::vector<Simple_Line> &lines, int pointIdx) {
    std::vector<int> neighborPts;
    std::vector<int> neighborIds;
    Simple_Line anchor_line = lines[pointIdx];
    for (int i = 0; i < lines.size(); ++i) {
        if (i != pointIdx) {
            Simple_Line target_line = lines[i];
            double angle = get_angle(anchor_line.slope, target_line.slope);
            // 重合部分占比: l2在l1形成的矩形中的长度/l1的长度
            double coverage_ratio = calc_box_coverages(anchor_line, target_line);

            // 两线长度
            double len1 = anchor_line.line_length;
            double len2 = target_line.line_length;

            // 计算两个直线最近点之间间距
            double d1 = euclideanDistance(anchor_line.start, target_line.start);
            double d2 = euclideanDistance(anchor_line.start, target_line.end);
            double d3 = euclideanDistance(anchor_line.end, target_line.start);
            double d4 = euclideanDistance(anchor_line.end, target_line.end);
            double min_dist = std::min(std::min(d1, d2), std::min(d3, d4));
            double max_dist = std::max(std::max(d1, d2), std::max(d3, d4));

            // 同一条边的不同长度的线: 长边压住短边的情况
            if (angle < _epsilon_angle && coverage_ratio > 0) {
                neighborPts.push_back(i);
                neighborIds.push_back(lines[i].id);
                continue;
            }

            // 一条边被分为多段且不存在互压的情况
            if (angle < _epsilon_angle && coverage_ratio <= 0.0) {
                Simple_Point l1_s{anchor_line.sx, anchor_line.sy};
                Simple_Point l2_s{target_line.sx, target_line.sy};
                Simple_Point l1_e{anchor_line.ex, anchor_line.ey};
                Simple_Point l2_e{target_line.ex, target_line.ey};
                std::vector<std::pair<Simple_Point, Simple_Point>> pp{{l1_s, l2_s}, {l1_s, l2_e}, {l1_e, l2_s}, {l1_e, l2_e}};
                std::vector<double> dist_list;
                dist_list.push_back(euclideanDistance(l1_s, l2_s));
                dist_list.push_back(euclideanDistance(l1_s, l2_e));
                dist_list.push_back(euclideanDistance(l1_e, l2_s));
                dist_list.push_back(euclideanDistance(l1_e, l2_e));
                auto min_it = std::min_element(dist_list.begin(), dist_list.end());
                int min_index = std::distance(dist_list.begin(), min_it);

                auto max_it = std::max_element(dist_list.begin(), dist_list.end());
                int max_index = std::distance(dist_list.begin(), max_it);
                // 两线四个端点间的最小距离
                double min_dist = dist_list[min_index];

                Simple_Line max_dist_line{pp[min_index].first.x, pp[min_index].first.y, pp[min_index].second.x, pp[min_index].second.y, -1};

                double a1 = get_angle(max_dist_line.slope, anchor_line.slope);
                double a2 = get_angle(max_dist_line.slope, target_line.slope);
                double min_angle = std::min(a1, a2);
                if (min_angle < _epsilon_angle * 2 && min_dist < _r) {
                    neighborPts.push_back(i);
                    neighborIds.push_back(lines[i].id);
                    continue;
                }
            }

            // 角度较大, 且有一点交叉, 两个线长度差不多, 最近点之间间距较小
            if (angle < 45. && coverage_ratio > 0. && std::abs(len1 - len2) < _r * 0.5) {
                if (min_dist < _r * 0.5 && max_dist < _r * 2) {
                    neighborPts.push_back(i);
                    neighborIds.push_back(lines[i].id);
                    continue;
                }
            }
        }
    }

    return neighborPts;
}

void LDBSCAN::expandCluster(std::vector<Simple_Line> &lines, int pointIdx, int clusterIdx) {
    lines[pointIdx].cluster = clusterIdx;
    std::vector<int> neighborPts = regionQuery(lines, pointIdx);

    if (neighborPts.size() < _minPts)
        return;

    for (int i = 0; i < neighborPts.size(); ++i) {
        int nextPointIdx = neighborPts[i];
        if (!lines[nextPointIdx].visited) {
            lines[nextPointIdx].visited = true;
            expandCluster(lines, nextPointIdx, clusterIdx);
        }
        if (lines[nextPointIdx].cluster == -1) {
            // 未定义类别的点，将其标记为当前簇
            lines[nextPointIdx].cluster = clusterIdx;
        }
    }
}

void LDBSCAN::ldbscan(std::vector<Simple_Line> &lines) {
    int clusterIdx = 0;
    for (int i = 0; i < lines.size(); ++i) {
        if (!lines[i].visited) {
            lines[i].visited = true;
            std::vector<int> neighborPts = regionQuery(lines, i);
            expandCluster(lines, i, ++clusterIdx);
        }
    }
}

bool compareSegments(const Simple_Line &s1, const Simple_Line &s2) { return s1.line_length > s2.line_length; }

std::vector<Simple_Line> LDBSCAN::findLongestSegments(std::vector<Simple_Line> lines, std::vector<int> labels, std::vector<int> &longest_line_index) {
    std::unordered_map<int, std::vector<Simple_Line>> categoryMap;
    std::unordered_map<int, double> angleMap;
    for (int i = 0; i < lines.size(); ++i) {
        categoryMap[labels[i]].push_back(lines[i]);
    }
    std::vector<Simple_Line> longestLines;
    for (auto &entry : categoryMap) {
        std::vector<Simple_Line> cur_class_lines = entry.second;
        // 对每个类别中的线段按长度从大到小排序
        sort(cur_class_lines.begin(), cur_class_lines.end(), compareSegments);
        // 将最长的线段放入结果集合
        longestLines.push_back(cur_class_lines[0]);
        longest_line_index.push_back(cur_class_lines[0].id);
    }
    return longestLines;
}

std::vector<Simple_Line> LDBSCAN::findBestSegments(std::vector<Simple_Line> lines, std::vector<int> &bestLines_index, Plane_Obj_Vectorization loopBoj) {
    std::unordered_map<int, std::vector<Simple_Line>> categoryMap;
    std::unordered_map<int, double> angleMap;
    for (int i = 0; i < lines.size(); ++i) {
        if (lines[i].line_length == 0.0) {
            std::cout << "Error: line length is 0" << std::endl;
            // continue;
        }
        std::vector<double> anglelist;
        // 第i类别的放在一个索引的vector中LoopArrow中的左右邻居
        categoryMap[lines[i].cluster].push_back(lines[i]);
        // 通过当前线的类别找到loop的该类最优线id
        int cls_id_in_loop = loopBoj._cls_id[lines[i].cluster];
        if (!loopBoj._joints.count(cls_id_in_loop)) {
            std::cout << "Error: joints not found for cls_id_in_loop: " << cls_id_in_loop << std::endl;
            continue;
        }
        int next_line_id = loopBoj._joints[cls_id_in_loop].first;
        if (!loopBoj._prv_id.count(cls_id_in_loop)) {
            std::cout << "Error: joints not found for cls_id_in_loop: " << cls_id_in_loop << std::endl;
            continue;
        }
        int prv_line_id = loopBoj._prv_id[cls_id_in_loop];
        double max_angle = 0.0;
        for (int j = 0; j < lines.size(); ++j) {
            // 当前线和两个邻居的夹角
            if (lines[j].id == next_line_id || lines[j].id == prv_line_id) {
                anglelist.push_back(get_angle(lines[i].slope, lines[j].slope));
            }
        }
        double sum = 0.0;
        sum += anglelist[0];
        sum += anglelist[1];
        angleMap[lines[i].id] = sum / 2.0;
    }
    std::vector<Simple_Line> bestLines;
    for (auto &entry : categoryMap) {
        std::vector<Simple_Line> cur_class_lines = entry.second;
        // 对每个类别中的线段按长度排序
        sort(cur_class_lines.begin(), cur_class_lines.end(), compareSegments);
        // 以每类的最长线为单位1
        double max_length = cur_class_lines[0].line_length + 1e-10;
        // 记录得分最高的线分数和索引
        double max_score = -100;
        double max_score_line_index = -1;
        for (int i = 0; i < cur_class_lines.size(); ++i) {
            double cur_angle = 0.0;
            // 当前线的角度分数和长度分数
            // FIXME: 这里角度取和其他线的最大值是不正确的，却复杂越容易出错，直行箭头应该可以，这里的本意是这条线应该和相邻的两条线夹角最大
            if (angleMap.count(cur_class_lines[i].id)) {
                cur_angle = angleMap[cur_class_lines[i].id];
            }
            double angle_score = cur_angle / 90.0;
            double length_score = cur_class_lines[i].line_length / max_length;
            double cur_line_score = angle_score * _angle_weight + length_score * _length_weight;
            if (std::isnan(cur_line_score)) {
                std::cout << "Error: cur_line_score is nan" << std::endl;
            }
            // 记录最高分数
            if (cur_line_score >= max_score) {
                max_score = cur_line_score;
                max_score_line_index = i;
            }
        }

        // 将最长的线段放入结果集合
        bestLines.push_back(cur_class_lines[max_score_line_index]);
        bestLines_index.push_back(cur_class_lines[max_score_line_index].id);
    }
    return bestLines;
}
