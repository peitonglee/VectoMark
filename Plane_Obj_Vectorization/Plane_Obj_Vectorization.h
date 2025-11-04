#ifndef LOOP_ARROW_H
#define LOOP_ARROW_H
#include "simple_data.h"

using uint64 = unsigned long long;

// 状态数据结构
struct State {
    double cost;           // 累计成本
    int current;           // 当前节点索引
    uint64 mask;           // 访问掩码（64位）
    std::vector<int> path; // 路径记录

    bool operator<(const State &other) const {
        return cost > other.cost; // 最小堆
    }
};

// 状态记录（使用复合键）
struct StateKey {
    int current;
    uint64 mask;

    bool operator==(const StateKey &other) const { return current == other.current && mask == other.mask; }
};

struct Hash {
    size_t operator()(const StateKey &k) const { return std::hash<int>()(k.current) ^ std::hash<uint64>()(k.mask); }
};

class Plane_Obj_Vectorization {
  public:
    Plane_Obj_Vectorization() {}
    ~Plane_Obj_Vectorization() {}

    void loop_arrow(std::vector<Simple_Line> processedLines);
    void get_intersections(std::unordered_map<int, Simple_Line> lines);
    std::tuple<std::vector<Simple_Point>, std::unordered_map<int, int>, std::vector<int>> process_lines(const std::vector<Simple_Line> &lines);
    std::pair<std::vector<int>, double> find_optimal_path(const std::vector<Simple_Point> &points, const std::unordered_map<int, int> &pair_map);
    std::pair<std::vector<int>, double> AStar(const std::vector<Simple_Point> &points, const std::unordered_map<int, int> &pair_map);
    std::pair<std::vector<int>, double> find_approximate_path(const std::vector<Simple_Point> &points, const std::unordered_map<int, int> &pair_map);
    void update_endPoint(std::unordered_map<int, Simple_Line> lines, int cur_id, int target_id, Simple_Point &inter, int &mask, int type);
    void set_A_star();
    std::map<int, std::vector<double>> _m_line_coord_points;
    std::vector<std::pair<int, int>> _parallel_line_pairs;
    std::vector<std::pair<int, int>> _vertical_line_pairs;
    // 有环的链表
    std::unordered_map<int, std::pair<int, std::pair<Simple_Point, Simple_Point>>> _joints;
    std::vector<std::pair<int, std::vector<double>>> _figure_cycle_info;
    std::unordered_map<int, int> _prv_id;
    std::unordered_map<int, int> _cls_id;
    bool _loop_state = false;

  private:
    double _mean_length;
    double _min_length;
    double _max_length;
    bool _A_star_flag = false;
};
#endif // LOOP_ARROW_H