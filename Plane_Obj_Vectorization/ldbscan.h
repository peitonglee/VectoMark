#ifndef LDBSCAN_H
#define LDBSCAN_H

#include "simple_data.h"
#include "Plane_Obj_Vectorization.h"
#include <vector>

class LDBSCAN {
  public:
    LDBSCAN() {}
    ~LDBSCAN() {}

    double euclideanDistance(Simple_Point p1, Simple_Point p2);

    Simple_Point findIntersection(Simple_Segment l1, Simple_Segment l2);

    bool isPointInsideRectangle(Simple_Point p, std::vector<Simple_Point> rect, double area);

    bool doIntersect(Simple_Segment l1, Simple_Segment l2);

    double calc_overlap_area_dist(Simple_Line l1, Simple_Line l2);

    double calc_box_coverages(Simple_Line l1, Simple_Line l2);

    std::vector<int> regionQuery(const std::vector<Simple_Line> &lines, int pointIdx);

    void expandCluster(std::vector<Simple_Line> &lines, int pointIdx, int clusterIdx);

    void ldbscan(std::vector<Simple_Line> &lines);

    void set_minimum_length(double min_length) {
        if (min_length < 0.15 && min_length > 0.2) {
            // 根据经验值设定
            min_length = 0.15;
        }
        _r = min_length;
    }
    double get_minimum_length() { return _r; }

    std::vector<Simple_Line> findLongestSegments(std::vector<Simple_Line> lines, std::vector<int> labels, std::vector<int> &longest_line_index);

    std::vector<Simple_Line> findBestSegments(std::vector<Simple_Line> lines, std::vector<int> &bestLines_index, Plane_Obj_Vectorization loopBoj);

    double get_angle(double k1, double k2);
    double _box_r_scale_ratio = 0.6;

  private:
    double _epsilon_clover = 0.33;
    double _epsilon_angle = 6;
    int _minPts = 0;
    double _r = 0.0;
    double _angle_weight = 0.9;
    double _length_weight = 0.1;
};

#endif // LDBSCAN_H