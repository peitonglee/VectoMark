#ifndef SIMPLE_DATA_H
#define SIMPLE_DATA_H
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

struct Simple_Point {
    double x, y, id;
    Simple_Point() : x(0), y(0), id(0) {}
    Simple_Point(double _x, double _y) : x(_x), y(_y), id(0) {}
    Simple_Point(double _x, double _y, double _id) : x(_x), y(_y), id(_id) {}
};

struct Simple_Segment {
    Simple_Point start, end;
};

struct Simple_Line {
    // 端点坐标
    double sx, sy, ex, ey, mx, my;
    Simple_Point start, end;
    double line_length;
    // 类别
    int cluster;
    // 访问标记
    bool visited;
    // 外接矩形的宽
    double _r;
    // 斜率、x轴的夹角的角度
    double slope;
    double angle;
    // 和x轴的夹角的sin和cos值
    double sinValue;
    double cosValue;
    // bbox: 最小外接矩形的四个顶点坐标
    std::vector<double> bbox;
    std::vector<Simple_Point> bbox_lines;
    double bbox_area;

    int id;

    Simple_Line() : sx(0), sy(0), ex(0), ey(0), cluster(-1), visited(false) {}
    Simple_Line(double _sx, double _sy, double _ex, double _ey, int _id) : sx(_sx), sy(_sy), ex(_ex), ey(_ey), cluster(-1), visited(false), id(_id) {
        if (_sx > _ex) {
            double temp = sx;
            sx = ex;
            ex = temp;
            temp = sy;
            sy = ey;
            ey = temp;
        }
        start = {sx, sy};
        end = {ex, ey};
        mx = (sx + ex) / 2.0;
        my = (sy + ey) / 2.0;
        line_length = sqrt(pow(start.x - end.x, 2) + pow(start.y - end.y, 2));
        slope = calculateSlope(sx, sy, ex, ey);
        angle = calculateAngle(slope);
        sinValue = sin(angle * M_PI / 180.0);
        cosValue = cos(angle * M_PI / 180.0);
    }

    void calc_box(double r) {
        _r = r;
        bbox = std::vector<double>{
            sx - r * sinValue, sy + r * cosValue, sx + r * sinValue, sy - r * cosValue, ex - r * sinValue, ey + r * cosValue, ex + r * sinValue, ey - r * cosValue,
        };
        bbox_lines.push_back({sx - r * sinValue, sy + r * cosValue}); // 左上
        bbox_lines.push_back({sx + r * sinValue, sy - r * cosValue}); // 左下
        bbox_lines.push_back({ex + r * sinValue, ey - r * cosValue}); // 右下
        bbox_lines.push_back({ex - r * sinValue, ey + r * cosValue}); // 右上
        bbox_area = 2 * r * line_length;
    }

    double calculateSlope(double x1, double y1, double x2, double y2) { return (y2 - y1) / (x2 - x1); }

    double calculateAngle(double slope) { return atan(slope) * 180.0 / M_PI; }
};

struct Line_Creater {
    std::map<int, std::vector<double>> lines; // sx, sy, ex, ey, id
    double minimum_length = 100;              // 最小线段长度

    std::vector<Simple_Line> gen_lines(float scale_ratio) {
        std::vector<Simple_Line> result;

        for (auto &line : lines) {
            auto id = line.first;
            auto info = line.second;
            Simple_Line sl(info[0], info[1], info[2], info[3], id);
            result.push_back(sl);
            if (sl.line_length < minimum_length) {
                minimum_length = sl.line_length;
            }
        }
        // 这里本来是设置最短线的长度*scale_ratio, 但是最短线可能太短以至于相近的线无法聚类, 这里根据经验设置为0.15
        minimum_length = 0.15;
        for (auto &line : result) {
            line.calc_box(minimum_length * scale_ratio);
        }
        return result;
    }
};
#endif
