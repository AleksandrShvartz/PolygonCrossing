#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

struct Line {
    int a;
    int b;
    int c;
    Line() :a(0), b(0), c(0) {

    }
    Line(int a, int b, int c) :a(a), b(b), c(c) {
    }
    int sign(std::pair<int, int> p) {
        auto value = dist(p);
        return  (value > 0L) ? 1 : ((value < 0L) ? -1 : 0);
    }
    long long dist(std::pair<int, int> p) {
        long long lx = p.first;
        long long ly = p.second;
        long long v = a * lx + b * ly + c;
        return v;
    }
};


struct Edge {
    Edge() {}
    Edge(std::pair<int, int> a, std::pair<int, int> b) :a(a), b(b) {}
    std::pair<int, int> a;
    std::pair<int, int> b;

    std::pair<double, double> findIntersection(Line& l, std::vector<std::pair<double, double>>& answer) {
        long long A = (long long)b.second - a.second; // A = y2 - y1.
        long long B = (long long)a.first - b.first;  // B = x1 - x2.
        // C = y1 * (x2 - x1) - (y2 - y1) * x1
        long long C = a.second * ((long long)b.first - a.first) - ((long long)b.second - a.second) * a.first;
        // |A1 B1|
        // |A2 B2|
        long long det = l.a * B - A * l.b;
        // |A1 C1|
        // |A2 C2|
        long long det_y = l.a * C - A * l.c;
        // |C1 B1|
        // |C2 B2|
        long long det_x = l.c * B - l.b * C;
        //lines are parallel 
        if (det == 0L) {
            // check if lines matches     
            if (det_x == 0L && det_y == 0L) {
                answer.push_back(a);
                answer.push_back(b);
                answer.push_back(std::make_pair(INT_MAX, INT_MAX));   // trash
                return a;
            }

            else {
                throw "Error: there should be an intersection point";
            }
        }
        double x = -(double)det_x / det;
        double y = -(double)det_y / det;

        answer.push_back({ x,y });
        return { x,y };
    }
};

struct ShiftedArr {
    std::vector<std::pair<int, int>>& arr;
    int shift;
    ShiftedArr(std::vector<std::pair<int, int>>& a, int sh) :arr(a), shift(sh) {}
    // return shifted element.
    std::pair<int, int>& operator[](int i) {
        return arr[(i + arr.size() + shift) % arr.size()];
    }
    int ternaryMinSearch(int l, int r, Line& line) {
        // near needed point function value is negative
        int sign = line.sign((*this)[l]);
        int l1, l2, sign_l1, sign_l2;
        while (l < r) {
            l1 = l + (r - l) / 3;
            l2 = r - (r - l) / 3;
            sign_l1 = line.sign((*this)[l1]) == sign ? 1 : -1;
            sign_l2 = line.sign((*this)[l2]) == sign ? 1 : -1;
            if (abs(line.dist((*this)[l1])) * sign_l1 == abs(line.dist((*this)[l2])) * sign_l2) {
                l1++;
                sign_l1 = line.sign((*this)[l1]) == sign ? 1 : -1;
            }
            if (l1 == l2) {
                // if points are neigbhours and with same distance
                int min_ind = l;
                int sign_l = line.sign((*this)[l]) == sign ? 1 : -1;
                auto min = abs(line.dist((*this)[l])) * sign_l;
                for (int i = l; i <= r; ++i) {
                    int sign_i = line.sign((*this)[i]) == sign ? 1 : -1;
                    if ((abs(line.dist((*this)[i])) * sign_i) < min) {
                        min = (abs(line.dist((*this)[i])) * sign_i);
                        min_ind = i;
                    }
                }
                return min_ind;
            }
            if ((abs(line.dist((*this)[l1])) * sign_l1) > (abs(line.dist((*this)[l2])) * sign_l2)) {
                l = l1 + 1;
            }
            else {
                r = l2 - 1;
            }
        }

        return l;
    }
    int binSearch(int l, int r, Line& line) {
        int mid = (l + r) / 2;
        while (l < r - 1) {
            mid = (l + r) / 2;
            if (line.sign((*this)[l]) != line.sign((*this)[mid]))
                r = mid;
            else
                l = mid;
        }
        return l;
    }

};
/*
*  Return vector of all intersection points,
*   if line consists with edge return vector with size > 4
*/
std::vector<std::pair<double, double>> solve(
    std::vector<std::pair<int, int>>& verts,
    Line line,
    int min,
    int max,
    int min_y_ind,
    int max_y_ind
) {
    int shift = line.b == 0 ? max_y_ind : min;
    ShiftedArr poly(verts, shift);
    std::vector<std::pair<double, double>> answer;
    std::vector<int> answer_edges;

    int l = 0;
    int r = poly.arr.size();
    bool is_first_half = false;
    // check if line || Oy
    if (line.b == 0) {
        long long value = (line.a * (long long)poly[0].first + line.c); // (a*x_1+c)
        auto sign = [](long long x) {return x > 0 ? 1 : x < 0 ? -1 : 0; };
        if (sign(value) * sign(line.a) != 1) {             // right part
            r = min_y_ind;
            if (r <= 0)
                r = poly.arr.size();
            is_first_half = true;
        }
        else {//left part
            l =  min_y_ind - max_y_ind;
        }
    }
    else {
        // std::cout << "Init Left point: " << poly[l].first << " " << poly[l].second <<"index: "<<l<< "\n";
        //std::cout << "Init Right point: " << poly[r].first << " " << poly[r].second << "index: " << r << "\n";
        //deside which part of polygon we should choose

        long long value_l = -(line.a * (long long)poly[0].first + line.c); // -(a*x_0+c)
        long long value_r = line.b * (long long)poly[0].second; //   y*b

        if (((value_l < value_r) && line.b > 0) || ((value_l > value_r) && line.b < 0)) { // bottom part
            l = max - min;
            // std::cout << poly[1].second << "\n";
            //  std::cout <<"Bottom" << "\n";
        }
        else {  // top part
            r = max - min;
            if (r <= 0)
                r = poly.arr.size();
            is_first_half = true;
            //std::cout << "Top" << "\n";
        }
    }
    //  if intersection point exists
    auto find_answer = [&](int left, int right) {
        int n = poly.arr.size();
        left = (left % n + n) % n;
        right = (right % n + n) % n;
        if (left > right) {
            right += n;
        }
        if (line.sign(poly[left]) * line.sign(poly[right]) == 1) {

            return;
        }
        int v = poly.binSearch(left, right, line);
        Edge edge(poly[v], poly[v + 1]);
        auto point = edge.findIntersection(line, answer);
        //check if point is vertex
        if (floor(point.first) == point.first && floor(point.second) == point.second) {
            std::pair<int, int> int_point = std::make_pair(point.first + 0.5 - (point.first < 0),
                point.second + 0.5 - (point.second < 0));
            if (int_point == poly[v]) {
                answer_edges.push_back(v);
            }
            if (int_point == poly[v + 1]) {
                answer_edges.push_back(v + 1);
            }
        }

    };


    //ternary search to find middle point (which devide to parts with one answer maximum) 
    // it is a minimum by siged distance where sign is "+" if point on the [l] and [r] side and "-" if opposite 
    if (line.sign(poly[l]) == 0) {
        find_answer(l, l);
        //  if (line.sign(poly[left]) * line.sign(poly[right]) == 1)
        find_answer(l + 1, l - 1 + poly.arr.size());

    }
    else {
        int minimum_of_the_fun = poly.ternaryMinSearch(l, r, line);
        //std::cout << "First point: " << poly[0].first << " " << poly[0].second << "\n";
        // std::cout << "Ternary point: " <<poly[].first<<" "<<poly[v].second << "index: " << v << "\n";
        //std::cout << " Left point: " << poly[l].first << " " << poly[l].second << "index: " << l << "\n";
        //std::cout << " Right point: " << poly[r].first << " " << poly[r].second << "index: " << r << "\n";
        //two or one binary searches to find answer
        int sign_l = line.sign(poly[l]);
        int sign_r = line.sign(poly[r]);
        int mid = minimum_of_the_fun;
        int sign_mid = line.sign(poly[mid]);

        if (sign_mid == 0) {
            find_answer(mid, mid);
            // if(line.sign(poly[mid-1])* line.sign(poly[mid + 1]) != 1)
            find_answer(mid + 1, mid - 1 + poly.arr.size());
        }
        /*else {
            if (sign_l * sign_mid != 1) {
                find_answer(l, mid);
            }
            if (sign_r * sign_mid != 1) {
                find_answer(mid, r);
            }
        }
        if (sign_l * sign_mid != 1) {
            find_answer(l, mid);
            find_answer(mid + 1, l - 1);
        }
        else if (sign_r * sign_mid != 1) {
                find_answer(mid, r);
                find_answer(r+1, mid-1);}*/
        if (sign_l * sign_mid != 1) {
            find_answer(l, mid);
        }
        if (sign_r * sign_mid != 1) {
            find_answer(mid, r);
        }
        if (answer.size() == 1 && sign_l * sign_r != 1) { // can be another point
 /*           int ind = (line.b == 0 ? min_y_ind : max - min);*/
            if (is_first_half) {
                //check bottom 
                find_answer(r, l);
            }
            else {
                find_answer(0, l);
            }
        }
    }
    //check parallel   
    if (answer_edges.size() == 2 && (abs((answer_edges[0] - answer_edges[1])) == 1 ||
        abs((answer_edges[0] - answer_edges[1])) == poly.arr.size() - 1)) {
        answer.push_back(poly[answer_edges[0]]);
        answer.push_back(poly[answer_edges[1]]);
        answer.push_back(std::make_pair(INT_MAX, INT_MAX));
    }
    return answer;
}

int main() {
    std::string line;
    std::ifstream in("input.txt");
    std::ofstream out("output.txt");
    if (in.is_open())
    {
        int n = -1;
        in >> n;
        std::pair<std::pair<int, int>, int> min_x;
        std::pair<std::pair<int, int>, int> max_x;
        std::pair<std::pair<int, int>, int> min_y;
        std::pair<std::pair<int, int>, int> max_y;
        std::vector<std::pair<int, int>> verts;
        verts.reserve(n);
        for (int i = 0; i != n; ++i) {
            int x, y;
            in >> x >> y;
            verts.emplace_back(std::make_pair(x, y));
            if (i) {
                min_x = min(min_x, { {x,y},i });
                max_x = max(max_x, { {x,y},i });
                min_y = (min_y.first.second > y) ? std::make_pair(std::make_pair(x, y), i) : min_y;
                max_y = (max_y.first.second < y) ? std::make_pair(std::make_pair(x, y), i) : max_y;
            }
            else {
                min_x = { {x,y},i };
                max_x = { {x,y},i };
                min_y = { {x,y},i };
                max_y = { {x,y},i };
            }
        }

        int m = -1;
        in >> m;
        std::vector<Line> lines;
        lines.reserve(m);
        for (int i = 0; i != m; ++i)
        {
            int a, b, c;
            in >> a >> b >> c;
            lines.emplace_back(Line(a, b, c));
        }
        if (out.is_open()) {
            for (auto line : lines) {
                auto res = solve(verts, line, min_x.second, max_x.second, min_y.second, max_y.second);             // this is answer
                res.erase(std::unique(res.begin(), res.end()), res.end());
                if (res.size() > 2) {
                    // Infinity
                    out << -1 << "\n";
                    out << res[0].first << ' ' << res[0].second << ' ';
                    out << res[1].first << ' ' << res[1].second;

                    out << '\n';
                }
                else {
                    out << res.size() << '\n';
                    for (auto coord : res) {
                        out << coord.first << ' ' << coord.second << ' ';
                    }
                    out << '\n';
                }
            }
            out.close();
        }
        else {
            std::cout << "Error with output file\n";
        }
        in.close();
    }
    else {
        std::cout << "Error with input file\n";
    }

    return 0;
}
