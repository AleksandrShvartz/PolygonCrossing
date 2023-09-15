#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
struct Line {
    int a;
    int b;
    int c;
    Line():a(0),b(0),c(0) {

    }
    Line(int a, int b, int c) :a(a), b(b), c(c) {
    }
    int sign(std::pair<int, int> p) {
        long long lx = p.first;
        long long ly = p.second;
        long long v = a * lx + b * ly + c;
        return  v > 0L ? 1 : v < 0L ? -1 : 0;
    }
    long long dist(std::pair<int, int> p) {
        long long lx = p.first;
        long long ly = p.second;
        long long v = a * lx + b * ly + c;
        return v;
    }
};
struct Edge {
    Edge(){}
    Edge(std::pair<int, int> a, std::pair<int, int> b):a(a),b(b){}
    std::pair<int, int> a;
    std::pair<int, int> b;

    std::pair<int,int> findIntersection(Line& l, std::vector<int>& answer) {
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
               answer.push_back(a.first);
               answer.push_back(a.second);
               answer.push_back(b.first);
               answer.push_back(b.second);
               answer.push_back(b.first);
               answer.push_back(b.second);
               return a;
            }
           
            else {
                throw "Error: there should be an intersection point";
            }
        }
        int x = - det_x / det;
        int y = - det_y / det;

        answer.push_back(x);
        answer.push_back(y);
        return {x,y};
    }
};

struct ShiftedArr {
    std::vector<std::pair<int, int>>& arr;
    int shift;
    ShiftedArr(std::vector<std::pair<int, int>>& a, int sh) :arr(a), shift(sh) {}
    // return shifted element.
    std::pair<int, int>& operator[](int i) {
        return arr[((long long)i + shift) % arr.size()];
    }
    int ternaryMinSearch(int l, int r, Line& line, int sign) {
        int l1, l2, sign_l1,sign_l2;
        while (l < r) {
            l1 = l + (r - l) / 3;
            l2 = r - (r - l) / 3;
            sign_l1 = line.sign((*this)[l1]) == sign ? 1 : -1;
            sign_l2 = line.sign((*this)[l2]) == sign ? 1 : -1;
            if ((abs(line.dist((*this)[l1])) * sign_l1) == (abs(line.dist((*this)[l2])) * sign_l2)) {
                l1++;
                sign_l1 = line.sign((*this)[l1]) == sign ? 1 : -1;
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
    std::pair<int, int> binSearch(int l, int r, Line& line) {
        int mid = (l + r) / 2;
        while (l < r - 1) {
            mid = (l + r) / 2;
            if (line.sign((*this)[l]) != line.sign((*this)[mid]))
                r = mid;
            else
                l = mid;
        }
        return { l,r };
    }

};
/*
*  Return vector of all intersection points,
*   if line consists with edge return vector with size > 4 
*/
std::vector<int> solve(
    std::vector<std::pair<int, int>>& verts, 
    Line line,
    std::pair<int,int>& min,
    std::pair<int, int>& max
    ) {

    ShiftedArr poly(verts, min.first);
    std::vector<int> answer;
    std::vector<std::pair<int, std::pair<int,int>>> answer_edges;
    // check if line || Ox
    if (line.a == 0) {
       //todo
    }
    //check if line || Oy
    if (line.b == 0) {
        //todo 
    }

    int l = 0;
    int r = poly.arr.size();
    bool is_top = false;
    std::cout << "Init Left point: " << poly[l].first << " " << poly[l].second <<"index: "<<l<< "\n";
    std::cout << "Init Right point: " << poly[r].first << " " << poly[r].second << "index: " << r << "\n";
        //deside which part of polygon we should choose 
        long long value_l = -(line.a * (long long) poly[min.first].first + line.c); // -(a*x_0+c)
        long long value_r =   line.b * (long long) poly[min.first].second; //   y*b
        if (((value_l < value_r) && line.b > 0) || ((value_l > value_r) && line.b < 0)) { // bottom part
            l = max.first;
           // std::cout << poly[1].second << "\n";
            std::cout <<"Bottom" << "\n";
        }
        else {  // top part
            r = max.first;
            if (r == 0)
                r = poly.arr.size();
            is_top = true;
            std::cout << "Top" << "\n";
        }
        //ternary search to find middle point (which devide to parts with one answer maximum) 
        // it is a minimum by siged distance where sign is "+" if point on the [l] and [r] side and "-" if opposite 
        auto v = poly.ternaryMinSearch(l, r, line, line.sign(poly[l]));
        std::cout << "First point: " << poly[0].first << " " << poly[0].second << "\n";
        std::cout << "Ternary point: " <<poly[v].first<<" "<<poly[v].second << "index: " << v << "\n";
        std::cout << " Left point: " << poly[l].first << " " << poly[l].second << "index: " << l << "\n";
        std::cout << " Right point: " << poly[r].first << " " << poly[r].second << "index: " << r << "\n";
        //two or one binary searches to find answer
        int sign_l = line.sign(poly[l]);
        int sign_r = line.sign(poly[r]);
        int mid  = v; // todo check if l==r in search
        int sign_mid = line.sign(poly[mid]);
        //  if intersection point exists
        auto find_answer = [&](int left, int right) {
            auto edge = poly.binSearch(left, right, line);
            Edge segment(poly[edge.first], poly[edge.second]);
            int size = answer.size();
            auto point = segment.findIntersection(line, answer);
            if (point == poly[edge.first]) {
                answer_edges.push_back(std::make_pair(edge.first, point));
            }
            if (point == poly[edge.second]) {
                answer_edges.push_back(std::make_pair(edge.second, point));
            }
        };
        if (sign_l * sign_mid != 1) {
            find_answer(l, mid);
        }
        if (sign_r * sign_mid != 1) {
            find_answer(mid, r);
        }
      
        if (answer.size() == 2 && sign_l *sign_r != 1) { // can be another point
            if (is_top) {
               //check bottom 
                find_answer(max.first, poly.arr.size() - 1);         
            }
            else {
                find_answer(0, max.first);             
            }
        }
        //check parallel
        if (answer_edges.size() == 2 && abs(answer_edges[0].first - answer_edges[1].first) == 1) {
            answer.push_back(answer_edges[0].second.first);
            answer.push_back(answer_edges[0].second.second);
        }
    return answer;
}
int main() {
    std::string line;  
    std::ifstream in("input.txt"); 
    std::ofstream out("out.txt");
    if (in.is_open())
    {
        int n = -1;
        in >> n;
        std::pair<int, int> min_x;
        std::pair<int, int> max_x;
        std::vector<std::pair<int, int>> verts;
        verts.reserve(n);
        for (int i = 0; i != n; ++i) {
            int x, y;
            in >> x >> y;
            verts.emplace_back(std::make_pair(x, y));
            if (i) {
                max_x = max_x.second >= x ? max_x : std::make_pair(i, x);
                min_x = min_x.second <= x ? min_x : std::make_pair(i, x);
            }
            else {
                max_x = std::make_pair(0, x);
                min_x = std::make_pair(0, x);
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
            lines.emplace_back(Line(a,b,c));
        }
        if (out.is_open()) {
            for (auto line : lines) {
                auto res = solve(verts, line, min_x, max_x);

                if (res.size() > 4) {
                    out << -1 << "\n";
                    out << res[0] << ' ' << res[1] << ' ';
                    for (int i = 2; i < res.size(); i+=2) {
                        if (res[0] != res[i] || res[1] != res[i + 1]) {
                            out << res[i] << ' ' << res[i+1] << ' ';
                            break;
                        }
                    }
                    out << '\n';
                }
                else {
                       if (res.size() == 4 && res[0] == res[2] && res[1] == res[3]) {
                            out << 1 << '\n';
                            out << res[0] << ' ' << res[1] << '\n';
                            continue;
                        }
                    out << res.size() / 2 << '\n';
                    for (auto coord : res) {
                        out << coord << ' ';
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