#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <random>
#include <complex>
struct Line {
    int a;
    int b;
    int c;
    Line() :a(0), b(0), c(0) {

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
// for testing
//struct pt {
//    pt(int x, int y) : x(x), y(y) {}
//    int x, y;
//};
//
//bool cmp(pt a, pt b) {
//    return a.x < b.x || a.x == b.x && a.y < b.y;
//}
//
//bool cw(pt a, pt b, pt c) {
//    return a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y) < 0;
//}
//
//bool ccw(pt a, pt b, pt c) {
//    return a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y) > 0;
//}
//
//void convex_hull(std::vector<pt>& a) {
//    if (a.size() == 1)  return;
//    sort(a.begin(), a.end(), &cmp);
//    pt p1 = a[0], p2 = a.back();
//    std::vector<pt> up, down;
//    up.push_back(p1);
//    down.push_back(p1);
//    for (size_t i = 1; i < a.size(); ++i) {
//        if (i == a.size() - 1 || cw(p1, a[i], p2)) {
//            while (up.size() >= 2 && !cw(up[up.size() - 2], up[up.size() - 1], a[i]))
//                up.pop_back();
//            up.push_back(a[i]);
//        }
//        if (i == a.size() - 1 || ccw(p1, a[i], p2)) {
//            while (down.size() >= 2 && !ccw(down[down.size() - 2], down[down.size() - 1], a[i]))
//                down.pop_back();
//            down.push_back(a[i]);
//        }
//    }
//    a.clear();
//    for (size_t i = 0; i < up.size(); ++i)
//        a.push_back(up[i]);
//    for (size_t i = down.size() - 2; i > 0; --i)
//        a.push_back(down[i]);
//}
//
//void run(Line& line, std::vector<pt> points, std::string out_name);
//
//bool test(int number) {
//    int n = rand() % 100;
//    std::vector<pt> points;
// 
//    std::vector<std::pair<int, int>> poly;
//    for (int i = 0; i < n; ++i) {
//        points.push_back(pt(rand()%1000, rand()%1000));
//       
//    }
//    convex_hull(points);
//    n = points.size();
//    const int k = 2;     // doesnot
//    for (int i = 0; i < n; ++i) {
//        points[i].x *= k;
//        points[i].y *= k;
//        poly.push_back({ points[i].x, points[i].y });
//    }
//    bool res = true;
//    // test 2 and infinity intersections 
//    for (int i = 0; i < n; ++i)
//    {
//        for (int j = i + 1; j < n - 1; ++j)
//        {
//            int t1 = 1;
//            int t2 = 1 ;
//            int x1 = (points[i].x + points[i + 1].x) / k * t1;
//            int y1 = (points[i].y + points[i + 1].y) / k * t1;
//            int x2 = (points[j].x + points[j + 1].x) / k * t2;
//            int y2 = (points[j].y + points[j + 1].y) / k * t2;
//            int a = y1 - y2;
//            int b = x2 - x1;
//            int c = -x1*a-y1*b;
//            Line line(a, b, c);
//            auto name = std::string("out") + std::to_string(number)+"_line"+std::to_string(i)+std::to_string(j) + ".txt";
//            std::ofstream inp(std::string("input") + std::to_string(number) + "_line" + std::to_string(i) + std::to_string(j)+ ".txt");
//            if (inp.is_open()) {
//                inp << n << '\n';
//                for (int i = 0; i < n; ++i) {
//                    inp << points[i].x << ' ' << points[i].y<<'\n';
//                }
//                inp << 1 <<'\n';
//                inp << line.a << ' ' << line.b << ' ' << line.c<<'\n';
//                inp.close();
//            }
//            run(line , points,name);
//            std::ifstream in(name);
//            int num = -5,x_1,x_2,y_1,y_2;
//            try {
//                in >> num;
//                in >> x_1 >> y_1 >> x_2 >> y_2;
//                res&= x_1 == x1 && y_1 == y1 && x_2 == x2 && y_2 == y2 || 
//                    x_1 == x2 && y_1 == y2 && x_2 == x1 && y_2 == y1;
//            }
//            catch (std::exception& e) {
//                return false;
//            }
//
//        }
//    }
//    for (int i = 0; i < n; ++i)
//    {
//        for (int j = i + 1; j < n - 1; ++j)
//        {
//            int t1 = 1;
//            int t2 = 1;
//            int x1 = (points[i].x + points[i + 1].x) / k * t1;
//            int y1 = (points[i].y + points[i + 1].y) / k * t1;
//            int x2 = (points[j].x + points[j + 1].x) / k * t2;
//            int y2 = (points[j].y + points[j + 1].y) / k * t2;
//            int a = y1 - y2;
//            int b = x2 - x1;
//            int c = -x1 * a - y1 * b;
//            Line line(a, b, c);
//            auto name = std::string("out") + std::to_string(number) + "_line" + std::to_string(i) + std::to_string(j) + ".txt";
//            std::ofstream inp(std::string("input") + std::to_string(number) + "_line" + std::to_string(i) + std::to_string(j) + ".txt");
//            if (inp.is_open()) {
//                inp << n << '\n';
//                for (int i = 0; i < n; ++i) {
//                    inp << points[i].x << ' ' << points[i].y << '\n';
//                }
//                inp << 1 << '\n';
//                inp << line.a << ' ' << line.b << ' ' << line.c << '\n';
//                inp.close();
//            }
//            run(line, points, name);
//            std::ifstream in(name);
//            int num = -5, x_1, x_2, y_1, y_2;
//            try {
//                in >> num;
//                in >> x_1 >> y_1 >> x_2 >> y_2;
//                res &= x_1 == x1 && y_1 == y1 && x_2 == x2 && y_2 == y2 ||
//                    x_1 == x2 && y_1 == y2 && x_2 == x1 && y_2 == y1;
//            }
//            catch (std::exception& e) {
//                return false;
//            }
//
//        }
//    }
//
//    return res;
//}


struct Edge {
    Edge(){}
    Edge(std::pair<int, int> a, std::pair<int, int> b):a(a),b(b){}
    std::pair<int, int> a;
    std::pair<int, int> b;

    std::pair<double,double> findIntersection(Line& l, std::vector<std::pair<double, double>>& answer) {
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
        double x = - (double)det_x / det;
        double y = - (double)det_y / det;

        answer.push_back({ x,y });
        return {x,y};
    }
};

struct ShiftedArr {
    std::vector<std::pair<int, int>>& arr;
    int shift;
    ShiftedArr(std::vector<std::pair<int, int>>& a, int sh) :arr(a), shift(sh) {}
    // return shifted element.
    std::pair<int, int>& operator[](int i) {
        return arr[((long long)i + arr.size() + shift) % arr.size()];
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

    ShiftedArr poly(verts, line.b == 0 ? max_y_ind : min);
    std::vector<std::pair<double,double>> answer;
    std::vector<int> answer_edges;
   
    int l = 0;
    int r = poly.arr.size();
    bool is_first_half = false;
    // check if line || Ox
       if (line.b == 0) {          
           long long value = (line.a * (long long)poly[l].first + line.c); // (a*x_1+c)
           auto sign = [](long long x) {return x > 0 ? 1 : x < 0 ? -1 : 0; };
           if (sign(value) * sign(line.a) != 1) {             // right part
               r = min_y_ind;
               if (r == 0)
                   r = poly.arr.size();
               is_first_half = true;
           }
           else {//left part
               l = min_y_ind;
           }
       }
       else {
           // std::cout << "Init Left point: " << poly[l].first << " " << poly[l].second <<"index: "<<l<< "\n";
            //std::cout << "Init Right point: " << poly[r].first << " " << poly[r].second << "index: " << r << "\n";
                //deside which part of polygon we should choose

           long long value_l = -(line.a * (long long)poly[min].first + line.c); // -(a*x_0+c)
           long long value_r = line.b * (long long)poly[min].second; //   y*b

           if (((value_l < value_r) && line.b > 0) || ((value_l > value_r) && line.b < 0)) { // bottom part
               l = max;
               // std::cout << poly[1].second << "\n";
              //  std::cout <<"Bottom" << "\n";
           }
           else {  // top part
               r = max;
               if (r == 0)
                   r = poly.arr.size();
               is_first_half = true;
               //std::cout << "Top" << "\n";
           }
       }
        //  if intersection point exists
        auto find_answer = [&](int left, int right) {
            if (line.sign(poly[left]) * line.sign(poly[right]) == 1) {

                return;
            }
            int v = poly.binSearch(left, right, line);
            Edge edge(poly[v], poly[v + 1]);
            auto point = edge.findIntersection(line, answer);
            //check if point is vertex
            std::pair<int, int> int_point = std::make_pair(point.first + 0.5 - (point.first < 0),
                point.second + 0.5 - (point.second < 0));
            if (int_point == poly[v]) {
                answer_edges.push_back(v);
            }
            if (int_point == poly[v+1]) {
                answer_edges.push_back(v+1);
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
            int minimum_of_the_fun = poly.ternaryMinSearch(l, r, line, line.sign(poly[l]));
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
            if (sign_l * sign_mid != 1) {
                find_answer(l, mid);
            }
            if (sign_r * sign_mid != 1) {
                find_answer(mid, r);
            }

            if (answer.size() == 1 && sign_l * sign_r != 1) { // can be another point
                int ind = (line.b == 0 ? min_y_ind : max);
                if (is_first_half) {
                    //check bottom 
                    find_answer(ind, poly.arr.size() - 1);
                }
                else {
                    find_answer(0, ind);
                }
            }
        }
        //check parallel   
        if (answer_edges.size() == 2 && (abs((answer_edges[0] - answer_edges[1])) == 1 ||
            abs((answer_edges[0] - answer_edges[1])) == poly.arr.size()-1)) {
            answer.push_back(poly[answer_edges[0]]);
            answer.push_back(poly[answer_edges[1]]);
            answer.push_back(std::make_pair(INT_MAX, INT_MAX));
        }
    return answer;
}
//for testing
//void run(Line& line, std::vector<pt> points, std::string out_name) {
//    std::ofstream out(out_name);
//    std::pair<std::pair<int, int>, int> min_x;
//    std::pair<std::pair<int, int>, int> max_x;
//    std::vector<std::pair<int, int>> poly;
//    for (int i = 0; i != points.size(); ++i) {
//        int x=points[i].x, y=points[i].y;
//        poly.push_back({ x,y });
//        if (i) {
//            min_x = min(min_x, { {x,y},i });
//            max_x = max(max_x, { {x,y},i });
//        }
//        else {
//            min_x = { {x,y},i };
//            max_x = { {x,y},i };
//        }
//    }
//    std::vector<std::pair<int, int>> res = solve(poly, line, min_x.second, max_x.second);             // this is answer
//    res.erase(std::unique(res.begin(), res.end()), res.end());
//    if (res.size() > 2) {
//        // Infinity
//        out << -1 << "\n";
//        out << res[0].first << ' ' << res[0].second << ' ';
//        out << res[1].first << ' ' << res[1].second;
//
//        out << '\n';
//    }
//    else {
//        out << res.size() << '\n';
//        for (auto coord : res) {
//            out << coord.first << ' ' << coord.second << ' ';
//        }
//        out << '\n';
//    }
//    out.close();
//}
//void main() {
//    std::cout<<test(239);
//}
int main() {
    std::string line;  
    std::ifstream in("input.txt"); 
    std::ofstream out("out.txt");
    if (in.is_open())
    {
        int n = -1;
        in >> n;
        std::pair<std::pair<int,int>, int> min_x;
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
                min_y = min(min_y, { {x,y},i });
                max_y = max(max_y, { {x,y},i });
            }
            else {
                min_x =  { {x,y},i };
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
            lines.emplace_back(Line(a,b,c));
        }
        if (out.is_open()) {
            for (auto line : lines) {
                auto res= solve(verts, line, min_x.second, max_x.second, min_y.second, max_y.second);             // this is answer
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
