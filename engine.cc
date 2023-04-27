#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <list>
#include <cmath>
#include <algorithm>
#include "vector3d.h"
#include "vector3d.cc"
#include <limits>
#include "ZBuffer.h"
#include "Point2D.h"
#include "Figure.h"
#include "Face.h"

/*Classes, namespaces and typedefs*/
const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
bool test_is_on = false;

img::Color vectorToColor(std::vector<double> kleur){
    img::Color to_return = img::Color(kleur[0]*255, kleur[1]*255, kleur[2]*255);
    return to_return;
}
double to_radialen(double graden){
    return graden*pi/180;
}
inline int roundToInt(double d)
{
    return static_cast<int>(std::round(d));
}
class Line2D{
public:
    Point2D a;
    Point2D b;
    img::Color color;
    double z1;
    double z2;
    Line2D(){};
    Line2D(Point2D ap, Point2D bp, img::Color kleur){
        a=ap;
        b=bp;
        color = kleur;
    }
    Line2D(Point2D ap, Point2D bp, img::Color kleur, double z_1, double z_2){
        a=ap;
        b=bp;
        color = kleur;
        z1 = z_1;
        z2 = z_2;
    }

};
using Lines2D = std::vector<Line2D>;

typedef std::list<Figure> Figures3D;

/*Functions*/

std::pair<double,double> getMinimum(const Lines2D &lines){
    double x = lines[0].a.x;
    double y = lines[0].a.y;
    for(auto lijn: lines){
        if (std::isinf(lijn.a.x) || std::isinf(-lijn.a.x)){
            lijn.a.x = 0;
            std::cout << "Inf detected \n";
        }
        if (std::isinf(lijn.a.y) || std::isinf(-lijn.a.y)){
            lijn.a.y = 0;
            std::cout << "Inf detected \n";
        }
        if (lijn.a.x < x) x = lijn.a.x;
        if (lijn.b.x < x) x = lijn.b.x;
        if (lijn.a.y < y) y = lijn.a.y;
        if (lijn.b.y < y) y = lijn.b.y;
    }
    std::pair<double,double> to_return = std::make_pair(x,y);
    return to_return;
}
std::pair<double,double> getMaximum(const Lines2D &lines){
    double x = lines[0].a.x;
    double y = lines[0].a.y;
    for(auto lijn: lines){
        if (lijn.a.x > x) x = lijn.a.x;
        if (lijn.b.x > x) x = lijn.b.x;
        if (lijn.a.y > y) y = lijn.a.y;
        if (lijn.b.y > y) y = lijn.b.y;
    }
    std::pair<double,double> to_return = std::make_pair(x,y);
    return to_return;
}

std::vector<Face> triangulate(const Face& face){
    std::vector<Face> to_return;
    int n = face.point_indexes.size();
    for(int i = 1; i < n - 1; i++){
        to_return.push_back(Face({face.point_indexes[0], face.point_indexes[i], face.point_indexes[i + 1]}));
    }
    return to_return;
}
//template <typename T>
//unsigned long findINdex(const std::vector<T> &vec, const T &toFind) {
//    unsigned long i;
//    for (auto e : vec) {
//        if (e == toFind) { return i; }
//        i++;
//    }
//    std::cerr << "Element not found" << std::endl;
//    return 0;
//}

img::EasyImage draw3DLines(const Lines2D &lines, const int size, img::Color background_color, const ini::Configuration &configuration){
    double x_min = getMinimum(lines).first;
    double y_min = getMinimum(lines).second;
    double x_max = getMaximum(lines).first;
    double y_max = getMaximum(lines).second;
    // Bereken x_range en y_range
    double x_range = x_max - x_min;
    double y_range = y_max - y_min;
    // Bereken imagex
    int imagex = static_cast<int>(std::round(size*x_range/std::max(x_range, y_range)));
    // Bereken imagey
    int imagey = static_cast<int>(std::round(size*y_range/std::max(x_range, y_range)));
    // Bereken d
    double d = 0.95*imagex/x_range;
    // Bereken DCx
    double dcx = d*(x_min+x_max)/2;
    // Bereken DCy
    double dcy = d*(y_min+y_max)/2;
    //dx
    double dx = imagex/2 - dcx;
    //dy
    double dy = imagey/2 - dcy;
    if(imagey < 1) imagey = 1;
    if(imagex < 1) imagex = 1;
    img::EasyImage to_return(imagex, imagey, background_color);
    for(auto lijn: lines){
        // Vermenigvuldig alle punten met d
        lijn.a.x = d*lijn.a.x + dx;
        lijn.a.y = d*lijn.a.y + dy;
        lijn.b.x = d*lijn.b.x + dx;
        lijn.b.y = d*lijn.b.y + dy;
        if(lijn.a.x != lijn.b.x || lijn.a.y != lijn.b.y){
            to_return.draw_zbuf_line(roundToInt(lijn.a.x), roundToInt(lijn.a.y), roundToInt(lijn.a.z), roundToInt(lijn.b.x),
                                     roundToInt(lijn.b.y), lijn.b.z, lijn.color);
        }
    }
    return to_return;
}
img::EasyImage draw2DLines(const Lines2D &lines, const int size, img::Color background_color, const ini::Configuration &configuration){
    // Bereken x_min, x_max, y_min, y_max;
    double x_min = getMinimum(lines).first;
    double y_min = getMinimum(lines).second;
    double x_max = getMaximum(lines).first;
    double y_max = getMaximum(lines).second;
    // Bereken x_range en y_range
    double x_range = x_max - x_min;
    double y_range = y_max - y_min;
    // Bereken imagex
    int imagex = roundToInt(std::round(size*x_range/std::max(x_range, y_range)));
    // Bereken imagey
    int imagey = roundToInt(std::round(size*y_range/std::max(x_range, y_range)));
    // Bereken d
    double d = 0.95*imagex/x_range;
    // Bereken DCx
    double dcx = d*(x_min+x_max)/2.0;
    // Bereken DCy
    double dcy = d*(y_min+y_max)/2.0;
    //dx
    double dx = imagex/2.0 - dcx;
    //dy
    double dy = imagey/2.0 - dcy;
    if(imagey < 1) imagey = 1;
    if(imagex < 1) imagex = 1;
    img::EasyImage to_return(imagex, imagey, background_color);
    for(auto lijn: lines){
        // Vermenigvuldig alle punten met d
        lijn.a.x = d*lijn.a.x + dx;
        lijn.a.y = d*lijn.a.y + dy;
        lijn.b.x = d*lijn.b.x + dx;
        lijn.b.y = d*lijn.b.y + dy;
        if(lijn.a.x != lijn.b.x || lijn.a.y != lijn.b.y){
            to_return.draw_line(roundToInt(lijn.a.x), roundToInt(lijn.a.y), roundToInt(lijn.b.x), roundToInt(lijn.b.y), lijn.color);
        }
    }
    return to_return;
}
Lines2D drawLSystem(const LParser::LSystem2D &l_system, const ini::Configuration &configuration){
    std::vector<std::pair<double,double>> positie_stack;
    std::vector<double> hoek_stack;
    std::set<char> alfabet = l_system.get_alphabet();
    // Gehardcoded kleur
    img::Color kleur = vectorToColor(configuration["2DLSystem"]["color"]);
    // Gehardcoded grootte van de stap
    double stap_grootte = 100.0;
    Lines2D to_return;
    // Get oorspronkelijke string
    std::string begin = l_system.get_initiator();
    // String to_return
    std::string res = "";
    // Loop alle strings over (char per char) om tot finale string te bekomen
    int iterations = l_system.get_nr_iterations();
    while (iterations > 0){
        for(char sym: begin){
            if(alfabet.find(sym) != alfabet.end()) res += l_system.get_replacement(sym);
            else res += sym;
        }
        begin = res;
        res = "";
        iterations--;
    }
    // Set oorspronkelijke parameters
    double huidige_hoek = l_system.get_starting_angle()*pi/180;
    double hoek = l_system.get_angle()*pi/180;
    std::pair<double,double> positie = std::make_pair(0,0);
    for(char sym:begin){
        // Als symbool in alfabet zit:
        if(alfabet.find(sym) != alfabet.end()){
            // Voor elke stap voeg een nieuwe lijn toe aan de Lines 2D
            std::pair<double,double> co1 = positie;
            std::pair<double,double> co2 = std::make_pair(positie.first + stap_grootte*std::cos(huidige_hoek), positie.second + stap_grootte*std::sin(huidige_hoek));
            Point2D p1 = Point2D(co1);
            Point2D p2 = Point2D(co2);
            Line2D lijn = Line2D(p1,p2,kleur);
            // Voeg de lijn als we lijn als dat nodig is:
            if (l_system.draw(sym)) to_return.push_back(lijn);
            // Update positie
            positie = co2;
        }
        // Als niet:
        else{
            if(sym == '+') huidige_hoek += hoek;
            else if(sym == '-') huidige_hoek -= hoek;
            else if(sym == '(') {
                positie_stack.push_back(positie);
                hoek_stack.push_back(huidige_hoek);
            }
            else if(sym == ')') {
                positie = positie_stack[positie_stack.size()-1];
                huidige_hoek = hoek_stack[hoek_stack.size()-1];
                positie_stack.pop_back();
                hoek_stack.pop_back();
            }
            // In andere gevallen skip
        }
    }
    if(test_is_on) std::cout << begin << std::endl;
    return to_return;
}
Figure draw3DLSystem(const LParser::LSystem3D &l_system, const ini::Configuration &figConfig){
    std::vector<Vector3D> positie_stack;
    std::vector<Vector3D> H_stack;
    std::vector<Vector3D> L_stack;
    std::vector<Vector3D> U_stack;
    std::set<char> alfabet = l_system.get_alphabet();
    double hoek = to_radialen(l_system.get_angle());
    // Gehardcoded grootte van de stap
    double stap_grootte = 1;
    Figure to_return;
    // Get oorspronkelijke string
    std::string begin = l_system.get_initiator();
    /////////////////
    Vector3D H = Vector3D::vector(1, 0, 0);
    Vector3D L = Vector3D::vector(0, 1, 0);
    Vector3D U = Vector3D::vector(0, 0, 1);
    Vector3D positie = Vector3D::point(0,0,0);
    std::string res = "";
    int iterations = l_system.get_nr_iterations();
    while (iterations > 0){
        for(char sym: begin){
            if(alfabet.find(sym) != alfabet.end()) res += l_system.get_replacement(sym);
            else res += sym;
        }
        begin = res;
        res = "";
        iterations--;
    }

    for(char sym:begin){
        // Als symbool in alfabet zit:
        if(alfabet.find(sym) != alfabet.end()){
            // Voor elke stap voeg een nieuwe lijn toe aan de Lines 2D
            Vector3D co1 = positie;
            Vector3D co2 = Vector3D::point(positie.x + H.x, positie.y + H.y, positie.z + H.z);
            // Voeg de lijn als dat nodig is:
            if (l_system.draw(sym)) {
                to_return.points.push_back(co1);
                to_return.points.push_back(co2);
                to_return.faces.push_back(Face({static_cast<int>(to_return.points.size()-2), static_cast<int>(to_return.points.size()-1)}));
            }
            // Update positie
            positie = co2;
        }
            // Als niet:
        else{
            if(sym == '+' ||sym == '-') {
                double delta = hoek;
                if(sym == '-') delta *= -1;
                Vector3D Ht = H;
                H = H*std::cos(delta) + L*std::sin(delta);
                L = -Ht*std::sin(delta) + L*std::cos(delta);
                H.normalise();
                L.normalise();
            }
            else if(sym == '^' || sym == '&'){
                double delta = hoek;
                if(sym == '&') delta *= -1;
                Vector3D Ht = H;
                H = H*std::cos(delta) + U*std::sin(delta);
                U = -Ht*std::sin(delta) + U*std::cos(delta);
                H.normalise();
                U.normalise();
            }
            else if(sym == '\\' || sym == '/'){
                double delta = hoek;
                if(sym == '/') delta *= -1;
                Vector3D Lt = L;
                L = L*std::cos(delta) - U*std::sin(delta);
                U = Lt*std::sin(delta) + U*std::cos(delta);
                U.normalise();
                L.normalise();
            }
            else if(sym == '|'){
                H = -H;
                L = -L;
            }
            else if(sym == '(') {
                positie_stack.push_back(positie);
                H_stack.push_back(H);
                L_stack.push_back(L);
                U_stack.push_back(U);
            }
            else if(sym == ')') {
                positie = positie_stack[positie_stack.size()-1];
                H = H_stack[H_stack.size()-1];
                L = L_stack[L_stack.size()-1];
                U = U_stack[U_stack.size()-1];
                positie_stack.pop_back();
                H_stack.pop_back();
                L_stack.pop_back();
                U_stack.pop_back();
            }
            // In andere gevallen skip
        }
    }
    if(test_is_on) {
        std::cout << begin << std::endl;
    }
    return to_return;
}
Matrix scaleFigure(const double scale){
    Matrix to_return;
    for(auto i:{1,2,3}){
        to_return(i,i) = to_return(i,i)*scale;
    }
    return to_return;
}
Matrix rotateX(const double angle){
    Matrix to_return;
    double cangle = to_radialen(angle);
    to_return(1,1) = 1;
    to_return(2,2) = std::cos(cangle);
    to_return(3,3) = std::cos(cangle);
    to_return(2,3) = std::sin(cangle);
    to_return(3,2) = -std::sin(cangle);
    return to_return;
}
Matrix rotateY(const double cangle){
    Matrix to_return;
    double angle = to_radialen(cangle);
    to_return(1,1) = std::cos(angle);
    to_return(1,3) = -std::sin(angle);
    to_return(3,1) = std::sin(angle);
    to_return(3,3) = std::cos(angle);
    return to_return;
}
Matrix rotateZ(const double cangle){
    Matrix to_return;
    double angle = to_radialen(cangle);
    to_return(1,1) = std::cos(angle);
    to_return(2,1) = -std::sin(angle);
    to_return(1,2) = std::sin(angle);
    to_return(2,2) = std::cos(angle);
    return to_return;
}
Matrix translate(const Vector3D &vector){
    Matrix to_return;
    to_return(4,1) = vector.x;
    to_return(4,2) = vector.y;
    to_return(4,3) = vector.z;
    return to_return;
}
void toPolar(const Vector3D &point, double &theta, double &phi, double &r){
    theta = std::atan2(point.y, point.x);
    r = std::sqrt(std::pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    //if(r != 0) phi = std::acos(r);
    if(r != 0) phi = std::acos((point.z)/r);
    else phi = 0;
}
Matrix eyePointTrans(const Vector3D &eyepoint){
    double theta, phi, r;
    toPolar(eyepoint, theta, phi, r);
    Matrix to_return;
    to_return(1,1) = -std::sin(theta);
    to_return(1,2) = -std::cos(theta)*std::cos(phi);
    to_return(1,3) = std::sin(phi)*std::cos(theta);
    to_return(1,4) = 0;

    to_return(2,1) = std::cos(theta);
    to_return(2,2) = -std::sin(theta)*std::cos(phi);
    to_return(2,3) = std::sin(theta)*std::sin(phi);
    to_return(2,4) = 0;

    to_return(3,1) = 0;
    to_return(3,2) = std::sin(phi);
    to_return(3,3) = std::cos(phi);
    to_return(3,4) = 0;

    to_return(4,1) = 0;
    to_return(4,2) = 0;
    to_return(4,3) = -r;
    to_return(4,4) = 1;


    //Matrix to_return = rotateZ(pi/2 + theta)* rotateX(phi)* translate(Vector3D::point(0, 0, -r));
    return to_return;
}
void applyTransformation(Figures3D &figs, const Matrix &m){
    for(auto fig: figs){
        for(auto point: fig.points) point = point * m;
    }
}
void applyTransformation(Figure &fig, const Matrix &m){
        for(Vector3D & point: fig.points) {
            point = point*m;
        }
}
Point2D doProjection(const Vector3D &point, const double d){
    Point2D to_return;
    to_return.x = (d*point.x)/(-point.z);
    to_return.y = (d*point.y)/(-point.z);
    if(point.z == 0) {
        to_return.x = 0;
        to_return.y = 0;
    }
    to_return.z = point.z;
    return to_return;
}
Vector3D findMiddle(Vector3D a, Vector3D b){
    return Vector3D::point((a.x+b.x)/2, (a.y+b.y)/2, (a.z+b.z)/2);
}
Lines2D doProjection(const Figure & figuur){
    Lines2D to_return;
        std::vector<Point2D> currentPoints;
        for(auto point: figuur.points){
            currentPoints.push_back(doProjection(point, 1));
        }
        // Now we have to work with faces. We add a new line with color of the figure (line goes from prev to current point)
        for(auto face: figuur.faces){
            int currentIndex = face.point_indexes[0];
            int nextIndex = face.point_indexes[1];
            if(face.point_indexes.size() == 2) to_return.push_back(Line2D(currentPoints[currentIndex],currentPoints[nextIndex], figuur.color, currentPoints[currentIndex].z, currentPoints[nextIndex].z));
            else {
                int i = 0;
                while(i < face.point_indexes.size()-1){
                    to_return.push_back(Line2D(currentPoints[face.point_indexes[i]],currentPoints[face.point_indexes[i+1]], figuur.color, currentPoints[face.point_indexes[i]].z,currentPoints[face.point_indexes[i+1]].z));
                    i++;
                }
                to_return.push_back(Line2D(currentPoints[face.point_indexes[face.point_indexes.size()-1]],currentPoints[face.point_indexes[0]], figuur.color, currentPoints[face.point_indexes[face.point_indexes.size()-1]].z,currentPoints[face.point_indexes[0]].z));
            }
        }
    return to_return;
}

Figure createIcosahedron(){
    Figure figuur;
    figuur.points.push_back(Vector3D::point(0,0, std::sqrt(5)/2));
    for(int l = 2; l < 7; l++){
        figuur.points.push_back(Vector3D::point(std::cos(2*pi*(l-2)/5), std::sin(2*pi*(l-2)/5), 0.5));
    }
    for(int l = 7; l < 12; l++){
        figuur.points.push_back(Vector3D::point(std::cos((pi/5)+((l-7)*2*pi/5)), std::sin((pi/5)+((l-7)*2*pi/5)), -0.5));
    }
    figuur.points.push_back(Vector3D::point(0,0, -std::sqrt(5)/2));

    Face f1 = Face({0,1,2});
    Face f2;
    f2.point_indexes.push_back(0);
    f2.point_indexes.push_back(2);
    f2.point_indexes.push_back(3);
    Face f3;
    f3.point_indexes.push_back(0);
    f3.point_indexes.push_back(3);
    f3.point_indexes.push_back(4);
    Face f4;
    f4.point_indexes.push_back(0);
    f4.point_indexes.push_back(4);
    f4.point_indexes.push_back(5);
    Face f5;
    f5.point_indexes.push_back(0);
    f5.point_indexes.push_back(5);
    f5.point_indexes.push_back(1);
    Face f6;
    f6.point_indexes.push_back(1);
    f6.point_indexes.push_back(6);
    f6.point_indexes.push_back(2);
    Face f7;
    f7.point_indexes.push_back(2);
    f7.point_indexes.push_back(6);
    f7.point_indexes.push_back(7);
    Face f8;
    f8.point_indexes.push_back(2);
    f8.point_indexes.push_back(7);
    f8.point_indexes.push_back(3);
    Face f9;
    f9.point_indexes.push_back(3);
    f9.point_indexes.push_back(7);
    f9.point_indexes.push_back(8);
    Face f10;
    f10.point_indexes.push_back(3);
    f10.point_indexes.push_back(8);
    f10.point_indexes.push_back(4);
    Face f11;
    f11.point_indexes.push_back(4);
    f11.point_indexes.push_back(8);
    f11.point_indexes.push_back(9);
    Face f12;
    f12.point_indexes.push_back(4);
    f12.point_indexes.push_back(9);
    f12.point_indexes.push_back(5);
    Face f13;
    f13.point_indexes.push_back(5);
    f13.point_indexes.push_back(9);
    f13.point_indexes.push_back(10);
    Face f14;
    f14.point_indexes.push_back(5);
    f14.point_indexes.push_back(10);
    f14.point_indexes.push_back(1);
    Face f15;
    f15.point_indexes.push_back(1);
    f15.point_indexes.push_back(10);
    f15.point_indexes.push_back(6);
    Face f16;
    f16.point_indexes.push_back(11);
    f16.point_indexes.push_back(7);
    f16.point_indexes.push_back(6);
    Face f17;
    f17.point_indexes.push_back(11);
    f17.point_indexes.push_back(8);
    f17.point_indexes.push_back(7);
    Face f18;
    f18.point_indexes.push_back(11);
    f18.point_indexes.push_back(9);
    f18.point_indexes.push_back(8);
    Face f19;
    f19.point_indexes.push_back(11);
    f19.point_indexes.push_back(10);
    f19.point_indexes.push_back(9);
    Face f20;
    f20.point_indexes.push_back(11);
    f20.point_indexes.push_back(6);
    f20.point_indexes.push_back(10);

    figuur.faces.push_back(f1);
    figuur.faces.push_back(f2);
    figuur.faces.push_back(f3);
    figuur.faces.push_back(f4);
    figuur.faces.push_back(f5);
    figuur.faces.push_back(f6);
    figuur.faces.push_back(f7);
    figuur.faces.push_back(f8);
    figuur.faces.push_back(f9);
    figuur.faces.push_back(f10);
    figuur.faces.push_back(f11);
    figuur.faces.push_back(f12);
    figuur.faces.push_back(f13);
    figuur.faces.push_back(f14);
    figuur.faces.push_back(f15);
    figuur.faces.push_back(f16);
    figuur.faces.push_back(f17);
    figuur.faces.push_back(f18);
    figuur.faces.push_back(f19);
    figuur.faces.push_back(f20);
    return figuur;
}
void herschaalPuntenBal(Vector3D& punt){
    double r = std::sqrt(std::pow(punt.x, 2) + std::pow(punt.y, 2) + std::pow(punt.z, 2));
    punt.x /= r;
    punt.y /= r;
    punt.z /= r;
}
img::EasyImage generate_image(const ini::Configuration &configuration)
{
    int size = configuration["General"]["size"];
    std::vector<double> vecback = configuration["General"]["backgroundcolor"];
    img::Color backgroundcolor = vectorToColor(vecback);
    img::EasyImage to_return(size,size, backgroundcolor);
    std::string type = configuration["General"]["type"];
    if(type == "2DLSystem"){
        LParser::LSystem2D l_systeem;
        std::ifstream input_stream(configuration["2DLSystem"]["inputfile"]);
        input_stream >> l_systeem;
        input_stream.close();
        Lines2D lijst = drawLSystem(l_systeem, configuration);
        to_return = draw2DLines(lijst, configuration["General"]["size"], vectorToColor(configuration["General"]["backgroundcolor"]), configuration);
    }
    else if(type == "Wireframe" || type == "ZBufferedWireframe" || type == "ZBuffering"){
        int aantalf = configuration["General"]["nrFigures"];
        Lines2D toDraw;
        std::vector<Figure> alle_figuren;
        // Iterate over all figures
        for(int numb = 0; numb < aantalf; numb++) {

            auto figConfig = configuration["Figure" + std::to_string(numb)];
            std::string typefig = figConfig["type"];
            Matrix scale = scaleFigure(figConfig["scale"]);
            Matrix X = rotateX(figConfig["rotateX"]);
            Matrix Y = rotateY(figConfig["rotateY"]);
            Matrix Z = rotateZ(figConfig["rotateZ"]);
            std::vector<double> centerv = figConfig["center"];
            Vector3D center = Vector3D::point(centerv[0], centerv[1], centerv[2]);
            Matrix T = translate(center);
            //Get all transformation matrices
            std::vector<double> eyevec = configuration["General"]["eye"];
            Vector3D eye = Vector3D::point(eyevec[0], eyevec[1], eyevec[2]);
            Matrix eyeTransf = eyePointTrans(eye);
            Matrix finalTrans = scale * X * Y * Z * T * eyeTransf;
            img::Color kleur = vectorToColor(figConfig["color"]);

            Figure figuur;
            if(typefig.find("Fractal") != typefig.back()){
                typefig = typefig.substr(typefig.find("Fractal") + 7, typefig.size()-7);
            }
            if (typefig == "LineDrawing") {
                int aantalp = figConfig["nrPoints"];
                std::vector<Vector3D> points;
                // Add all points
                for (int p = 0; p < aantalp; p++) {
                    std::vector<double> vecpoint = figConfig["point" + std::to_string(p)];
                    Vector3D point = Vector3D::point(vecpoint[0], vecpoint[1], vecpoint[2]);
                    points.push_back(point);
                }
                // Add all faces
                int aantall = figConfig["nrLines"];
                std::vector<Face> faces;
                for (int l = 0; l < aantall; l++) {
                    Face face;
                    std::vector<int> fp = figConfig["line" + std::to_string(l)];
                    face.point_indexes = fp;
                    faces.push_back(face);
                }
                // Initialise figure
                figuur.points = points;
                figuur.color = kleur;
                figuur.faces = faces;

//                // Use the finalTrans matrix
//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if (typefig == "Cube") {

                figuur.cube();
                figuur.color = kleur;

//                // Use the finalTrans matrix
//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if (typefig == "Tetrahedron") {

                figuur.tetrahedron();
                figuur.color = kleur;

//                // Use the finalTrans matrix
//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());

            }
            else if (typefig == "Octahedron") {

                figuur.octahedron();
                figuur.color = kleur;

//                // Use the finalTrans matrix
//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());

            }
            else if (typefig == "Icosahedron") {
                figuur = createIcosahedron();

                figuur.color = kleur;

//                // Use the finalTrans matrix
//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());

            }
            else if (typefig == "Dodecahedron") {
                figuur.dodecahedron();
                figuur.color = kleur;

//                // Use the finalTrans matrix
//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if (typefig == "Cylinder") {

                figuur.color = kleur;
                int n = figConfig["n"];
                double height = figConfig["height"];

                for (int ind = 0; ind < n; ind++) {
                    figuur.points.push_back(Vector3D::point(std::cos(2 * pi * ind / n), std::sin(2 * pi * ind / n), 0));
                    figuur.faces.push_back(Face({ind, (ind + 1) % n, n + (ind + 1) % (n), n + ind}));
                }
                for (int ind = 0; ind < n; ind++) {
                    figuur.points.push_back(
                            Vector3D::point(std::cos(2 * pi * ind / n), std::sin(2 * pi * ind / n), height));
                }
//                std::vector<int> intsLastFace;
//                for(int ind = n - 1; ind >= 0; ind--) intsLastFace.push_back(ind);
//                figuur.faces.push_back(Face(intsLastFace));

                // Voeg bovenvlak toe
                Face bovenvlak;
                for (int ind = 0; ind < n; ind++) {
                    bovenvlak.point_indexes.push_back(ind + n);
                }
                // Voeg ondervlak toe
                Face ondervlak;
                for (int ind = n - 1; ind >= 0; ind--) {
                    ondervlak.point_indexes.push_back(ind);
                }
                figuur.faces.push_back(bovenvlak);
                figuur.faces.push_back(ondervlak);
//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if (typefig == "Cone") {

                figuur.color = kleur;
                int n = figConfig["n"];
                double height = figConfig["height"];

                for (int ind = 0; ind < n; ind++) {
                    figuur.points.push_back(Vector3D::point(std::cos(2 * pi * ind / n), std::sin(2 * pi * ind / n), 0));
                    figuur.faces.push_back(Face({ind, (ind + 1) % n, n}));
                }
                figuur.points.push_back(Vector3D::point(0, 0, height));
                std::vector<int> intsLastFace;
                for (int ind = n - 1; ind >= 0; ind--) intsLastFace.push_back(ind);
                figuur.faces.push_back(Face(intsLastFace));

//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if (typefig == "Sphere") {
                figuur = createIcosahedron();

                figuur.color = kleur;
                int n = figConfig["n"];
                while (n > 0) {
                    std::vector<Face> newFaces;
                    std::vector<Vector3D> newPoints;
                    // Collect newFaces and newPoints
                    for (auto face: figuur.faces) {
                        Vector3D A = figuur.points[face.point_indexes[0]];
                        Vector3D B = figuur.points[face.point_indexes[1]];
                        Vector3D C = figuur.points[face.point_indexes[2]];
                        Vector3D D = findMiddle(A, B);
                        // ??? Swap E and F
                        Vector3D E = findMiddle(B, C);
                        Vector3D F = findMiddle(C, A);
                        int indexA = newPoints.size();
                        int indexB = newPoints.size() + 1;
                        int indexC = newPoints.size() + 2;
                        int indexD = newPoints.size() + 3;
                        int indexE = newPoints.size() + 4;
                        int indexF = newPoints.size() + 5;
                        for (auto letter: {A, B, C, D, E, F}) newPoints.push_back(letter);
                        std::vector<std::vector<Vector3D>> driehoeken = {{A, D, F},
                                                                         {B, E, D},
                                                                         {C, F, E},
                                                                         {D, E, F}};
                        Face f1 = Face({indexA, indexD, indexF});
                        Face f2 = Face({indexB, indexE, indexD});
                        Face f3 = Face({indexC, indexF, indexE});
                        Face f4 = Face({indexD, indexE, indexF});
                        for (auto f: {f1, f2, f3, f4}) {
                            newFaces.push_back(f);
                        }
                    }
                    // Use the finalTrans matrix

                    figuur.faces = newFaces;
                    figuur.points = newPoints;
                    n--;
                }
                // Herschaal alle punten
                for (auto &p: figuur.points) herschaalPuntenBal(p);

//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if (typefig == "Torus") {

                figuur.color = kleur;
                int n = figConfig["n"];
                double r = figConfig["r"];
                double R = figConfig["R"];
                int m = figConfig["m"];

                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        // Is dat juiste volgorde???
                        double u = 2 * i * pi / n;
                        double v = 2 * j * pi / m;
                        Vector3D p = Vector3D::point((R + r * std::cos(v)) * std::cos(u),
                                                     (R + r * std::cos(v)) * std::sin(u), r * std::sin(v));
                        figuur.points.push_back(p);
                        // Ind of the point i,j = i*m + j,
                        Face f = Face({i * m + j, ((i + 1) % n) * m + j, ((i + 1) % n) * m + (j + 1) % m,
                                       i * m + (j + 1) % m});
                        figuur.faces.push_back(f);
                    }
                }
//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if (typefig == "3DLSystem") {
                // Read l_system
                // drawLSystem new function, that gives lines back
                LParser::LSystem3D l_systeem;
                std::ifstream input_stream(figConfig["inputfile"]);
                input_stream >> l_systeem;
                input_stream.close();
                figuur = draw3DLSystem(l_systeem, configuration);
                figuur.color = kleur;
//                applyTransformation(figuur, finalTrans);
//                // Do projection
//                Lines2D to_add = doProjection(figuur);
//                // Insert getted lines
//                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            if (type == "ZBuffering") {
                std::vector<Face> faces;
                // Triangulate all faces of figure
                for (const auto &f: figuur.faces) {
                    std::vector<Face> to_add = triangulate(f);
                    faces.insert(faces.end(), to_add.begin(), to_add.end());
                }
                figuur.faces = faces;
                //alle_figuren.push_back(figuur);
            }
            std::string typefig_full = figConfig["type"];
            // Fractalen
            if (typefig_full.find("Fractal") != typefig_full.back()) {
                int nrIteration = figConfig["nrIterations"];
                double fractalScale = figConfig["fractalScale"];
                std::vector<Figure> fractals;
                figuur.generateFractal(fractals, nrIteration, fractalScale);
                // Save all fractals
                for (const auto& fig: fractals) {
                    alle_figuren.push_back(fig);
                }
                if(nrIteration == 0){
                    alle_figuren.push_back(figuur);
                }
            }
            else{
                alle_figuren.push_back(figuur);
            }
            for (auto &fi: alle_figuren) {
                // Use the finalTrans matrix
                applyTransformation(fi, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(fi);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
        }
        if(type == "ZBufferedWireframe") to_return = draw3DLines(toDraw, size, vectorToColor(configuration["General"]["backgroundcolor"]), configuration);
        else if (type == "Wireframe") to_return = draw2DLines(toDraw, size, vectorToColor(configuration["General"]["backgroundcolor"]), configuration);
        else if (type == "ZBuffering") {
            // Bereken all shit waarden
            double x_min = getMinimum(toDraw).first;
            double y_min = getMinimum(toDraw).second;
            double x_max = getMaximum(toDraw).first;
            double y_max = getMaximum(toDraw).second;
            // Bereken x_range en y_range
            double x_range = x_max - x_min;
            double y_range = y_max - y_min;
            // Bereken imagex
            double imagex = size*x_range/std::max(x_range, y_range);
            // Bereken imagey
            double imagey = size*y_range/std::max(x_range, y_range);
            to_return = img::EasyImage(static_cast<int>(std::round(imagex)), static_cast<int>(std::round(imagey)), backgroundcolor);
            double d = 0.95*imagex/x_range;
            // Bereken DCx
            double dcx = d*(x_min+x_max)/2.0;
            // Bereken DCy
            double dcy = d*(y_min+y_max)/2.0;
            //dx
            double dx = imagex/2.0 - dcx;
            //dy
            double dy = imagey/2.0 - dcy;
            if(imagey < 1) imagey = 1;
            if(imagex < 1) imagex = 1;
            for(auto fig:alle_figuren){
                for(auto fac: fig.faces){
                    int A = fac.point_indexes[0];
                    int B = fac.point_indexes[1];
                    int C = fac.point_indexes[2];
                    to_return.draw_zbuf_triag(fig.points[A], fig.points[B], fig.points[C],
                                              d, dx, dy, fig.color);
                }
            }
        }
    }
	return to_return;
}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                if (fin.peek() == std::istream::traits_type::eof()) {
                                    std::cout << "Ini file appears empty. Does '" <<
                                    fileName << "' exist?" << std::endl;
                                    continue;
                                }
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
