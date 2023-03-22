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

/*Classes, namespaces and typedefs*/
const double pi = 3.141592653589793238462653589;
bool test_is_on = false;

img::Color vectorToColor(std::vector<double> kleur){
    img::Color to_return = img::Color(kleur[0]*255, kleur[1]*255, kleur[2]*255);
    return to_return;
}
double to_radialen(double graden){
    return graden*pi/180;
}
class Point2D{
public:
    Point2D(){};
    Point2D(double xy, double yy){
        x=xy;
        y=yy;
    }
    Point2D(std::pair<double,double> co){
        x = co.first;
        y = co.second;
    }
    Point2D(double co1, double co2, double co3){
        x = co1;
        y = co2;
        z = co3;
    }
    double x;
    double y;
    double z;
};
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
class Face
{
public:
    Face(){};
    Face(std::vector<int> inds){
        point_indexes = inds;
    }
    //De indexen refereren naar
    //punten in de ‘points’ vector
    //van de Figure-klasse
    std::vector<int> point_indexes;
};
class Figure
{
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    img::Color color;
};
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
    int imagex = std::round(size*x_range/std::max(x_range, y_range));
    // Bereken imagey
    int imagey = std::round(size*y_range/std::max(x_range, y_range));
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
            to_return.draw_line(lijn.a.x, lijn.a.y, lijn.b.x, lijn.b.y, lijn.color);
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
                H = H*std::cos(delta) + L*std::sin(delta);
                L = -H*std::sin(delta) + L*std::cos(delta);
                H.normalise();
                L.normalise();
            }
            else if(sym == '^' || sym == '&'){
                double delta = hoek;
                if(sym == '&') delta *= -1;
                H = H*std::cos(delta) + U*std::sin(delta);
                U = -H*std::sin(delta) + U*std::cos(delta);
                H.normalise();
                U.normalise();
            }
            else if(sym == '\\' || sym == '/'){
                double delta = hoek;
                if(sym == '/') delta *= -1;
                L = L*std::cos(delta) - U*std::sin(delta);
                U = L*std::sin(delta) + U*std::cos(delta);
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
    if(test_is_on) std::cout << begin << std::endl;
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
    img::EasyImage to_return(size,size);
    std::string type = configuration["General"]["type"];
    if(type == "2DLSystem"){
        LParser::LSystem2D l_systeem;
        std::ifstream input_stream(configuration["2DLSystem"]["inputfile"]);
        input_stream >> l_systeem;
        input_stream.close();
        Lines2D lijst = drawLSystem(l_systeem, configuration);
        to_return = draw2DLines(lijst, configuration["General"]["size"], vectorToColor(configuration["General"]["backgroundcolor"]), configuration);
    }
    else if(type == "Wireframe" || type == "ZBufferedWireframe"){
        int aantalf = configuration["General"]["nrFigures"];
        Lines2D toDraw;
        for(int numb = 0; numb < aantalf; numb++){

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

            if(typefig == "LineDrawing") {
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
                Figure figuur;
                figuur.points = points;
                figuur.color = kleur;
                figuur.faces = faces;

                // Use the finalTrans matrix
                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if(typefig == "Cube"){
                Figure figuur;

                figuur.points.push_back(Vector3D::point(1,-1,-1));
                figuur.points.push_back(Vector3D::point(-1,1,-1));
                figuur.points.push_back(Vector3D::point(1,1,1));
                figuur.points.push_back(Vector3D::point(-1,-1,1));
                figuur.points.push_back(Vector3D::point(1,1,-1));
                figuur.points.push_back(Vector3D::point(-1,-1,-1));
                figuur.points.push_back(Vector3D::point(1,-1,1));
                figuur.points.push_back(Vector3D::point(-1,1,1));

                Face f1;
                f1.point_indexes.push_back(0);
                f1.point_indexes.push_back(4);
                f1.point_indexes.push_back(2);
                f1.point_indexes.push_back(6);
                Face f2;
                f2.point_indexes.push_back(4);
                f2.point_indexes.push_back(1);
                f2.point_indexes.push_back(7);
                f2.point_indexes.push_back(2);
                Face f3;
                f3.point_indexes.push_back(1);
                f3.point_indexes.push_back(5);
                f3.point_indexes.push_back(3);
                f3.point_indexes.push_back(7);
                Face f4;
                f4.point_indexes.push_back(5);
                f4.point_indexes.push_back(0);
                f4.point_indexes.push_back(6);
                f4.point_indexes.push_back(3);
                Face f5;
                f5.point_indexes.push_back(6);
                f5.point_indexes.push_back(2);
                f5.point_indexes.push_back(7);
                f5.point_indexes.push_back(3);
                Face f6;
                f6.point_indexes.push_back(0);
                f6.point_indexes.push_back(5);
                f6.point_indexes.push_back(1);
                f6.point_indexes.push_back(4);

                figuur.faces.push_back(f1);
                figuur.faces.push_back(f2);
                figuur.faces.push_back(f3);
                figuur.faces.push_back(f4);
                figuur.faces.push_back(f5);
                figuur.faces.push_back(f6);

                figuur.color = kleur;

                // Use the finalTrans matrix
                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if(typefig == "Tetrahedron"){
                Figure figuur;

                figuur.points.push_back(Vector3D::point(1,-1,-1));
                figuur.points.push_back(Vector3D::point(-1,1,-1));
                figuur.points.push_back(Vector3D::point(1,1,1));
                figuur.points.push_back(Vector3D::point(-1,-1,1));

                Face f1;
                f1.point_indexes.push_back(0);
                f1.point_indexes.push_back(1);
                f1.point_indexes.push_back(2);
                Face f2;
                f2.point_indexes.push_back(1);
                f2.point_indexes.push_back(3);
                f2.point_indexes.push_back(2);
                Face f3;
                f3.point_indexes.push_back(0);
                f3.point_indexes.push_back(3);
                f3.point_indexes.push_back(1);
                Face f4;
                f4.point_indexes.push_back(0);
                f4.point_indexes.push_back(2);
                f4.point_indexes.push_back(3);

                figuur.faces.push_back(f1);
                figuur.faces.push_back(f2);
                figuur.faces.push_back(f3);
                figuur.faces.push_back(f4);

                figuur.color = kleur;

                // Use the finalTrans matrix
                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());

            }
            else if(typefig == "Octahedron"){
                Figure figuur;

                figuur.points.push_back(Vector3D::point(1,0,0));
                figuur.points.push_back(Vector3D::point(0,1,0));
                figuur.points.push_back(Vector3D::point(-1,0,0));
                figuur.points.push_back(Vector3D::point(0,-1,0));
                figuur.points.push_back(Vector3D::point(0,0,-1));
                figuur.points.push_back(Vector3D::point(0,0,1));

                Face f1;
                f1.point_indexes.push_back(0);
                f1.point_indexes.push_back(1);
                f1.point_indexes.push_back(5);
                Face f2;
                f2.point_indexes.push_back(1);
                f2.point_indexes.push_back(2);
                f2.point_indexes.push_back(5);
                Face f3;
                f3.point_indexes.push_back(2);
                f3.point_indexes.push_back(3);
                f3.point_indexes.push_back(5);
                Face f4;
                f4.point_indexes.push_back(3);
                f4.point_indexes.push_back(0);
                f4.point_indexes.push_back(5);
                Face f5;
                f5.point_indexes.push_back(1);
                f5.point_indexes.push_back(0);
                f5.point_indexes.push_back(4);
                Face f6;
                f6.point_indexes.push_back(2);
                f6.point_indexes.push_back(1);
                f6.point_indexes.push_back(4);
                Face f7;
                f7.point_indexes.push_back(3);
                f7.point_indexes.push_back(2);
                f7.point_indexes.push_back(4);
                Face f8;
                f8.point_indexes.push_back(0);
                f8.point_indexes.push_back(3);
                f8.point_indexes.push_back(4);

                figuur.faces.push_back(f1);
                figuur.faces.push_back(f2);
                figuur.faces.push_back(f3);
                figuur.faces.push_back(f4);
                figuur.faces.push_back(f5);
                figuur.faces.push_back(f6);
                figuur.faces.push_back(f7);
                figuur.faces.push_back(f8);

                figuur.color = kleur;

                // Use the finalTrans matrix
                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());

            }
            else if(typefig == "Icosahedron"){
                Figure figuur = createIcosahedron();

                figuur.color = kleur;

                // Use the finalTrans matrix
                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());

            }
            else if(typefig == "Dodecahedron"){
                std::vector<Vector3D> icoPoints;
                Figure figuur;
                icoPoints.push_back(Vector3D::point(0,0, std::sqrt(5)/2));
                for(int l = 2; l < 7; l++){
                    icoPoints.push_back(Vector3D::point(std::cos(2*pi*(l-2)/5), std::sin(2*pi*(l-2)/5), 0.5));
                }
                for(int l = 7; l < 12; l++){
                    icoPoints.push_back(Vector3D::point(std::cos((pi/5)+((l-7)*2*pi/5)), std::sin((pi/5)+((l-7)*2*pi/5)), -0.5));
                }
                icoPoints.push_back(Vector3D::point(0,0, -std::sqrt(5)/2));
                // Here i'd like to give a point from old faces
                std::vector<std::vector<int>> icoInd = {{0,1,2}, {0,2,3}, {0,3,4},
                                                        {0,4,5}, {0,5,1}, {1,6,2},
                                                        {2,6,7}, {2,7,3}, {3,7,8},
                                                        {3,8,4}, {4,8,9}, {4,9,5},
                                                        {5,9,10}, {5,10,1}, {1,10,6},
                                                        {11,7,6}, {11,8,7}, {11,9,8},
                                                        {11,10,9}, {11,6,10}};
                for(auto point: icoInd){
                    figuur.points.push_back(Vector3D::point((icoPoints[point[0]].x + icoPoints[point[1]].x + icoPoints[point[2]].x)/3,
                                    (icoPoints[point[0]].y + icoPoints[point[1]].y + icoPoints[point[2]].y)/3,
                                    (icoPoints[point[0]].z + icoPoints[point[1]].z + icoPoints[point[2]].z)/3));
                }
                Face f1 = Face({0,1,2,3,4});
                Face f2 = Face({0,5,6,7,1});
                Face f3 = Face({1,7,8,9,2});
                Face f4 = Face({2,9,10,11,3});
                Face f5 = Face({3,11,12,13,4});
                Face f6 = Face({4,13,14,5,0});
                Face f7 = Face({19,18,17,16,15});
                Face f8;
                f8.point_indexes.push_back(19);
                f8.point_indexes.push_back(14);
                f8.point_indexes.push_back(13);
                f8.point_indexes.push_back(12);
                f8.point_indexes.push_back(18);
                Face f9;
                f9.point_indexes.push_back(18);
                f9.point_indexes.push_back(12);
                f9.point_indexes.push_back(11);
                f9.point_indexes.push_back(10);
                f9.point_indexes.push_back(17);
                Face f10;
                f10.point_indexes.push_back(17);
                f10.point_indexes.push_back(10);
                f10.point_indexes.push_back(9);
                f10.point_indexes.push_back(8);
                f10.point_indexes.push_back(16);
                Face f11;
                f11.point_indexes.push_back(16);
                f11.point_indexes.push_back(8);
                f11.point_indexes.push_back(7);
                f11.point_indexes.push_back(6);
                f11.point_indexes.push_back(15);
                Face f12;
                f12.point_indexes.push_back(15);
                f12.point_indexes.push_back(6);
                f12.point_indexes.push_back(5);
                f12.point_indexes.push_back(14);
                f12.point_indexes.push_back(19);

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

                figuur.color = kleur;

                // Use the finalTrans matrix
                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if(typefig == "Cylinder"){
                Figure figuur;

                figuur.color = kleur;
                int n = figConfig["n"];
                double height = figConfig["height"];

                for(int ind = 0; ind < n; ind++){
                    figuur.points.push_back(Vector3D::point(std::cos(2*pi*ind/n), std::sin(2*pi*ind/n), 0));
                    figuur.faces.push_back(Face({ind, (ind+1)%n, n + (ind+1)%(n),  n + ind}));
                }
                for(int ind = 0; ind < n; ind++){
                    figuur.points.push_back(Vector3D::point(std::cos(2*pi*ind/n), std::sin(2*pi*ind/n), height));
                }
//                std::vector<int> intsLastFace;
//                for(int ind = n - 1; ind >= 0; ind--) intsLastFace.push_back(ind);
//                figuur.faces.push_back(Face(intsLastFace));

                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if(typefig == "Cone"){
                Figure figuur;

                figuur.color = kleur;
                int n = figConfig["n"];
                double height = figConfig["height"];

                for(int ind = 0; ind < n; ind++){
                    figuur.points.push_back(Vector3D::point(std::cos(2*pi*ind/n), std::sin(2*pi*ind/n), 0));
                    figuur.faces.push_back(Face({ind, (ind+1)%n, n}));
                }
                figuur.points.push_back(Vector3D::point(0,0,height));
                std::vector<int> intsLastFace;
                for(int ind = n - 1; ind >= 0; ind--) intsLastFace.push_back(ind);
                figuur.faces.push_back(Face(intsLastFace));

                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if(typefig == "Sphere") {
                Figure figuur = createIcosahedron();

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
                for(auto &p:figuur.points) herschaalPuntenBal(p);

                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if(typefig == "Torus"){
                Figure figuur;

                figuur.color = kleur;
                int n = figConfig["n"];
                double r = figConfig["r"];
                double R = figConfig["R"];
                int m = figConfig["m"];

                for(int i = 0; i < n; i++){
                    for(int j = 0; j < m; j++){
                        // Is dat juiste volgorde???
                        double u = 2*i*pi/n;
                        double v = 2*j*pi/m;
                        Vector3D p = Vector3D::point((R+r*std::cos(v))*std::cos(u),(R+r*std::cos(v))*std::sin(u),r*std::sin(v));
                        figuur.points.push_back(p);
                        // Ind of the point i,j = i*m + j,
                        Face f = Face({i*m + j, ((i+1)%n)*m + (j+1)%m, i*m + (j+1)%m});
                        figuur.faces.push_back(f);
                    }
                }
                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            else if(typefig == "3DLSystem"){
                // Read l_system
                // drawLSystem new function, that gives lines back
                LParser::LSystem3D l_systeem;
                std::ifstream input_stream(figConfig["inputfile"]);
                input_stream >> l_systeem;
                input_stream.close();
                Figure figuur = draw3DLSystem(l_systeem, configuration);
                figuur.color = kleur;
                applyTransformation(figuur, finalTrans);
                // Do projection
                Lines2D to_add = doProjection(figuur);
                // Insert getted lines
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());
            }
            }
        if(type == "ZBufferedWireframe") to_return = draw3DLines(toDraw, size, vectorToColor(configuration["General"]["backgroundcolor"]), configuration);
        else to_return = draw2DLines(toDraw, size, vectorToColor(configuration["General"]["backgroundcolor"]), configuration);

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
