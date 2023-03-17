#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <list>
#include <cmath>
#include "vector3d.h"
#include "vector3d.cc"

/*Classes, namespaces and typedefs*/
const double pi = 3.141592653589793238462653589;
bool test_is_on = false;
img::Color vectorToColor(std::vector<double> kleur){
    img::Color to_return = img::Color(kleur[0]*255, kleur[1]*255, kleur[2]*255);
    return to_return;
}
//class Color{
//public:
//    Color(){};
//    Color(double red, double green, double blue){
//        red = red;
//        green = green;
//        blue = blue;
//    }
//    double red;
//    double green;
//    double blue;
//};
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
    double x;
    double y;
};
class Line2D{
public:
    Point2D a;
    Point2D b;
    img::Color color;
    Line2D(){};
    Line2D(Point2D ap, Point2D bp, img::Color kleur){
        a=ap;
        b=bp;
        color = kleur;
    }
};
using Lines2D = std::vector<Line2D>;
class Face
{
public:
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
    return to_return;
}
Lines2D doProjection(const Figures3D & figuren){
    Lines2D to_return;
    for(auto figuur: figuren){
        std::vector<Point2D> currentPoints;
        for(auto point: figuur.points){
            currentPoints.push_back(doProjection(point, 1));
        }
        // Now we have to work with faces. We add a new line with color of the figure (line goes from prev to current point)
        for(auto face: figuur.faces){
            int currentIndex = face.point_indexes[0];
            int nextIndex = face.point_indexes[1];
            while(nextIndex < face.point_indexes.size()){
                to_return.push_back(Line2D(currentPoints[currentIndex],currentPoints[nextIndex], figuur.color));
                currentIndex++;
                nextIndex++;
            }
        }
    }
    return  to_return;
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
            if(face.point_indexes.size() == 2) to_return.push_back(Line2D(currentPoints[currentIndex],currentPoints[nextIndex], figuur.color));
            else if (face.point_indexes.size() == 3) {
                // ??? Change sequence
                to_return.push_back(Line2D(currentPoints[0],currentPoints[1], figuur.color));
                to_return.push_back(Line2D(currentPoints[1],currentPoints[2], figuur.color));
                to_return.push_back(Line2D(currentPoints[2],currentPoints[0], figuur.color));
            }
            else{
                std::cout << "Fout aantal punten in de face \n";
            }
        }
    return  to_return;
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
    }else if(type == "Wireframe"){
        int aantalf = configuration["General"]["nrFigures"];
        Lines2D toDraw;
        for(int i = 0; i < aantalf; i++){
            auto figConfig = configuration["Figure" + std::to_string(i)];
            std::string typefig = figConfig["type"];
            if(typefig == "LineDrawing") {
                //Get all transformation matrices
                Matrix scale = scaleFigure(figConfig["scale"]);
                Matrix X = rotateX(figConfig["rotateX"]);
                Matrix Y = rotateY(figConfig["rotateY"]);
                Matrix Z = rotateZ(figConfig["rotateZ"]);
                std::vector<double> centerv = figConfig["center"];
                Vector3D center = Vector3D::point(centerv[0], centerv[1], centerv[2]);
                Matrix T = translate(center);
                std::vector<double> eyevec = configuration["General"]["eye"];
                Vector3D eye = Vector3D::point(eyevec[0], eyevec[1], eyevec[2]);
                Matrix eyeTransf = eyePointTrans(eye);
                Matrix finalTrans = scale * X * Y * Z * T * eyeTransf;

                img::Color kleur = vectorToColor(figConfig["color"]);
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
                toDraw.insert(toDraw.end(), to_add.begin(), to_add.end());}
            else if(typefig == "Cube"){
                /*
                 * Kubus (Cube)
                • Tetrahedron
                • Octahedron
                • Icosahedron
                • Dodecahedron
                • Cilinder (Cylinder)
                • Kegel (Cone)
                • Bol (Sphere)
                • Torus*/
            }
            else if(typefig == "Tetrahedron"){

            }
            else if(typefig == "Octahedron"){

            }
            else if(typefig == "Icosahedron"){

            }
            else if(typefig == "Dodecahedron"){

            }
            else if(typefig == "Cylinder"){

            }
            else if(typefig == "Cone"){

            }
            else if(typefig == "Sphere"){

            }
            else if(typefig == "Torus"){

            }
            else if(typefig == "3DLSystem"){

            }
            }
            to_return = draw2DLines(toDraw, size, vectorToColor(configuration["General"]["backgroundcolor"]), configuration);
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
