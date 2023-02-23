#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <list>
#include <cmath>
const double pi = 3.1415;
class Color{
public:
    Color(){};
    Color(double red, double green, double blue){
        red = red;
        green = green;
        blue = blue;
    }
    double red;
    double green;
    double blue;
};
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
    Color color;
    Line2D(){};
    Line2D(Point2D ap, Point2D bp, Color kleur){
        a=ap;
        b=bp;
        color = kleur;
    }
};
using Lines2D = std::vector<Line2D>;
std::pair<double,double> getMinimum(const Lines2D &lines){
    double x = lines[0].a.x;
    double y = lines[0].a.y;
    for(auto lijn: lines){
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

img::EasyImage draw2DLines(const Lines2D &lines, const int size){
    // Bereken x_min, x_max, y_min, y_max;
    double x_min = getMinimum(lines).first;
    double y_min = getMinimum(lines).second;
    double x_max = getMaximum(lines).first;
    double y_max = getMaximum(lines).second;
    // Bereken x_range en y_range
    int x_range = std::round(x_max - x_min);
    int y_range = std::round(y_max - y_min);
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
    img::EasyImage to_return(imagex, imagey);
    for(auto lijn: lines){
        // Vermenigvuldig alle punten met d
        lijn.a.x = std::round(d*lijn.a.x + dx);
        lijn.a.y = std::round(d*lijn.a.y + dy);
        lijn.b.x = std::round(d*lijn.b.x + dx);
        lijn.b.y = std::round(d*lijn.b.y + dy);
        img::Color trans_color = img::Color(100, 100, 100);
        to_return.draw_line(lijn.a.x, lijn.a.y, lijn.b.x, lijn.b.y, trans_color);
    }
    return to_return;
}
Lines2D drawLSystem(const LParser::LSystem2D &l_system){
    std::set<char> alfabet = l_system.get_alphabet();
    // Gehardcoded kleur
    Color kleur = Color(150,150,150);
    // Gehardcoded grootte van de stap
    int stap_grootte = 30;
    Lines2D to_return;
    // Get oorspronkelijke string
    std::string begin = l_system.get_initiator();
    // String to_return
    std::string res = begin;
    // Loop alle strings over (char per char)
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
    std::pair<int,int> positie = std::make_pair(0,0);
    for(char sym:begin){
        // Als symbool in alfabet zit:
        if(alfabet.find(sym) != alfabet.end()){
            // Voor elke stap voeg een nieuwe lijn toe aan de Lines 2D
            std::pair<double,double> co1 = positie;
            std::pair<double,double> co2 = std::make_pair(positie.first + stap_grootte*std::cos(huidige_hoek), positie.second + std::sin(huidige_hoek));
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
            // In andere gevallen skip
        }
    }
    return to_return;
}
img::EasyImage generate_image(const ini::Configuration &configuration)
{
    // Deel 1
    Point2D a = Point2D(0,0);
    Point2D b = Point2D(80,40);
    Color kleur = Color(150,150,150);
    Line2D lijn1 = Line2D(a,b, kleur);
    Point2D a2 = Point2D(0,0);
    Point2D b2 = Point2D(40,80);
    Color kleur2 = Color(150,150,150);
    Line2D lijn2 = Line2D(a2,b2, kleur2);
    Lines2D lijst = Lines2D();
    lijst.push_back(lijn1);
    lijst.push_back(lijn2);
    // ???: Volgende lijnen in aparte functie steken (en alle andere dingen in de juiste klassen)
    img::EasyImage to_return = draw2DLines(lijst, 30);

    LParser::LSystem2D l_systeem;
    std::ifstream input_stream(configuration["2DLSystem"]["inputfile"]);
    input_stream >> l_systeem;
    input_stream.close();
    lijst = drawLSystem(l_systeem);
    to_return = draw2DLines(lijst, 200);
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
