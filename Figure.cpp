//
// Created by dasha on 25.04.2023.
//

#include "Figure.h"
#include "Face.h"
#include "vector3d.h"
#include <cmath>
#include <fstream>
#include "obj_parser.h"
#include "Texture.h"
#include "utils.h"
#define _USE_MATH_DEFINES

Figure::Figure() {}
//Figure::Figure(Figure& fig){
//    points = std::vector<Vector3D> (fig.points);
//    faces = std::vector<Face> (fig.faces);
//    color = fig.color;
//}
void Figure::scaleFigure(const double scale){
    Matrix m;
    for(auto i:{1,2,3}){
        m(i,i) = m(i,i)*scale;
    }
    // Apply scale
    for(Vector3D & point: points) {
        point = point*m;
    }
}
void Figure::translate(const Vector3D &vector){
    Matrix to_return;
    to_return(4,1) = vector.x;
    to_return(4,2) = vector.y;
    to_return(4,3) = vector.z;
    // Apply scale
    for(Vector3D & point: points) {
        point = point*to_return;
    }
}
Figure Figure::copyFigure(const Figure& fig){
    Figure fract_fig;
    // Copy constructor
    fract_fig.color = fig.color;
    fract_fig.faces = std::vector<Face>(fig.faces);
    fract_fig.points = std::vector<Vector3D>(fig.points);
    fract_fig.diffuseReflection = std::vector<double>(fig.diffuseReflection);
    fract_fig.fullAmbientReflection = std::vector<double>(fig.fullAmbientReflection);
    fract_fig.specularReflection = std::vector<double>(fig.specularReflection);
    fract_fig.reflectionCoefficient = fig.reflectionCoefficient;
    fract_fig.textureNrs = fig.textureNrs;
    fract_fig.cubeMapNr = fig.cubeMapNr;
    return fract_fig;
}
void Figure::generateFractal(std::vector<Figure> & fractal, const int nr_iterations, const double scale) const{
    std::vector<Figure> figures_to_fractalise;

    // Copy constructor
    Figure figure = copyFigure(*this);

    figures_to_fractalise.push_back(figure);
    for(int n = 0; n < nr_iterations; n++) {
        for(const auto& fig:figures_to_fractalise) {
            for (int i = 0; i < fig.points.size(); i++) {
                Figure fract_fig = copyFigure(fig);
                // Scale figure
                fract_fig.scaleFigure(1 / scale);
                // Calculate vector
                Vector3D to_move = fig.points[i] - fract_fig.points[i];
                // Move
                fract_fig.translate(to_move);
                // Add figure to fractal
                fractal.push_back(fract_fig);
                // Do it for more iterations
            }
        }
        figures_to_fractalise = fractal;
        fractal.clear();
    }
    fractal = figures_to_fractalise;
}
void Figure::cube(){
    points.push_back(Vector3D::point(1,-1,-1));
    points.push_back(Vector3D::point(-1,1,-1));
    points.push_back(Vector3D::point(1,1,1));
    points.push_back(Vector3D::point(-1,-1,1));
    points.push_back(Vector3D::point(1,1,-1));
    points.push_back(Vector3D::point(-1,-1,-1));
    points.push_back(Vector3D::point(1,-1,1));
    points.push_back(Vector3D::point(-1,1,1));

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

    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
}
void Figure::tetrahedron(){
    points.push_back(Vector3D::point(1,-1,-1));
    points.push_back(Vector3D::point(-1,1,-1));
    points.push_back(Vector3D::point(1,1,1));
    points.push_back(Vector3D::point(-1,-1,1));

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

    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
}
void Figure::octahedron(){
    points.push_back(Vector3D::point(1,0,0));
    points.push_back(Vector3D::point(0,1,0));
    points.push_back(Vector3D::point(-1,0,0));
    points.push_back(Vector3D::point(0,-1,0));
    points.push_back(Vector3D::point(0,0,-1));
    points.push_back(Vector3D::point(0,0,1));

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

    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
    faces.push_back(f7);
    faces.push_back(f8);
}
void Figure::dodecahedron() {
    std::vector<Vector3D> icoPoints;
    icoPoints.push_back(Vector3D::point(0,0, std::sqrt(5)/2));
    for(int l = 2; l < 7; l++){
        icoPoints.push_back(Vector3D::point(std::cos(2*M_PI*(l-2)/5), std::sin(2*M_PI*(l-2)/5), 0.5));
    }
    for(int l = 7; l < 12; l++){
        icoPoints.push_back(Vector3D::point(std::cos((M_PI/5)+((l-7)*2*M_PI/5)), std::sin((M_PI/5)+((l-7)*2*M_PI/5)), -0.5));
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
        points.push_back(Vector3D::point((icoPoints[point[0]].x + icoPoints[point[1]].x + icoPoints[point[2]].x)/3,
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

    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
    faces.push_back(f7);
    faces.push_back(f8);
    faces.push_back(f9);
    faces.push_back(f10);
    faces.push_back(f11);
    faces.push_back(f12);
}
void Figure::icosahedron(){
    points.push_back(Vector3D::point(0,0, std::sqrt(5)/2));
    for(int l = 2; l < 7; l++){
        points.push_back(Vector3D::point(std::cos(2*M_PI*(l-2)/5), std::sin(2*M_PI*(l-2)/5), 0.5));
    }
    for(int l = 7; l < 12; l++){
        points.push_back(Vector3D::point(std::cos((M_PI/5)+((l-7)*2*M_PI/5)), std::sin((M_PI/5)+((l-7)*2*M_PI/5)), -0.5));
    }
    points.push_back(Vector3D::point(0,0, -std::sqrt(5)/2));

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

    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
    faces.push_back(f7);
    faces.push_back(f8);
    faces.push_back(f9);
    faces.push_back(f10);
    faces.push_back(f11);
    faces.push_back(f12);
    faces.push_back(f13);
    faces.push_back(f14);
    faces.push_back(f15);
    faces.push_back(f16);
    faces.push_back(f17);
    faces.push_back(f18);
    faces.push_back(f19);
    faces.push_back(f20);
}

Figure Figure::parseObj(const std::string &src, std::vector<Texture *> &vector, Matrix eyeTransf, bool normalOn) {
    // Get the object
    obj::OBJFile obj_parser;
    std::ifstream input_stream(src);
    input_stream >> obj_parser;
    input_stream.close();
    obj::ObjectGroup object = obj_parser.get_object();

    // Get mtl file name
    std::string mtl_file = object.get_mtllib_file_name();
    // Get material name
    std::string material_name = object.get_mtl_name();
    int map_Ka = -1;
    int map_Kd = -1;
    int map_Ks = -1;
    if(material_name != ""){
        // Parse mtl file
        obj::MTLLibrary library;

        std::ifstream mtl_input_stream(mtl_file);
        mtl_input_stream >> library;
        mtl_input_stream.close();

        // TODO: do something met reflections
        // Gebruik naam van materiaaleigenschap als “entry”naam
        std::vector<double> amb = library[material_name]["Ka"].as_double_tuple_or_die();

        double m = library[material_name]["Ns"].as_double_or_default(1.0);

        // Read all textures
        std::string tex_map;
        if(library[material_name]["map_Ka"].as_string_if_exists(tex_map)){
            // Create texture and add the image
            Texture* texture = new Texture;
            texture->number = vector.size();
            texture->a = Vector3D::vector(0, 0, 0);
            texture->b = Vector3D::vector(0, 0, 0);
            texture->p = Vector3D::vector(0, 0, 0);
            img::EasyImage* textureImage = new img::EasyImage;
            std::ifstream fin(tex_map);
            fin >> *textureImage;
            fin.close();
            texture->image = textureImage;
            vector.push_back(texture);
            // Add texture number to face
            map_Ka = texture->number;
        }
        if(library[material_name]["map_Kd"].as_string_if_exists(tex_map)){
            // Create texture and add the image
            Texture* texture = new Texture;
            texture->number = vector.size();
            texture->a = Vector3D::vector(0, 0, 0);
            texture->b = Vector3D::vector(0, 0, 0);
            texture->p = Vector3D::vector(0, 0, 0);
            img::EasyImage* textureImage = new img::EasyImage;
            std::ifstream fin(tex_map);
            fin >> *textureImage;
            fin.close();
            texture->image = textureImage;
            vector.push_back(texture);
            // Add texture number to face
            map_Kd = texture->number;
        }
        if(library[material_name]["map_Ks"].as_string_if_exists(tex_map)){
            // Create texture and add the image
            Texture* texture = new Texture;
            texture->number = vector.size();
            texture->a = Vector3D::vector(0, 0, 0);
            texture->b = Vector3D::vector(0, 0, 0);
            texture->p = Vector3D::vector(0, 0, 0);
            img::EasyImage* textureImage = new img::EasyImage;
            std::ifstream fin(tex_map);
            fin >> *textureImage;
            fin.close();
            texture->image = textureImage;
            vector.push_back(texture);
            // Add texture number to face
            map_Ks = texture->number;
        }
    }

    // Create new figure
    Figure figure;

    // Get points coordinates
    std::vector<std::vector<double>> ObjPoints =  object.get_vertexes();
    // Set points
    for(std::vector<double> coos: ObjPoints){
        Vector3D point = Vector3D::point(coos[0], coos[1], coos[2]);
        figure.points.push_back(point);
    }

    // Get polygons (aka super faces)
    std::vector<obj::Polygon> polygons = object.get_polygons();
    // Get all uv coordinates
    std::vector<std::vector<double>> allUVs = object.get_texture_coordinates();
    // Get all normal vectors' coordinates
    std::vector<std::vector<double>> allNorms = object.get_vertex_normals();

    // Transform polygons in conventional faces
    for(obj::Polygon polygon: polygons){
        std::vector<int> indexes = polygon.get_indexes();
        for(int &index:indexes){
            index--;
        }
        // TODO: meerdere texturen, reflecties
        Face face = Face(indexes);
        // Get uv's of the face
        if(polygon.has_texture_indexes()){
            face.map_Ka = map_Ka;
            face.map_Kd = map_Kd;
            face.map_Ks = map_Ks;
            std::vector<int> uv_indexes = polygon.get_texture_indexes();
            for(int &index:uv_indexes){
                index--;
                face.uv.push_back(allUVs[index]);
            }
        }
        if(normalOn){
            // Get normal vector's of the face
            if(polygon.has_normal_indexes()){
                std::vector<int> norm_indexes = polygon.get_normal_indexes();
                for(int &index:norm_indexes){
                    // TODO: Or is it better to do with new?
                    index--;
                    Vector3D* pointNorm = new Vector3D;
                    *pointNorm = Vector3D::vector(allNorms[index][0], allNorms[index][1], allNorms[index][2]);
                    *pointNorm *= eyeTransf;
                    face.norm.push_back(pointNorm);
                }
            }
        }
        figure.faces.push_back(face);
    }

    return figure;
}
void Figure::cylinderZonderVlakken(int n, double height){
    for (int ind = 0; ind < n; ind++) {
        points.push_back(Vector3D::point(std::cos(2 * M_PI * ind / n), std::sin(2 * M_PI * ind / n), 0));
        faces.push_back(Face({ind, (ind + 1) % n, n + (ind + 1) % (n), n + ind}));
    }
    for (int ind = 0; ind < n; ind++) {
        points.push_back(
                Vector3D::point(std::cos(2 * M_PI * ind / n), std::sin(2 * M_PI * ind / n), height));
    }
}
void Figure::sphere(int n){
    icosahedron();
    while (n > 0) {
        std::vector<Face> newFaces;
        std::vector<Vector3D> newPoints;
        // Collect newFaces and newPoints
        for (auto face: faces) {
            Vector3D A = points[face.point_indexes[0]];
            Vector3D B = points[face.point_indexes[1]];
            Vector3D C = points[face.point_indexes[2]];
            Vector3D D = utils::findMiddle(A, B);
            // ??? Swap E and F
            Vector3D E = utils::findMiddle(B, C);
            Vector3D F = utils::findMiddle(C, A);
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
            for (const auto& f: {f1, f2, f3, f4}) {
                newFaces.push_back(f);
            }
        }
        // Use the finalTrans matrix

        faces = newFaces;
        points = newPoints;
        n--;
    }
    // Herschaal alle punten
    for (auto &p: points) herschaalPuntenBal(p);
}
void Figure::herschaalPuntenBal(Vector3D& punt){
    double r = std::sqrt(std::pow(punt.x, 2) + std::pow(punt.y, 2) + std::pow(punt.z, 2));
    punt.x /= r;
    punt.y /= r;
    punt.z /= r;
}

void Figure::buckyBall() {
    icosahedron();
}

void Figure::mengerSponge(){
    cube();
}
