/*
 * easy_image.cc
 * Copyright (C) 2011  Daniel van den Akker
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "easy_image.h"
#include "vector3d.h"
#include <algorithm>
#include <assert.h>
#include "Light.h"
#include "Figure.h"
#include "Face.h"
#include "utils.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include "Line2D.h"
#include "Texture.h"

#ifndef le32toh
#define le32toh(x) (x)
#endif


namespace
{
	//structs borrowed from wikipedia's article on the BMP file format
	struct bmpfile_magic
	{
			uint8_t magic[2];
	};

	struct bmpfile_header
	{
			uint32_t file_size;
			uint16_t reserved_1;
			uint16_t reserved_2;
			uint32_t bmp_offset;
	};
	struct bmp_header
	{
			uint32_t header_size;
			int32_t width;
			int32_t height;
			uint16_t nplanes;
			uint16_t bits_per_pixel;
			uint32_t compress_type;
			uint32_t pixel_size;
			int32_t hres;
			int32_t vres;
			uint32_t ncolors;
			uint32_t nimpcolors;
	};
	//copy-pasted from lparser.cc to allow these classes to be used independently from each other
	class enable_exceptions
	{
		private:
			std::ios& ios;
			std::ios::iostate state;
		public:
			enable_exceptions(std::ios& an_ios, std::ios::iostate exceptions) :
				ios(an_ios)
			{
				state = ios.exceptions();
				ios.exceptions(exceptions);
			}
			~enable_exceptions()
			{
				ios.exceptions(state);
			}
	};
	//helper function to convert a number (char, int, ...) to little endian
	//regardless of the endiannes of the system
	//more efficient machine-dependent functions exist, but this one is more portable
	template<typename T> T to_little_endian(T value)
	{
		//yes, unions must be used with caution, but this is a case in which a union is needed
		union
		{
				T t;
				uint8_t bytes[sizeof(T)];
		} temp_storage;

		for (uint8_t i = 0; i < sizeof(T); i++)
		{
			temp_storage.bytes[i] = value & 0xFF;
			value >>= 8;
		}
		return temp_storage.t;
	}

	template<typename T> T from_little_endian(T value)
	{
		//yes, unions must be used with caution, but this is a case in which a union is needed
		union
		{
				T t;
				uint8_t bytes[sizeof(T)];
		} temp_storage;
		temp_storage.t = value;
		T retVal = 0;

		for (uint8_t i = 0; i < sizeof(T); i++)
		{
			retVal = (retVal << 8) | temp_storage.bytes[sizeof(T) - i - 1];
		}
		return retVal;
	}

}
img::Color::Color() :
	blue(0), green(0), red(0)
{
}
img::Color::Color(uint8_t r, uint8_t g, uint8_t b) :
	blue(b), green(g), red(r)
{
}
img::Color::~Color()
{
}

img::UnsupportedFileTypeException::UnsupportedFileTypeException(std::string const& msg) :
	message(msg)
{
}
img::UnsupportedFileTypeException::UnsupportedFileTypeException(const UnsupportedFileTypeException &original)
: std::exception(original)
, message(original.message)
{
}
img::UnsupportedFileTypeException::~UnsupportedFileTypeException() throw ()
{
}
img::UnsupportedFileTypeException& img::UnsupportedFileTypeException::operator=(UnsupportedFileTypeException const& original)
{
	this->message = original.message;
	return *this;
}
const char* img::UnsupportedFileTypeException::what() const throw ()
{
	return message.c_str();
}

img::EasyImage::EasyImage() :
	width(0), height(0), bitmap()
{
}

img::EasyImage::EasyImage(unsigned int _width, unsigned int _height, Color color) :
	width(_width), height(_height), bitmap(width * height, color), buf(_width, _height)
{
}

img::EasyImage::EasyImage(EasyImage const& img) :
	width(img.width), height(img.height), bitmap(img.bitmap)
{
}

img::EasyImage::~EasyImage()
{
    std::cout << "The image is deleted" << std::endl;
	bitmap.clear();
}

img::EasyImage& img::EasyImage::operator=(img::EasyImage const& img)
{
	width = img.width;
	height = img.height;
	bitmap.assign(img.bitmap.begin(),img.bitmap.end());
	return (*this);
}

unsigned int img::EasyImage::get_width() const
{
	return width;
}

unsigned int img::EasyImage::get_height() const
{
	return height;
}

void img::EasyImage::clear(Color color)
{
	for (std::vector<Color>::iterator i = bitmap.begin(); i != bitmap.end(); i++)
	{
		*i = color;
	}
}

img::Color& img::EasyImage::operator()(unsigned int x, unsigned int y)
{
	assert(x < this->width);
	assert(y < this->height);
	return bitmap.at(x * height + y);
}

img::Color const& img::EasyImage::operator()(unsigned int x, unsigned int y) const
{
	assert(x < this->width);
	assert(y < this->height);
	return bitmap.at(x * height + y);
}

void img::EasyImage::draw_line(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, Color color)
{
	if (x0 >= this->width || y0 >= this->height || x1 >= this->width || y1 > this->height) {
        std::stringstream ss;
        ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
        << this->width << " and height " << this->height;
        throw std::runtime_error(ss.str());
    }
	if (x0 == x1)
	{
		//special case for x0 == x1
		for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
		{
			(*this)(x0, i) = color;
		}
	}
	else if (y0 == y1)
	{
		//special case for y0 == y1
		for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
		{
			(*this)(i, y0) = color;
		}
	}
	else
	{
		if (x0 > x1)
		{
			//flip points if x1>x0: we want x0 to have the lowest value
			std::swap(x0, x1);
			std::swap(y0, y1);
		}
		double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
		if (-1.0 <= m && m <= 1.0)
		{
			for (unsigned int i = 0; i <= (x1 - x0); i++)
			{
				(*this)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
			}
		}
		else if (m > 1.0)
		{
			for (unsigned int i = 0; i <= (y1 - y0); i++)
			{
				(*this)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
			}
		}
		else if (m < -1.0)
		{
			for (unsigned int i = 0; i <= (y0 - y1); i++)
			{
				(*this)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
			}
		}
	}
}
void img::EasyImage::draw_zbuf_line(double x0, double y0, double z0, double x1, double y1, double z1, img::Color &color)
{
    if (x0 >= this->width || y0 >= this->height || x1 >= this->width || y1 > this->height) {
        std::stringstream ss;
        ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
           << this->width << " and height " << this->height;
        throw std::runtime_error(ss.str());
    }
    if(z0 == 0 || z1 == 0) std::cerr << "What I have to do if z0 pr z1 = 0?... Line: " << __LINE__ << std::endl;
    if (x0 == x1)
    {
        //special case for x0 == x1
        // Vert lijn
        for (unsigned int i = (unsigned int)std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            double een_z;
            if (y0 != y1) {
                een_z = ((std::max(y0,y1) - i) / (std::max(y0,y1) - std::min(y0,y1))) / z0 + (1 - ((std::max(y0,y1) - i) / (std::max(y0,y1) - std::min(y0,y1)))) / z1;
            } else{
                een_z = ((std::max(y0,y1) - i) + (1 - std::max(y0,y1) - i));
            }
            if (een_z < buf[i][static_cast<int>(std::lround(x0))]) {
                (*this)((unsigned int) x0, i) = color;
                buf[i][(unsigned int)(std::lround(x0))] = een_z;
            }
        }
    }
    else if (y0 == y1)
    {
        //special case for y0 == y1
        for (unsigned int i = (unsigned int)std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            double een_z;
            if (x0 != x1) {
                een_z = (((std::max(x0, x1) - i) / (std::max(x0, x1) - std::min(x0, x1))) / z0) + ((1 - ((std::max(x0, x1) - i) / (std::max(x0, x1) - std::min(x0, x1)))) / z1);
            } else {
                een_z = ((std::max(x0,x1) - i) + (1 - std::max(x0,x1) - i));
            }
            if (een_z < buf[(unsigned int)(std::lround(y0))][i]) {
                (*this)(i, (unsigned int) y0) = color;
                buf[(unsigned int)(std::lround(y0))][i] = een_z;
            }
        }
    }
    else
    {
        // Let's go from left to right!
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        // Calc rico
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                double een_z = (((x1 - x0 - i) / (x1 - x0)) / z0) + (1 - ((x1 - x0 - i) / (x1 - x0))) / z1;
                if (een_z < buf[(unsigned int) lround(y0 + m * i)][(unsigned int) lround(x0 + i)]) {
                    (*this)((unsigned int) x0 + i, (unsigned int) lround(y0 + m * i)) = color;
                    buf[(unsigned int) lround(y0 + m * i)][(unsigned int) lround(x0 + i)] = een_z;
                }
            }
        }
        else if (m > 1.0)
        {
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                double een_z = (((y1 - y0 - i) / (y1 - y0)) / z0) + (1 - ((y1 - y0 - i) / (y1 - y0))) / z1;
                if (een_z < buf[(unsigned int) y0 + i][(unsigned int) lround(x0 + (i / m))]) {
                    (*this)((unsigned int) lround(x0 + (i / m)), (unsigned int) y0 + i) = color;
                    buf[(unsigned int) y0 + i][(unsigned int) lround(x0 + (i / m))] = een_z;
                }
            }
        }
        else if (m < -1.0)
        {
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                double een_z = (((y0 - y1 - i) / (y0 - y1)) / z0) + (1 - ((y0 - y1 - i) / (y0 - y1))) / z1;
                if (een_z < buf[(int)y0 - i][(unsigned int) lround(x0 - (i / m))]) {
                    (*this)((unsigned int) lround(x0 - (i / m)), (unsigned int)y0 - i) = color;
                    buf[(int)y0 - i][(unsigned int) lround(x0 - (i / m))] = een_z;
                }
            }
        }
    }
}
double projecteerCo(double x, double z, double d, double dx){
    return ((d*x)/(-z)) + dx;
}
inline int roundToInt(double d)
{
    return static_cast<int>(std::round(d));
}
bool img::EasyImage::berekenZEnSteekInBuffer(unsigned int x, unsigned int y, double xg, double yg, double zg, double dzdx, double dzdy){
    /*
     * Geeft true terug als het punt getekend moet worden
     */
    double een_op_z = zg + (x - xg)*dzdx + (y - yg)*dzdy;
    if(buf[y][x] > een_op_z){
        buf[y][x] = een_op_z;
        return true;
    }
    return false;
}
void img::EasyImage::draw_zbuf_triag_line(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, img::Color &color, double xg, double yg, double zg, double dzdx, double dzdy){
    if (x0 >= this->width || y0 >= this->height || x1 >= this->width || y1 > this->height) {
        std::stringstream ss;
        ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
           << this->width << " and height " << this->height;
        throw std::runtime_error(ss.str());
    }
    if (x0 == x1)
    {
        //special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            if(berekenZEnSteekInBuffer(x0, i, xg, yg, zg, dzdx, dzdy)){
                (*this)(x0, i) = color;
            }
        }
    }
    else if (y0 == y1)
    {
        //special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            if(berekenZEnSteekInBuffer(i, y0, xg, yg, zg, dzdx, dzdy)) {
                (*this)(i, y0) = color;
            }
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                if(berekenZEnSteekInBuffer(x0 + i, (unsigned int) round(y0 + m * i), xg, yg, zg, dzdx, dzdy)) {
                    (*this)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                }
            }
        }
        else if (m > 1.0)
        {
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                if(berekenZEnSteekInBuffer((unsigned int) round(x0 + (i / m)), y0 + i, xg, yg, zg, dzdx, dzdy)) {
                    (*this)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                }
            }
        }
        else if (m < -1.0)
        {
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                if(berekenZEnSteekInBuffer((unsigned int) round(x0 - (i / m)), y0 - i, xg, yg, zg, dzdx, dzdy)) {
                    (*this)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                }
            }
        }
    }
}
#include <vector>
#include <memory>

// TODO: embed this code to utils

// Function to perform deep copy of MyClass
std::unique_ptr<Figure> DeepCopy(const Figure& obj) {
    return std::make_unique<Figure>(obj);
}

// Function to make deep copy of vector of pointers
std::vector<std::unique_ptr<Figure>> DeepCopyVector(const std::vector<std::unique_ptr<Figure>>& vec) {
    std::vector<std::unique_ptr<Figure>> copiedVec;
    copiedVec.reserve(vec.size());

    for (const auto& ptr : vec) {
        copiedVec.push_back(DeepCopy(*ptr));
    }

    return copiedVec;
}

void img::EasyImage::fillShadowBuffers(const std::vector<Figure>& figures, std::vector<Light*> &lights, double d, double dx, double dy) const{
    for(Light* light: lights){
        // Do copy of all figures
        std::vector<Figure> lightFigures = std::vector<Figure>(figures);
        utils::applyTransformation(lightFigures, light->eye);
        Lines2D light_toDraw;
        utils::doProjection(lightFigures, light_toDraw);
        // Set dx, dy and d
        utils::getDxDyD(light_toDraw, light->shadowMask.size(), light->dx, light->dy, light->d);
        for(const auto& fig:lightFigures){
            for(const auto& fac: fig.faces){
                int A = fac.point_indexes[0];
                int B = fac.point_indexes[1];
                int C = fac.point_indexes[2];
                shadow_zbuf_triag(fig.points[A], fig.points[B], fig.points[C], light);
                }
            }
    }
}
void img::EasyImage::shadow_zbuf_triag(const Vector3D& A, const Vector3D& B, const Vector3D& C, Light* light) const{
    // Bereken points
    Point2D a = Point2D(projecteerCo(A.x, A.z, light->d, light->dx), projecteerCo(A.y, A.z, light->d, light->dy));
    Point2D b = Point2D(projecteerCo(B.x, B.z, light->d, light->dx), projecteerCo(B.y, B.z, light->d, light->dy));
    Point2D c = Point2D(projecteerCo(C.x, C.z, light->d, light->dx), projecteerCo(C.y, C.z, light->d, light->dy));
    // Pixels tot driehoek
    double y_max = std::max(a.y, std::max(b.y, c.y));
    double y_min = std::min(a.y, std::min(b.y, c.y));
    // Bereken xg, yg, zg
    double xg = (a.x + b.x + c.x)/3.0;
    double yg = (a.y + b.y + c.y)/3.0;
    double zg = (1.0/(3.0*A.z)) + (1.0/(3.0*B.z)) + (1.0/(3.0*C.z));
    for(int y_i = roundToInt(y_min+0.5); y_i <= roundToInt(y_max-0.5); y_i++){
        //double y_i = static_cast<double> (y_ik);
        double x_lab = std::numeric_limits<double>::infinity();
        double x_lac = std::numeric_limits<double>::infinity();
        double x_lbc = std::numeric_limits<double>::infinity();
        double x_rab = -std::numeric_limits<double>::infinity();
        double x_rac = -std::numeric_limits<double>::infinity();
        double x_rbc = -std::numeric_limits<double>::infinity();
        // Splits in 3 gevallen
        // AB
        if((y_i - a.y)*(y_i - b.y) <= 0 && a.y != b.y){
            double x_i = b.x + ((a.x - b.x)*(y_i - b.y)/(a.y - b.y));
            x_lab = x_i;
            x_rab = x_i;
        }
        // BC
        if((y_i - b.y)*(y_i - c.y) <= 0 && b.y != c.y){
            double x_i = c.x + ((b.x - c.x)*(y_i - c.y)/(b.y - c.y));
            x_lbc = x_i;
            x_rbc = x_i;
        }
        // AC
        if((y_i - a.y)*(y_i - c.y) <= 0 && a.y != c.y){
            double x_i = c.x + ((a.x - c.x)*(y_i - c.y)/(a.y - c.y));
            x_lac = x_i;
            x_rac = x_i;
        }
        // Bereken u-vector
        Vector3D u = B - A;
        // Bereken v-vector
        Vector3D v = C - A;
        // Bereken w-vector (Vector3D::cross is voor watjes)
        Vector3D w = Vector3D::vector(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y-u.y*v.x);

        // Bereken k
        double k = w.x*A.x + w.y*A.y + w.z*A.z;
        // if (k <= 0) return; // Backface culling; 34 sec ipv 43
        if (k == 0) return;
        // Bereken dzdx
        double dzdx = w.x/(-light->d*k);
        // Bereken dzdy
        double dzdy = w.y/(-light->d*k);
        // "Draw" everything line per line
        if(std::round(std::min(x_lab, std::min(x_lac, x_lbc))) != std::numeric_limits<double>::infinity() &&
           std::round(std::max(x_rab, std::max(x_rac, x_rbc))) != -std::numeric_limits<double>::infinity()) {
            int x_r = roundToInt(std::min(x_lab, std::min(x_lac, x_lbc)) + 0.5);
            int x_l = roundToInt(std::max(x_rab, std::max(x_rac, x_rbc)) - 0.5);
            // Geef ook dzdx en dzdy mee
            double Z = zg + (x_r-xg)*dzdx + (y_i-yg)*dzdy;

            for (int x = x_r; x <= x_l; x++) {
                // TODO: wat als shadowMask is kleiner dan deze waarden?
                if(light->shadowMask[y_i][x] > Z){
                    light->shadowMask[y_i][x] = Z;
                    if(Z > light->shadowMask.max) light->shadowMask.max = Z;
                    if(Z < light->shadowMask.min) light->shadowMask.min = Z;
                }else{
                    //std::cout << "";
                }
                Z += dzdx;
            }

        }
    }
}
img::Color vectorToColor2(std::vector<double> kleur){
    img::Color to_return = img::Color(kleur[0]*255, kleur[1]*255, kleur[2]*255);
    return to_return;
}
void img::EasyImage::draw_zbuf_triag(const Vector3D& A, const Vector3D& B, const Vector3D& C, double d, double dx, double dy, std::vector<double>& fullAmbientRef, std::vector<double>& diffuseRef, std::vector<double>& specularRef, double& refCoeff, const std::vector<Light*>& lights, const Vector3D& eyeCamera, const Matrix &eyeTransf, bool shadowOn, const std::vector<Texture*>& textures, const std::vector<Texture*>& fig_cubemap_textures, const std::vector<std::vector<double> >& uv, std::vector<Vector3D* > norm, const double& cubeSize){
    // Bereken points
    Point2D a = Point2D(projecteerCo(A.x, A.z, d, dx), projecteerCo(A.y, A.z, d, dy));
    Point2D b = Point2D(projecteerCo(B.x, B.z, d, dx), projecteerCo(B.y, B.z, d, dy));
    Point2D c = Point2D(projecteerCo(C.x, C.z, d, dx), projecteerCo(C.y, C.z, d, dy));
    // Pixels tot driehoek
    double y_max = std::max(a.y, std::max(b.y, c.y));
    double y_min = std::min(a.y, std::min(b.y, c.y));
    // Bereken xg, yg, zg
    double xg = (a.x + b.x + c.x)/3.0;
    double yg = (a.y + b.y + c.y)/3.0;
    double zg = (1.0/(3.0*A.z)) + (1.0/(3.0*B.z)) + (1.0/(3.0*C.z));
    for(int y_i = roundToInt(y_min+0.5); y_i <= roundToInt(y_max-0.5); y_i++){
        //double y_i = static_cast<double> (y_ik);
        double x_lab = std::numeric_limits<double>::infinity();
        double x_lac = std::numeric_limits<double>::infinity();
        double x_lbc = std::numeric_limits<double>::infinity();
        double x_rab = -std::numeric_limits<double>::infinity();
        double x_rac = -std::numeric_limits<double>::infinity();
        double x_rbc = -std::numeric_limits<double>::infinity();
        // Splits in 3 gevallen
        // AB
        if((y_i - a.y)*(y_i - b.y) <= 0 && a.y != b.y){
            double x_i = b.x + ((a.x - b.x)*(y_i - b.y)/(a.y - b.y));
            x_lab = x_i;
            x_rab = x_i;
        }
        // BC
        if((y_i - b.y)*(y_i - c.y) <= 0 && b.y != c.y){
            double x_i = c.x + ((b.x - c.x)*(y_i - c.y)/(b.y - c.y));
            x_lbc = x_i;
            x_rbc = x_i;
        }
        // AC
        if((y_i - a.y)*(y_i - c.y) <= 0 && a.y != c.y){
            double x_i = c.x + ((a.x - c.x)*(y_i - c.y)/(a.y - c.y));
            x_lac = x_i;
            x_rac = x_i;
        }
        // Bereken u-vector
        Vector3D u = B - A;
        // Bereken v-vector
        Vector3D v = C - A;
        // Bereken w-vector (Vector3D::cross is voor watjes)
        Vector3D w = Vector3D::vector(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y-u.y*v.x);

        // Bereken n
        Vector3D n = w/w.length();
        n.normalise();
        // Bereken totale diffuse en specular component van alle lichten
        std::vector<double> fullDifCo = {0,0,0};
        std::vector<double> fullSpecCo = {0,0,0};
        for(Light* light:lights){
            // Go further if some components are not zero
            if(light->ldVector.length() != 0 && light->inf){
                // Bereken l
                Vector3D l = -light->ldVector / light->ldVector.length();
                l.normalise();
                double cos_alpha = n.x * l.x + n.y * l.y + n.z * l.z;
                if(cos_alpha > 0) {
                    fullDifCo[0] += diffuseRef[0] * light->diffuseLight[0] * cos_alpha;
                    fullDifCo[1] += diffuseRef[1] * light->diffuseLight[1] * cos_alpha;
                    fullDifCo[2] += diffuseRef[2] * light->diffuseLight[2] * cos_alpha;
                }
            }
        }

        // Bereken k
        double k = w.x*A.x + w.y*A.y + w.z*A.z;
        // if (k <= 0) return; // Backface culling; 34 sec ipv 43
        if (k == 0) return;
        // Bereken dzdx
        double dzdx = w.x/(-d*k);
        // Bereken dzdy
        double dzdy = w.y/(-d*k);
        // Draw everything line per line
        if(std::round(std::min(x_lab, std::min(x_lac, x_lbc))) != std::numeric_limits<double>::infinity() &&
                    std::round(std::max(x_rab, std::max(x_rac, x_rbc))) != -std::numeric_limits<double>::infinity()) {
            int x_r = roundToInt(std::min(x_lab, std::min(x_lac, x_lbc)) + 0.5);
            int x_l = roundToInt(std::max(x_rab, std::max(x_rac, x_rbc)) - 0.5);
            // Geef ook dzdx en dzdy mee
            double Z = zg + (x_r-xg)*dzdx + (y_i-yg)*dzdy;

            for (int x = x_r; x <= x_l; x++) {
                if(x < 0 || y_i < 0) continue;
                if(y_i > buf.size()) continue;
                if(x > buf[y_i].size()) continue;
                if(buf[y_i][x] > Z){
                    std::vector<double> resulted_color = {fullAmbientRef[0] + fullDifCo[0] + fullSpecCo[0], fullAmbientRef[1] + fullDifCo[1] + fullSpecCo[1],
                                                          fullAmbientRef[2] + fullDifCo[2] + fullSpecCo[2]};
                    double x_e = -(x - dx)/(d*Z);
                    double y_e = -(y_i - dy)/(d*Z);
                    for(Light* light:lights){
                        // Bepaal positie in eye co
                        Vector3D eyeCo = Vector3D::vector(- (x - dx)/(Z*d), - (y_i - dy)/(Z*d), 1/Z);
                        // Specular inf light
                        if(light->inf && !(light->specularLight[0] == 0 && light->specularLight[1] == 0 && light->specularLight[2] == 0)){
                            Vector3D l = -light->ldVector / light->ldVector.length();
                            l.normalise();
                            double cos_alpha = n.x * l.x + n.y * l.y + n.z * l.z;
                            Vector3D r = 2*n*cos_alpha - l;
                            r.normalise();
                            Vector3D cameraVector = eyeCamera - eyeCo;
                            cameraVector.normalise();
                            double cos_beta = Vector3D::dot(r, cameraVector);
                            if (cos_beta > 0) {
                                resulted_color[0] +=
                                        specularRef[0] * light->specularLight[0] * std::pow(cos_beta, refCoeff);
                                resulted_color[1] +=
                                        specularRef[1] * light->specularLight[1] * std::pow(cos_beta, refCoeff);
                                resulted_color[2] +=
                                        specularRef[2] * light->specularLight[2] * std::pow(cos_beta, refCoeff);
                            }
                        }
                        // Point light
                        if(shadowOn){
                            if(!light->pointIsVisible(x, y_i, 1/Z, d, dx, dy, eyeTransf)) continue;
                        }
                        if(light->location.length() != 0 && !light->inf){
                            // Bereken l
                            Vector3D l = light->location - eyeCo;
                            l.normalise();
                            if(!norm.empty()){
                                n = Face::getNorm(A, B, C, norm[0], norm[1], norm[2], Vector3D::point(x_e, y_e, 1/Z));
                                n.normalise();
                            }
                            double cos_alpha = n.x * l.x + n.y * l.y + n.z * l.z;
                            Vector3D r = 2*n*cos_alpha - l;
                            r.normalise();
                            Vector3D cameraVector = eyeCamera - eyeCo;
                            cameraVector.normalise();
                            double cos_beta = Vector3D::dot(r, cameraVector);
                            if (cos_beta > 0) {
                                resulted_color[0] +=
                                        specularRef[0] * light->specularLight[0] * std::pow(cos_beta, refCoeff);
                                resulted_color[1] +=
                                        specularRef[1] * light->specularLight[1] * std::pow(cos_beta, refCoeff);
                                resulted_color[2] +=
                                        specularRef[2] * light->specularLight[2] * std::pow(cos_beta, refCoeff);
                            }
                            if(cos_alpha > 0 && cos_alpha > std::cos(light->spotAngle)) {
                                resulted_color[0] += diffuseRef[0] * light->diffuseLight[0] *
                                                     (1 - ((1 - cos_alpha) / (1 - std::cos(light->spotAngle))));
                                resulted_color[1] += diffuseRef[1] * light->diffuseLight[1] *
                                                     (1 - ((1 - cos_alpha) / (1 - std::cos(light->spotAngle))));
                                resulted_color[2] += diffuseRef[2] * light->diffuseLight[2] *
                                                     (1 - ((1 - cos_alpha) / (1 - std::cos(light->spotAngle))));
                            }
                        }
                    }
                    // Textuurcoordinaten of textuur op een willekeurige vlak
                    if(!textures.empty()){
                        // TODO: nog andere texturen
                        // Get world coos
                        if(!uv.empty()){
                            // Build matrix
                            Matrix to_inverse;
                            to_inverse(1,1) =  A.x;
                            to_inverse(1,2) =  A.y;
                            to_inverse(1,3) =  A.z;
                            to_inverse(2,1) =  B.x;
                            to_inverse(2,2) =  B.y;
                            to_inverse(2,3) =  B.z;
                            to_inverse(3,1) =  C.x;
                            to_inverse(3,2) =  C.y;
                            to_inverse(3,3) =  C.z;
                            Matrix inverted = Matrix::inv(to_inverse);
                            // Get 1-p-q, p, q
                            Vector3D coefs = Vector3D::point(x_e, y_e, 1/Z)*inverted;
                            double texture_u = coefs.x*uv[0][0] + coefs.y*uv[1][0]+ coefs.z*uv[2][0];
                            double texture_v = coefs.x*uv[0][1] + coefs.y*uv[1][1]+ coefs.z*uv[2][1];
                            int u_t = std::abs(roundToInt((textures[0]->image->width-1)*(std::fmod(texture_u,1))));
                            int v_t = std::abs(roundToInt((textures[0]->image->height-1)*(std::fmod(texture_v, 1))));
                            Color textureColor = (*(textures[0]->image))(u_t, v_t);
                            std::vector<double> texturefullAmbientRef = {0,0,0};
                            resulted_color[0] -= fullAmbientRef[0];
                            resulted_color[1] -= fullAmbientRef[1];
                            resulted_color[2] -= fullAmbientRef[2];
                            texturefullAmbientRef[0] += textureColor.red/255.0;
                            texturefullAmbientRef[1] += textureColor.green/255.0;
                            texturefullAmbientRef[2] += textureColor.blue/255.0;
                            for (Light *light: lights) {
                                resulted_color[0] += texturefullAmbientRef[0] * light->ambientLight[0];
                                resulted_color[1] += texturefullAmbientRef[1] * light->ambientLight[1];
                                resulted_color[2] += texturefullAmbientRef[2] * light->ambientLight[2];
                            }
                        }else{
                            Vector3D worldCo = Vector3D::point(x_e, y_e, 1/Z)*Matrix::inv(eyeTransf);
                            std::vector<double> texturefullAmbientRef = {0,0,0};
                            bool black = false;
                            for(Texture* texture: textures){
    //                            // Get uvm
                                Vector3D uvw = texture->getUVW(worldCo);
                                // Get color with uvm
                                if(uvw.x < 0 || uvw.x > 1 || uvw.y < 0 || uvw.y > 1) continue;
                                uvw.x = uvw.x;
                                uvw.y = uvw.y;
                                int u_t = roundToInt((texture->image->width-1)*uvw.x);
                                int v_t = roundToInt((texture->image->height-1)*uvw.y);
                                Color textureColor = (*(texture->image))(u_t, v_t);
                                // Remove it or not: that is the question
                                if(textureColor.red == 0 && textureColor.green == 0 && textureColor.blue == 0) {
                                    black = true;
                                }
                                texturefullAmbientRef[0] += textureColor.red/255.0;
                                texturefullAmbientRef[1] += textureColor.green/255.0;
                                texturefullAmbientRef[2] += textureColor.blue/255.0;
                            }
                            if(!(texturefullAmbientRef[0] == 0 && texturefullAmbientRef[1] == 0 && texturefullAmbientRef[2] == 0) || black) {
                                // Remove color gotten from ambient light
                                resulted_color[0] -= fullAmbientRef[0];
                                resulted_color[1] -= fullAmbientRef[1];
                                resulted_color[2] -= fullAmbientRef[2];
    //                         Multiply with lights
                                for (Light *light: lights) {
                                    resulted_color[0] += texturefullAmbientRef[0] * light->ambientLight[0];
                                    resulted_color[1] += texturefullAmbientRef[1] * light->ambientLight[1];
                                    resulted_color[2] += texturefullAmbientRef[2] * light->ambientLight[2];
                                }
                            }
                            black = false;
                        }
                    }
                    else if(cubeSize != -1){
                        Vector3D eyeCo = Vector3D::vector(- (x - dx)/(Z*d), - (y_i - dy)/(Z*d), 1/Z);
                        std::vector<double> ts;
                        double t_top = (cubeSize - eyeCo.y)/n.y;
                        ts.push_back(t_top);
                        double t_left = (-cubeSize - eyeCo.x)/n.x;
                        ts.push_back(t_left);
                        double t_front = (-cubeSize - eyeCo.z)/n.z;
                        ts.push_back(t_front);
                        double t_right = (cubeSize - eyeCo.x)/n.x;
                        ts.push_back(t_right);
                        double t_back = (cubeSize - eyeCo.z)/n.z;
                        ts.push_back(t_back);
                        double t_bottom = (-cubeSize - eyeCo.y)/n.y;
                        ts.push_back(t_bottom);
                        int cubemap_side_index = -1;
                        int iterative_index = -1;
                        double prev_t = *std::max_element(ts.begin(), ts.end());
                        for(double candidate_t: ts){
                            iterative_index++;
                            if(candidate_t < 0) continue;
                            if(candidate_t < prev_t) {
                                prev_t = candidate_t;
                                cubemap_side_index = iterative_index;
                            }
                        }
                        if(cubemap_side_index == -1) cubemap_side_index = 0;
                        Texture* texture = fig_cubemap_textures[cubemap_side_index];
                        //Vector3D worldCo = Vector3D::point(x_e, y_e, 1/Z)*Matrix::inv(eyeTransf);
                        std::vector<double> texturefullAmbientRef = {0,0,0};
                        bool black = false;
                        Vector3D uvw = texture->getUVW(eyeCo);
                        // Get color with uvm
                        int u_t = std::abs(roundToInt((texture->image->width-1)*(std::fmod(uvw.x,1.0))));
                        int v_t = std::abs(roundToInt((texture->image->height-1)*(std::fmod(uvw.y, 1.0))));
                        Color textureColor = (*(texture->image))(u_t, v_t);
                        // Remove it or not: that is the question
                        if(textureColor.red == 0 && textureColor.green == 0 && textureColor.blue == 0) {
                            black = true;
                        }
                        texturefullAmbientRef[0] += textureColor.red/255.0;
                        texturefullAmbientRef[1] += textureColor.green/255.0;
                        texturefullAmbientRef[2] += textureColor.blue/255.0;

                        if(!(texturefullAmbientRef[0] == 0 && texturefullAmbientRef[1] == 0 && texturefullAmbientRef[2] == 0) || black) {
                            // Remove color gotten from ambient light
                            resulted_color[0] -= fullAmbientRef[0];
                            resulted_color[1] -= fullAmbientRef[1];
                            resulted_color[2] -= fullAmbientRef[2];
                            //                         Multiply with lights
                            for (Light *light: lights) {
                                resulted_color[0] += texturefullAmbientRef[0] * light->ambientLight[0];
                                resulted_color[1] += texturefullAmbientRef[1] * light->ambientLight[1];
                                resulted_color[2] += texturefullAmbientRef[2] * light->ambientLight[2];
                            }
                        }
                        black = false;
                    }
                    if(resulted_color[0] > 1) resulted_color[0] = 1;
                    if(resulted_color[1] > 1) resulted_color[1] = 1;
                    if(resulted_color[2] > 1) resulted_color[2] = 1;

                    bitmap[x*height + y_i] = vectorToColor2(resulted_color);
                    buf[y_i][x] = Z;
                    //std::cout << "Pixel" <<std::endl;
                }else{
                    //std::cout << "";
                }
                Z += dzdx;
            }

        }
    }
}
std::ostream& img::operator<<(std::ostream& out, EasyImage const& image)
{

	//temporaryily enable exceptions on output stream
	enable_exceptions(out, std::ios::badbit | std::ios::failbit);
	//declare some struct-vars we're going to need:
	bmpfile_magic magic;
	bmpfile_header file_header;
	bmp_header header;
	uint8_t padding[] =
	{ 0, 0, 0, 0 };
	//calculate the total size of the pixel data
	unsigned int line_width = image.get_width() * 3; //3 bytes per pixel
	unsigned int line_padding = 0;
	if (line_width % 4 != 0)
	{
		line_padding = 4 - (line_width % 4);
	}
	//lines must be aligned to a multiple of 4 bytes
	line_width += line_padding;
	unsigned int pixel_size = image.get_height() * line_width;

	//start filling the headers
	magic.magic[0] = 'B';
	magic.magic[1] = 'M';

	file_header.file_size = to_little_endian(pixel_size + sizeof(file_header) + sizeof(header) + sizeof(magic));
	file_header.bmp_offset = to_little_endian(sizeof(file_header) + sizeof(header) + sizeof(magic));
	file_header.reserved_1 = 0;
	file_header.reserved_2 = 0;
	header.header_size = to_little_endian(sizeof(header));
	header.width = to_little_endian(image.get_width());
	header.height = to_little_endian(image.get_height());
	header.nplanes = to_little_endian(1);
	header.bits_per_pixel = to_little_endian(24);//3bytes or 24 bits per pixel
	header.compress_type = 0; //no compression
	header.pixel_size = pixel_size;
	header.hres = to_little_endian(11811); //11811 pixels/meter or 300dpi
	header.vres = to_little_endian(11811); //11811 pixels/meter or 300dpi
	header.ncolors = 0; //no color palette
	header.nimpcolors = 0;//no important colors

	//okay that should be all the header stuff: let's write it to the stream
	out.write((char*) &magic, sizeof(magic));
	out.write((char*) &file_header, sizeof(file_header));
	out.write((char*) &header, sizeof(header));

	//okay let's write the pixels themselves:
	//they are arranged left->right, bottom->top, b,g,r
	for (unsigned int i = 0; i < image.get_height(); i++)
	{
		//loop over all lines
		for (unsigned int j = 0; j < image.get_width(); j++)
		{
			//loop over all pixels in a line
			//we cast &color to char*. since the color fields are ordered blue,green,red they should be written automatically
			//in the right order
			out.write((char*) &image(j, i), 3 * sizeof(uint8_t));
		}
		if (line_padding > 0)
			out.write((char*) padding, line_padding);
	}
	//okay we should be done
	return out;
}
std::istream& img::operator>>(std::istream& in, EasyImage & image)
{
	enable_exceptions(in, std::ios::badbit | std::ios::failbit);
	//declare some struct-vars we're going to need
	bmpfile_magic magic;
	bmpfile_header file_header;
	bmp_header header;
	//a temp buffer for reading the padding at the end of each line
	uint8_t padding[] =
	{ 0, 0, 0, 0 };

	//read the headers && do some sanity checks
	in.read((char*) &magic, sizeof(magic));
	if (magic.magic[0] != 'B' || magic.magic[1] != 'M')
		throw UnsupportedFileTypeException("Could not parse BMP File: invalid magic header");
	in.read((char*) &file_header, sizeof(file_header));
	in.read((char*) &header, sizeof(header));
	if (le32toh(header.pixel_size) + le32toh(file_header.bmp_offset) != le32toh(file_header.file_size))
		throw UnsupportedFileTypeException("Could not parse BMP File: file size mismatch");
	if (le32toh(header.header_size) != sizeof(header))
		throw UnsupportedFileTypeException("Could not parse BMP File: Unsupported BITMAPV5HEADER size");
	if (le32toh(header.compress_type) != 0)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only uncompressed BMP files can be parsed");
	if (le32toh(header.nplanes) != 1)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only one plane should exist in the BMP file");
	if (le32toh(header.bits_per_pixel) != 24)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only 24bit/pixel BMP's are supported");
	//if height<0 -> read top to bottom instead of bottom to top
	bool invertedLines = from_little_endian(header.height) < 0;
	image.height = std::abs(from_little_endian(header.height));
	image.width = std::abs(from_little_endian(header.width));
	unsigned int line_padding = from_little_endian(header.pixel_size) / image.height - (3 * image.width);
	//re-initialize the image bitmap
	image.bitmap.clear();
	image.bitmap.assign(image.height * image.width, Color());
	//okay let's read the pixels themselves:
	//they are arranged left->right., bottom->top if height>0, top->bottom if height<0, b,g,r
	for (unsigned int i = 0; i < image.get_height(); i++)
	{
		//loop over all lines
		for (unsigned int j = 0; j < image.get_width(); j++)
		{
			//loop over all pixels in a line
			//we cast &color to char*. since the color fields are ordered blue,green,red, the data read should be written in the right variables
			if (invertedLines)
			{
				//store top-to-bottom
				in.read((char*) &image(j, image.height - 1 - i), 3 * sizeof(uint8_t));
			}
			else
			{
				//store bottom-to-top
				in.read((char*) &image(j, i), 3 * sizeof(uint8_t));
			}
		}
		if (line_padding > 0)
		{
			in.read((char*) padding, line_padding);
		}
	}
	//okay we're done
	return in;
}
