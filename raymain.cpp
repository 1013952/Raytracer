#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "matrices.h"
#include "bmphandling.h"
#include "geometry.h"

int main()
{

	const int width = 1024;
	const int height = 1024;
	
	
	const int screensize = 20;
	point3 screencent = {0, 5, 0};
	point3 right = {1, 0, 0};
	point3 down = {0, 0, 1};
	point3 forward = {0, 1, 0};

	double cameradistance = 20;
	point3 camera = screencent - forward * cameradistance;
 
	point3 c1 = {10, 200, 10};
	point3 c2 = {-3, 50, 0};
	point3 c3 = {5, 30, 0};
	point3 c4 = {0, -20, 0};
	colourt white = {0, 0.95, 0.95, 0.95};
	colourt black = {0, 0.05, 0.05, 0.05};
	colourt blue = {0, 0.2, 0.2, 0.95};
	colourt green = {0, 0.3, 0.95, 0.3};
	colourt red = {0, 0.95, 0.3, 0.3};
	colourt truered = {0, 0.8, 0.2, 0.2};
	colourt stone = {0, 0.3, 0.3, 0.3};

	std::vector<Shape*> shapes;
	shapes.push_back(new Sphere(c2, 10, white, 0));
	shapes.push_back(new Sphere(c3, 10, white, 0));
	//shapes.push_back(new Sphere(c1, 10, white, 0));
	shapes.push_back(new Plane(0, 0, 1, -10, stone));
	//shapes.push_back(new Plane(0, 0, -1, -200, stone));


	point3 pt = {50, -20, 30};
	point3 pt2 = {0, -30, 30};
	point3 pt3 = {-50, -20, 30};
	point3 sky = {0, 15, 50};
	std::vector<Light> lights;
	//Light light0(sky, white, true, 10);
	//Light light1(pt, white, false, 20);
	//Light light2(pt2, white, false, 20);
	//Light light3(pt3, white, false, 20);
	//lights.push_back(light0);
	//lights.push_back(light1);
	//lights.push_back(light2);
	//lights.push_back(light3);

	diffuseLight light4(pt, red, 25, 10, 20);
	light4.getuniques(lights);
	diffuseLight light5(pt3, blue, 25, 10, 20);
	light5.getuniques(lights);

	const colourt background = {0, 0, 0.825, 0.825};

	std::ofstream fout;
	fout.open("image.bmp", std::ios_base::out | std::ios_base::binary);
	int padSize = prepbmp(fout, width, height);

	unsigned char pad[3] = {0,0,0};
	double redpx;
	double greenpx;
	double bluepx;

	for(int i = 0; i < height; ++i){
		for(int j = 0; j < width; ++j)
		{
			point3 pixel = screencent + right * (- screensize/2 + screensize*(double) j/(double) width) + down * (- screensize/2 + screensize*(double) i/(double) height);
			
			redpx = 0; greenpx = 0; bluepx = 0;

			ray c = linetwopts(camera, pixel);
			//std::cout<<c.direction.x<<" "<<c.direction.y<<" "<<c.direction.z<<" "<<c.origin.x<<" "<<c.origin.y<<" "<<c.origin.z<<" "<<veclen(c.direction)<<std::endl;
			//system("pause");
			colourt temp, result;
			result = askallshapes(c, lights, shapes);
			redpx += result.red;
			bluepx += result.blue;
			greenpx += result.green;

			bmppixel(fout, redpx, greenpx, bluepx);
		}
		fout.write( (char*)pad, padSize );
	}

	return 0;
}