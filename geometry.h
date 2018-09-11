#ifndef _GEOMETRY
#define _GEOMETRY

#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "matrices.h"

const double epsilon = 0.001;

struct point3
{
	double x;
	double y;
	double z;
	point3 operator*(double lambda)
	{
		point3 b = {x*lambda, y*lambda, z*lambda};
		return b;
	}

	point3 operator*(point3 b)
	{
		point3 c = {x*b.x, y*b.y, z*b.z};
		return c;
	}

	point3 operator/(double lambda)
	{
		if(lambda!=0){
			point3 b = {x/lambda, y/lambda, z/lambda};
			return b;}
		else
			return *this;
	}

	point3 operator-(const point3 b){
		point3 c = {
			x - b.x, y - b.y, z - b.z
		};
		return c;
	}

	point3 operator+(const point3 b){
		point3 c = {
			x + b.x, y + b.y, z + b.z
		};
		return c;
	}

};



double scalarprod(point3 a, point3 b){
		return a.x*b.x + a.y*b.y + a.z*b.z;
	}


point3 crossprod(point3 a, point3 b){
	double x = a.y * b.z - a.z * b.y;
	double y = a.z * b.x - a.x * b.z;
	double z = a.x * b.y - a.y * b.x;
	point3 res = {x,y,z};
	return res;
}
double crossprodmod(point3 a, point3 b){
	double x = a.y * b.z - a.z * b.y;
	double y = a.z * b.x - a.x * b.z;
	double z = a.x * b.y - a.y * b.x;
	return sqrt(x*x + y*y + z*z);
}

double veclen(point3 a){
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

double pointdist(point3 a, point3 b){
	return veclen(b-a);
}


struct ray{ //x = x1+at, y = y1 + bt, z = z1 + ct
	point3 direction;
	point3 origin;
};

ray linetwopts(point3 a, point3 b)
{
	point3 c = (b - a)/veclen(b - a);
	ray res;
	res.direction = c;
	res.origin = a;
	return res;
}

struct plane{ //ax+by+cz=1;
	double a;
	double b;
	double c;
};

plane planefrom3pts(point3 a, point3 b, point3 c){
	plane res;
	mat34 newmatrix = {
			{{a.x, a.y, a.z, 1.0}, {b.x, b.y, b.z, 1.0}, {c.x, c.y, c.z, 1.0}}
		};
	newmatrix.rref();
	//TODO: if plane is not well-defined
	res.a = newmatrix.els[0][3];
	res.b = newmatrix.els[1][3];
	res.c = newmatrix.els[2][3];
	return res;
}

struct colourt{
	double t;
	double red;
	double green;
	double blue;
	colourt operator+(const colourt b)
	{
		colourt c = {t + b.t, red + b.red, green+b.green, blue+b.blue};
		return c;
	}
};
const colourt background = {-1, 0.125, 0.125, 0.25};

colourt operator*(colourt a, double b){
	colourt res = {0, a.red*b, a.green*b, a.blue*b};
	return res;
}

class Light{
private:
	point3 pos;
	colourt col;
	bool distant;
	double intensity;
public:
	Light(point3 apos, colourt acol, bool adist, double aint){
		pos = apos; col = acol; distant = adist; intensity = aint;
	}
	point3 get_pos(){return pos;}
	colourt get_col(){return col;}
	bool isdistantlight(){return distant;}
	double get_intensity(){return intensity;}
};

class diffuseLight : private Light{
private:
	double radius;
	int numpoints;
public:

	diffuseLight(point3 apos, colourt acol, double aint, double arad, int anumpoints) : Light(apos, acol, false, aint){
		radius = arad;
		numpoints = anumpoints;
	}

	void getuniques(std::vector<Light>& lights){

		int random = rand() % numpoints;
		

		double offset = 2.0/numpoints;
		double increment = M_PI * (3.0 - sqrt(5));

		for(int i = 0; i < numpoints; i++){
			double y = ((i*offset - 1) + offset/2);
			double r = sqrt(1 - pow(y, 2));

			double phi = ((i+random) % numpoints) * increment;
			double x = cos(phi) * r;
			double z = sin(phi) * r;

			point3 newpos = {x, y, z};
			newpos = newpos * radius + get_pos();

			lights.push_back(Light(newpos, get_col(), isdistantlight(), get_intensity()/numpoints));
		}
	}

};

struct shapeproperties{
	colourt col;
	double ka;
	double ks;
	double kd;
	double alpha;
	double albedo;
	double transparency;
	double density;
};

struct lightcolours{
	point3 is;
	point3 ia;
	point3 id;
};

const colourt ERR = {-1, 0, 0, 0};

class Shape{
public:
	virtual shapeproperties get_refconsts(point3) = 0;
	virtual colourt rayquery(ray, std::vector<Light>, std::vector<Shape*>, int) = 0;
	virtual bool cast_ray(ray) = 0;
	virtual bool obviousintersect(point3) = 0;
};

colourt Phong(shapeproperties ref_consts, point3 rayorigin, point3 target, lightcolours l_cols, point3 N, std::vector<Shape*> shapes, std::vector<Light> lights){
	
	colourt res;
	point3 illumination; 
	illumination = l_cols.ia*ref_consts.ka;

	res.red = illumination.x;
	res.green = illumination.y;
	res.blue = illumination.z;
	
	//Separately for each light
	for(int l = 0; l < lights.size(); l++){
		point3 L;
		Light light = lights[l];
 		
 		//Light intensity
 		//-----------------------------------------------------
 		double intensity;
 		if(light.isdistantlight()){
			L = light.get_pos()/veclen(light.get_pos());
			intensity = light.get_intensity();
		}
		else{
			L = (light.get_pos() - target)/veclen(light.get_pos() - target);
			intensity = light.get_intensity()/pow(pointdist(target, light.get_pos()), 2);
		}
		//------------------------------------------------------

		//Own colour
		//------------------------------------------------------
		point3 is = {light.get_col().red, light.get_col().green, light.get_col().blue};
		point3 id = (l_cols.id* is)*(ref_consts.albedo/3.14*light.get_intensity());
		//------------------------------------------------------

		

		//Geometry
		//------------------------------------------------------
		point3 R = N * 2 * scalarprod(L, N) - L;

		point3 V = (rayorigin - target)/veclen(rayorigin - target);

		point3 shadowstart = target + N * epsilon;
		ray shadowray = linetwopts(shadowstart, light.get_pos());
		//------------------------------------------------------


		//Visibility and shadows
		//------------------------------------------------------
		bool vis = true;
		for(int i = 0; i<shapes.size(); i++){
			vis = vis && !(shapes[i]->cast_ray(shadowray));
		}
		//------------------------------------------------------

		if(vis){
			illumination.x = 0; illumination.y = 0; illumination.z = 0;

			//Diffuse light
			if(scalarprod(L, N) > 0){
				illumination = illumination + id * scalarprod(L, N) * ref_consts.kd ;
			}

			//Specular light
			if(scalarprod(R,V) > 0)
				illumination = illumination + is * pow(scalarprod(R,V), ref_consts.alpha) * ref_consts.ks * light.get_intensity() * ref_consts.albedo / 3.14;

			res.red += illumination.x;
			res.green += illumination.y;
			res.blue += illumination.z;	
		}
	}
	return res;
}

colourt askallshapes(ray x, std::vector<Light> lights, std::vector<Shape*> shapes, int times = 0){
	colourt result, temp;
	result = background;
	for(int m = 0; m < shapes.size(); m++){
		temp = shapes[m]->rayquery(x, lights, shapes, times);
		if(temp.t>=0 && (result.t > temp.t || result.t <= 0)){
			result = temp;
		}
	}

	return result;
}

point3 refract(point3 I, point3 N, double density){
	double cosi;
	if(scalarprod(I, N) < -1)
		cosi = -1;
	if(scalarprod(I, N)>1)
		cosi = 1;
	else cosi = scalarprod(I,N);
	double etai = 1, etat = density;
	point3 n = N;
	if(cosi < 0) {cosi = -cosi;}
	else {std::swap(etai, etat); n = N*(-1);}
	double eta = etai/etat;
	double k = 1 - eta*eta*(1-cosi*cosi);
	if(k < 0){
		point3 null = {0, 0, 0};
		return null;
	}
	else
		return I*eta + n*(eta*cosi - sqrt(k));
}

double fresnel(point3 I, point3 N, double density){
	double cosi, kr;
	if(scalarprod(I, N) < -1)
		cosi = -1;
	if(scalarprod(I, N)>1)
		cosi = 1;
	else cosi = scalarprod(I,N);
    double etai = 1, etat = density; 
    if (cosi > 0) { std::swap(etai, etat); } 
    // Compute sini using Snell's law
    double sint = etai / etat * sqrtf(std::max(0.0, 1 - cosi * cosi)); 
    // Total internal reflection
    if (sint >= 1) { 
        kr = 1; 
    } 
    else { 
        float cost = sqrtf(std::max(0.0, 1 - sint * sint)); 
        cosi = fabsf(cosi); 
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
        kr = (Rs * Rs + Rp * Rp) / 2; 
    } 
    return kr;
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}
class Plane : public Shape{
private:
	point3 normal;
	double distance;
	// ax + by + cz = d
	shapeproperties ref_consts;
public:
	Plane(double na, double nb, double nc, double d, colourt x){
		normal.x = na; normal.y = nb; normal.z = nc;
		distance = d;
		ref_consts.col = x;
		ref_consts.ks = 0.2;
		ref_consts.kd = 0.8;
		ref_consts.ka = 0;
		ref_consts.alpha = 25;
		ref_consts.albedo = 0.15;
		ref_consts.transparency = 0;
		ref_consts.density = 1;
	}
	shapeproperties get_refconsts(point3 point){
		shapeproperties result = ref_consts;
		colourt black = {0, 0.05, 0.05, 0.05};
		if(((int) abs(floor(point.x))) % 2 == ((int) abs(floor(point.y))) % 2){
			result.col = black;
		}
		else{
			result.col.red = std::max(0.75*(100 - point.y)/100, 0.0);
			result.col.green = 0.75 - result.col.red;
		}

		return result;
	}
	bool cast_ray(ray x){
		if(scalarprod(normal, x.direction) == 0){
			if(scalarprod(normal, x.origin) == distance){
				return true;
			}
			else
				return false;
		}
		else{
			int t = (distance - scalarprod(normal, x.origin)) / scalarprod(normal, x.direction);
			return (t > 0);
		}
	}


	colourt rayquery(ray x, std::vector<Light> lights, std::vector <Shape*> shapes, int times = 0){

		colourt res;
		colourt col = ref_consts.col;
		if(scalarprod(normal, x.direction) == 0){
			if(scalarprod(normal, x.origin) == distance){
				res = col; res.t = 0;
			}
			else
				return ERR;
		}
		else{
			res = col;
			res.t = (distance - scalarprod(normal, x.origin))/scalarprod(normal, x.direction);
		}
				point3 rayorigin = x.origin;

		double t = res.t;

		point3 target = rayorigin + x.direction * t;

		if(pointdist(rayorigin, target) >= 500)
			return ERR;

		shapeproperties props = get_refconsts(target);
		col = props.col;

		lightcolours l_cols = {{1,1,1}, {0, 0, 0}, {col.red, col.green, col.blue}};
		
		
		point3 N = normal;
		point3 V = (rayorigin - target)/veclen(rayorigin - target);

		if(scalarprod(N, V) < 0){
			N.x *= -1; N.y *= -1; N.z *= -1;
		}
		
		res = Phong(props, rayorigin, target, l_cols, N, shapes, lights);
		res.t = t;

		if(pointdist(rayorigin, target) >= 500 || cast_ray(x) == false)
			return ERR;

		if(ref_consts.transparency == 0 || times >=3)
			return res;
		else{
			colourt oldres = res;
			point3 newtarget = target + N*epsilon;

			point3 R = N * 2 * scalarprod(x.direction*(-1), N) + x.direction;

			ray newray;
			newray.origin = newtarget;
			newray.direction = R;


			colourt newres = askallshapes(newray, lights, shapes, times+1);

			res = (oldres*(1-ref_consts.transparency)) + (newres*ref_consts.transparency);
			res.t = t;
			return res;
		}
	}

	bool obviousintersect(point3 p){
		return true;
	}
};

class Sphere : public Shape{
private:
	point3 centre;
	double radius;
	shapeproperties ref_consts;
public:
	Sphere(point3 a, double rad, colourt x, double tr)
	{
		centre = a; radius = rad;
		ref_consts.col = x;
		ref_consts.ks = 0.3; //specular coefficient
		ref_consts.kd = 0.7; //diffuse coefficient
		ref_consts.ka = 0.1; //ambient coefficient
		ref_consts.alpha = 200;
		ref_consts.albedo = 0.18;
		ref_consts.transparency = tr;
		ref_consts.density = 1.35;
	}

	shapeproperties get_refconsts(point3 point){
		return ref_consts;
	}
	bool cast_ray(ray x)
	{
		point3 rayvec = x.direction;
		double lenvec = veclen(rayvec);
		
		point3 rayorigin = x.origin;
		point3 ray1 = x.origin + x.direction;

		double distlinecent = crossprodmod(rayorigin - centre, ray1 - centre);
		distlinecent /= pointdist(rayorigin, ray1);

		if(pointdist(rayorigin, centre) <= radius || pointdist(ray1, centre) <= radius)
			return true;

		if(distlinecent > radius)
			return false;

		double dist1 = sqrt(radius*radius - distlinecent*distlinecent);
		double dist2 = pointdist(centre, rayorigin);
		double dist3 = sqrt(dist2*dist2 - distlinecent*distlinecent);
		double dist4 = pointdist(rayorigin, ray1);
		double t = - dist3 / dist4;
		if(t == 0) return true;
		point3 nearestpoint = x.origin + x.direction * t;

		if(pointdist(nearestpoint, centre) > radius)
		{
			nearestpoint = x.origin - x.direction * t;
		}


		point3 vec1 = rayorigin - nearestpoint;
		point3 vec2 = ray1 - nearestpoint;
		double a;
		if(vec1.x!=0)
			a = vec2.x / vec1.x;
		else if(vec1.y!=0) a = vec2.y / vec1.y;
		else if(vec1.z!=0) a = vec2.z / vec1.z;
		else return true;
		return a <= 1;

	}

	colourt rayquery(ray x, std::vector<Light> lights, std::vector<Shape*> shapes, int times = 0)
	{
		point3 rayorigin = x.origin;
		point3 ray1 = x.origin + x.direction;

		double distlinecent = crossprodmod(centre - rayorigin, centre - ray1);
		distlinecent /= pointdist(rayorigin, ray1);
		if(distlinecent > radius)
			return ERR;

		double dist1 = sqrt(radius*radius - distlinecent*distlinecent);
		double dist2 = pointdist(centre, rayorigin);
		double dist3 = sqrt(dist2*dist2 - distlinecent*distlinecent);
		double dist4 = pointdist(rayorigin, ray1);
		double t = (dist3 - dist1)/dist4;
		point3 target = x.origin + x.direction * t;

		shapeproperties props = get_refconsts(target);
		colourt col = props.col;

		lightcolours l_cols = {{1,1,1}, {0, 0, 0}, {col.red, col.green, col.blue}};
			
		point3 N = (target - centre)/veclen(target - centre);
		colourt res;
					
		res = Phong(props, rayorigin, target, l_cols, N, shapes, lights);
		res.t = t;

		if(pointdist(rayorigin, target) >= 500 || cast_ray(x) == false)
			return ERR;

		if(ref_consts.transparency == 0 || times >=10)
			return res;
		else{
			colourt oldres = res;
			point3 newtarget = target + N*epsilon;

			colourt reflectedcol = {0, 0, 0, 0};
			colourt refractedcol = {0, 0, 0, 0};
			double kr = fresnel(x.direction, N, props.density);

			if(kr > 0){
				point3 R = N * 2 * scalarprod(x.direction*(-1), N) + x.direction;

				ray reflectedray;
				reflectedray.origin = newtarget;
				reflectedray.direction = R;
				reflectedcol = askallshapes(reflectedray, lights, shapes, times+1);
			}
			if(1 - kr> 0){
				ray refractedray;
				refractedray.origin = newtarget;
				refractedray.direction = refract(x.direction, N, props.density);

				refractedray.origin = intersect(refractedray);
				point3 N_new = (refractedray.origin - centre)/veclen(refractedray.origin - centre);

				refractedray.origin = refractedray.origin + N_new*epsilon;
				refractedray.direction = refract(refractedray.direction, N_new, props.density);

				//std::cout<<refractedray.origin.x<<" "<<target.x<<std::endl;
				std::vector<Shape*> v;
				for(int i = 0; i<shapes.size(); i++){
					if(shapes[i]!=this) v.push_back(shapes[i]);
				}
				refractedcol = askallshapes(refractedray, lights, v, times+1);
			}

			res = oldres*(1 - props.transparency) + (reflectedcol*kr + refractedcol*(1-kr))*props.transparency;
			res.t = t;
			return res;
		}

	}

	point3 intersect(ray x){
		point3 rayorigin = x.origin;
		point3 ray1 = x.origin + x.direction;

		double distlinecent = crossprodmod(centre - rayorigin, centre - ray1);
		distlinecent /= pointdist(rayorigin, ray1);

		double dist1 = sqrt(radius*radius - distlinecent*distlinecent);
		double dist2 = pointdist(centre, rayorigin);
		double dist3 = sqrt(dist2*dist2 - distlinecent*distlinecent);
		double dist4 = pointdist(rayorigin, ray1);
		double t = (dist3 - dist1)/dist4;
		point3 target = x.origin + x.direction * t;
		return target;
	}

	bool obviousintersect(point3 pixel){
		point3 right = {1, 0, 0};
		point3 down = {0, 0, 1};
		point3 forward = {0, 1, 0};
		point3 proj = {centre.x, pixel.y, centre.z};
		return pointdist(proj, pixel) <= radius;
	}
};



#endif