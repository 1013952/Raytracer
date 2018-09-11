#ifndef _BMPHANDLING
#define _BMPHANDLING

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int prepbmp(std::ofstream & fout, int width, int height)
{
	// mimeType = "image/bmp";
	unsigned char file[14] = {
    	'B','M', // magic
    	0,0,0,0, // size in bytes
    	0,0, // app data
    	0,0, // app data
    	40+14,0,0,0 // start of data offset
	};
	unsigned char info[40] = {
    	40,0,0,0, // info hd size
    	0,0,0,0, // width
    	0,0,0,0, // heigth
    	1,0, // number color planes
    	24,0, // bits per pixel
    	0,0,0,0, // compression is none
    	0,0,0,0, // image bits size
    	0x13,0x0B,0,0, // horz resoluition in pixel / m
    	0x13,0x0B,0,0, // vert resolutions (0x03C3 = 96 dpi, 0x0B13 = 72 dpi)
    	0,0,0,0, // #colors in pallete
    	0,0,0,0, // #important colors
    	};
	int w=width;
	int h=height;

	int padSize  = (4-(w*3))%4;
	int sizeData = w*h*3 + h*padSize;
	int sizeAll  = sizeData + sizeof(file) + sizeof(info);

	file[2] = (unsigned char)( sizeAll    );
	file[3] = (unsigned char)( sizeAll>> 8);
	file[4] = (unsigned char)( sizeAll>>16);
	file[5] = (unsigned char)( sizeAll>>24);

	info[ 4] = (unsigned char)( w   );
	info[ 5] = (unsigned char)( w>> 8);
	info[ 6] = (unsigned char)( w>>16);
	info[ 7] = (unsigned char)( w>>24);

	info[ 8] = (unsigned char)( h    );
	info[ 9] = (unsigned char)( h>> 8);
	info[10] = (unsigned char)( h>>16);
	info[11] = (unsigned char)( h>>24);

	info[20] = (unsigned char)( sizeData    );
	info[21] = (unsigned char)( sizeData>> 8);
	info[22] = (unsigned char)( sizeData>>16);
	info[23] = (unsigned char)( sizeData>>24);

	fout.write( (char*)file, sizeof(file));
	fout.write( (char*)info, sizeof(info) );
	return padSize;
}

void bmppixel(std::ofstream & fout, double redpx, double greenpx, double bluepx){
	long red = lround( 255.0 * redpx);
        	if ( red < 0 ) red=0;
        	if ( red > 255 ) red=255;
       		long green = lround( 255.0 * greenpx);
        	if ( green < 0 ) green = 0;
        	if ( green > 255 ) green = 255;
        	long blue = lround( 255.0 * bluepx);
        	if ( blue < 0 ) blue = 0;
        	if ( blue > 255 ) blue = 255;

			unsigned char pixelc[3];
        	pixelc[0] = blue;
        	pixelc[1] = green;
        	pixelc[2] = red;

        	fout.write( (char*)pixelc, 3 );
}

#endif