// map2bmp
//=============================================================================
// Author: Michael Rieder (2013)
//
// Opens files produced by amr2map and produces BMP images
// Compile with -lm switch (link with math library)
//=============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

// structures
typedef struct mapdata_s {
	float *data;
	int nx;
	int ny;
} mapdata_t;

// BMP file header needs 2-byte alignment
#pragma pack(2)
typedef struct bmp_header_s {
	uint16_t bfType;
	uint32_t bfSize;
	uint16_t bfReserved1;
	uint16_t bfReserved2;
	uint32_t bfOffBits;
} bmp_header_t;

typedef struct bmp_infoheader_s {
  uint32_t biSize;
  int32_t  biWidth;
  int32_t  biHeight;
  uint16_t biPlanes;
  uint16_t biBitCount;
  uint32_t biCompression;
  uint32_t biSizeImage;
  int32_t  biXPelsPerMeter;
  int32_t  biYPelsPerMeter;
  uint32_t biClrUsed;
  uint32_t biClrImportant;
} bmp_infoheader_t;
#define BI_RGB 0
#pragma pack()

// global variables
const char *input_filename=0;
const char *output_filename=0;
char				output_fnstring[64];
int logarithmic_scale=0;
int do_min=0;
int do_max=0;
float minimum_value, maximum_value;
float gaussian_sigma=0.0;

// function declarations
void gaussian_blur( mapdata_t *mapdata, float sigma );
int open_mapfile( const char *mapfilename, mapdata_t *mapdata );
int save_bmp( mapdata_t *mapdata, const char *output_filename );
int parse_args( int argc, char *argv[] );
uint8_t heatmap_red( float heat );
uint8_t heatmap_green( float heat );
uint8_t heatmap_blue( float heat );
float	min( float a, float b );
float	max( float a, float b );


// function definitions

void gaussian_blur( mapdata_t *mapdata, float sigma ) {
  int imagesize;
  int i;
  float factor;

  imagesize = mapdata->nx * mapdata->ny; 
  factor = 0.5 / (sigma * sigma);
 
  for (i=0; i<imagesize;i++ ) {

    mapdata->data[i] = 0;
  }
}

int open_mapfile( const char *mapfilename, mapdata_t *mapdata )
{
	FILE	*fp;
	uint32_t 	header[4];
	int 	width, height, size;

	fp = fopen( mapfilename, "rb");

	if (fp){
		fread( header, sizeof(header), 1, fp );
		width = header[1];
		height = header[2];
		fread( &size, sizeof(size), 1, fp );
		if (size == width*height*sizeof(float)) {
			mapdata->data = malloc( size );
			mapdata->nx = width;
			mapdata->ny = height;
			fread( mapdata->data, size, 1, fp );
		}
		else {
			printf( "size mismatch in %s\n", mapfilename );
			return EXIT_FAILURE;
		}
	}
	else {
		printf( "error opening %s\n", mapfilename );
		return EXIT_FAILURE;
	}

  // apply gaussian filter
  if ( gaussian_sigma > 0.0 )
    gaussian_blur( mapdata, gaussian_sigma );

	fclose( fp );
	return EXIT_SUCCESS;
}

int save_bmp( mapdata_t *mapdata, const char *output_filename )
{
	bmp_header_t bmp_header;
	bmp_infoheader_t bmp_infoheader;
	uint8_t *bmpdata;
	int imagesize;
	int bmp_size;
	int bmp_linewidth;
	int i, k;
	float mapmin, mapmax, range;
	float heat;
	FILE *fp;

	imagesize = mapdata->nx * mapdata->ny;
	// each line has 4 byte padding
	bmp_linewidth = mapdata->nx*24/8;
	bmp_linewidth = bmp_linewidth + (4-bmp_linewidth%4)%4;
	
	bmp_size = bmp_linewidth * mapdata->ny;
	
	bmpdata = malloc( bmp_size );
	
	// produce image data
	// logarithmic scale
	if (logarithmic_scale == 1) {
		printf( "log " );
		for (i=0; i<imagesize;i++ )
			mapdata->data[i] = log10( mapdata->data[i] );
	}

	// min, max
	if (do_min == 1)
		mapmin = minimum_value;
	else {
		mapmin = mapdata->data[0];
		for (i=1;i<imagesize;i++)
			mapmin = min( mapmin, mapdata->data[i] );
	}
	if (do_max == 1)
		mapmax = maximum_value;
	else {
		mapmax = mapdata->data[0];
		for (i=1;i<imagesize;i++)
			mapmax = max( mapmax, mapdata->data[i] );
	}
	range = mapmax-mapmin;
	printf( "Min: %f Max: %f\n", mapmin, mapmax );
	
	// colors
	for (i=0;i<mapdata->ny;i++) {
		for (k=0;k<mapdata->nx;k++) {
			heat = mapdata->data[mapdata->nx*i+k] - mapmin;
			heat = heat / range;
			bmpdata[bmp_linewidth*i +3*k] = heatmap_blue( heat );
			bmpdata[bmp_linewidth*i +3*k+1] = heatmap_green( heat );
			bmpdata[bmp_linewidth*i +3*k+2] = heatmap_red( heat );
		}
	}

	// set file header
	bmp_header.bfType = 19778;
	bmp_header.bfSize = sizeof(bmp_infoheader_t)+sizeof(bmp_header_t)+bmp_size;
	bmp_header.bfReserved1 = 0;
	bmp_header.bfReserved2 = 0;
	bmp_header.bfOffBits = sizeof(bmp_infoheader_t)+sizeof(bmp_header_t);

	// set bmp info header
	bmp_infoheader.biSize = sizeof(bmp_infoheader_t);
	bmp_infoheader.biWidth = mapdata->nx;
	bmp_infoheader.biHeight = mapdata->ny;
	bmp_infoheader.biPlanes = 1;
	bmp_infoheader.biBitCount = 24;
	bmp_infoheader.biCompression = BI_RGB;
	bmp_infoheader.biSizeImage = bmp_size;
	bmp_infoheader.biXPelsPerMeter = 0;
	bmp_infoheader.biYPelsPerMeter = 0;
	bmp_infoheader.biClrUsed = 0;
	bmp_infoheader.biClrImportant = 0;

	// write to file
	fp = fopen( output_filename, "wb");

	if (fp){
		fwrite( &bmp_header, sizeof(bmp_header_t), 1, fp );
		fwrite( &bmp_infoheader, sizeof(bmp_infoheader_t), 1, fp );
		fwrite( bmpdata, bmp_size, 1, fp );
	}
	else {
		printf( "error creating %s\n", output_filename );
		return EXIT_FAILURE;
	}

  // done
	fclose( fp );
	free(bmpdata);
	return EXIT_SUCCESS;
}

int main (int argc, char *argv[])
{
	mapdata_t	mapdata;
	int ret;

	// check for float size - just in case...
	if (sizeof(float) != 4) {
		printf( "error: system float type %lu is incompatible to map file float type (4)!!\n", sizeof(float));
		return EXIT_FAILURE;
	}

	// get command line arguments
	if ( parse_args(argc,argv) == EXIT_FAILURE ) {
		printf( "Usage: %s <map.dat>\n", argv[0] );
		printf( "options:\n" );
		printf( "\t-log <0/1>\t logarithmic scale (default 0)\n" );
		printf( "\t-min <minimum>\t minimum value (default auto)\n" );
		printf( "\t-max <maximum>\t maximum value (default auto)\n" );
		printf( "\t-out <filename>\t output file (default input.bmp)\n" );
		return EXIT_FAILURE;
	}

	// open file
	if ( open_mapfile( input_filename, &mapdata ) == EXIT_FAILURE )
		return EXIT_FAILURE;

	printf( "%s, %ix%i", input_filename, mapdata.nx, mapdata.ny );
	// save to bmp
	ret = save_bmp( &mapdata, output_filename );

	// free the allocated map memory
	free( mapdata.data );

	if ( ret == EXIT_FAILURE ) {
		return EXIT_FAILURE;
	}	

	return EXIT_SUCCESS;
}

int parse_args( int argc, char *argv[] )
{
	int i;


	for (i=1;i<argc;i++) {
		if ( !strcmp(argv[i],"-log") ) {
			i++;
			sscanf( argv[i], "%i", &logarithmic_scale);
		}
		else if ( !strcmp(argv[i],"-min") ) {
			i++;
			sscanf( argv[i], "%f", &minimum_value);
			do_min = 1;
		}
		else if ( !strcmp(argv[i],"-max") ) {
			i++;
			sscanf( argv[i], "%f", &maximum_value);
			do_max = 1;
		}
		else if ( !strcmp(argv[i],"-out") ) {
			i++;
			output_filename = argv[i];
		}
		else if ( !strcmp(argv[i],"-gau") ) {
			i++;
			sscanf( argv[i], "%f", &gaussian_sigma);
		}
		else if ( argv[i][0] == '-' ) {
			printf( "unrecognized command: %s\n", argv[i] );
			return EXIT_FAILURE;
		}
		else input_filename = argv[i];
	}

	if (!input_filename)
		return EXIT_FAILURE;
	if (!output_filename) {
		strcpy( output_fnstring, input_filename );
		strcat( output_fnstring, ".bmp" );
		output_filename = output_fnstring;
	}
	return EXIT_SUCCESS;
}

float	min( float a, float b ) {
	if ( a < b )
		return a;
	else
		return b;
}

float	max( float a, float b ) {
	if ( a > b )
		return a;
	else
		return b;
}

uint8_t heatmap_red( float heat ) {

	if ( heat < 0.5 )
		return 0;
	else if ( heat < 0.75 )
		return 255 * (4*heat-2);
	else
		return 255;
}

uint8_t heatmap_green( float heat ) {

	if ( heat < 0 )
		return 0;
	else if ( heat < 0.25 )
		return 255 * 4*heat;
	else if ( heat < 0.75 )
		return 255;
	else if ( heat > 1 )
		return 0;
	else
		return 255 * (4 - 4*heat);
}

uint8_t heatmap_blue( float heat ) {

	if ( heat < 0.25 )
		return 255;
	else if ( heat < 0.5 )
		return 255 * (2 - 4*heat);
	else
		return 0;
}

