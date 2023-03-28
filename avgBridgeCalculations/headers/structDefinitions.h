#ifndef STRUCTDEFINITIONSHEUR_H
#define STRUCTDEFINITIONSHEUR_H

typedef struct boundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz;
	float xLength, yLength, zLength;
} BOUNDARY;

typedef struct trajectory
{
	int sino, atomType, ix, iy, iz;
	float x, y, z;
	int adsorbedID;
} TRAJECTORY;

typedef struct bondInfo
{
	float x1, y1, z1, x2, y2, z2;
	float xc, yc, zc;
	float xOrientationAngle;
	int index1, index2;
} BONDINFO;

typedef struct yDistribution
{
	float ylo, yhi;
	int count;
} YDIST;

typedef struct brigdesBin
{
	int count;
	float y1lo, y1hi, y2lo, y2hi;
} BRIDGESBIN;

typedef struct states
{
	float nBridges, nLoops, nDangles, nFreeChains;
} STATES;

typedef struct angleDistribution
{
	float angleLo, angleHi;
	int count;
} ANGLE_DISTRIBUTION;

#endif