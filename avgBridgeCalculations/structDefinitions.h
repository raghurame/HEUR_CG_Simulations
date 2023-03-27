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

typedef struct bondIndo
{
	float x1, y1, z1, x2, y2, z2;
	float xc, yc, zc;
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