#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

/*
This program computes the number of states, such as

(i) bridges
(ii) loops
(iii) dangles
(iv) free

INPUT ARGUMENTS:
~~~~~~~~~~~~~~~

{~} argv[0] = program
{~} argv[1] = input cluster data
{~} argv[2] = input data file

*/

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

typedef struct datafile_atoms
{
	int resNumber;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
	int ix, iy, iz;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2;
	int isBridge, isDangle, isFree, isLoop;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

int countNAtoms (int nAtoms, const char filename[])
{
	FILE *inputCluster;
	char *pipeString, lineString[5000];

	pipeString = (char *) malloc (1000 * sizeof (char));

	if (strstr (filename, ".xz")) {
		snprintf (pipeString, 1000, "xzcat %s", filename);
		inputCluster = popen (pipeString, "r"); }
	else {
		inputCluster = fopen (filename, "r"); }


	for (int i = 0; i < 4; ++i)
	{
		fgets (lineString, 1000, inputCluster);
	}

	sscanf (lineString, "%d\n", &nAtoms);
	printf("detected %d atoms from the input cluster file...\n", nAtoms);

	fclose (inputCluster);

	return nAtoms;
}

int countNTimesteps (int nTimesteps_cluster, const char filename[])
{
	FILE *inputCluster;
	char *pipeString, lineString[5000];

	pipeString = (char *) malloc (1000 * sizeof (char));

	if (strstr (filename, ".xz")) {
		snprintf (pipeString, 1000, "xzcat %s | grep \"ITEM: TIMESTEP\" | wc -l", filename);
		inputCluster = popen (pipeString, "r"); }
	else {
		snprintf (pipeString, 1000, "cat %s | grep \"ITEM: TIMESTEP\" | wc -l", filename);
		inputCluster = popen (pipeString, "r"); }

	fgets (lineString, 1000, inputCluster);
	sscanf (lineString, "%d\n", &nTimesteps_cluster);

	printf("Number of timesteps in input cluster file: %d...\n", nTimesteps_cluster);

	return nTimesteps_cluster;
}

DATAFILE_INFO readData (const char *dataFileName, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers)
{
	printf("Reading LAMMPS data file...\n");
	FILE *input;
	input = fopen (dataFileName, "r");

	int isAtomLine = 0, /*nAtoms = -1,*/ nAtomLine = 0;
	int isBondLine = 0, /*nBonds = -1,*/ nBondLine = 0;
	int isAngleLine = 0, /*nAngles = -1,*/ nAngleLine = 0;
	int isDihedralLine = 0, /*nDihedrals = -1,*/ nDihedralLine = 0;
	int isImproperLine = 0, /*nImpropers = -1,*/ nImproperLine = 0;
	int printHeaderInfo = 1;

	DATAFILE_INFO datafile;
	datafile.nAtoms = -1;
	datafile.nBonds = -1;
	datafile.nAngles = -1;
	datafile.nDihedrals = -1;
	datafile.nImpropers = -1;

	char lineString[1000];

	// DATA_ATOMS *atoms;
	// DATA_BONDS *bonds;
	// DATA_ANGLES *angles;
	// DATA_DIHEDRALS *dihedrals;
	// DATA_IMPROPERS *impropers;
	*atoms = NULL;
	*bonds = NULL;
	*angles = NULL;
	*dihedrals = NULL;
	*impropers = NULL;

	while ((fgets (lineString, 1000, input) != NULL))
	{
		if (strstr (lineString, "atoms"))
		{
			sscanf (lineString, "%d \n", &datafile.nAtoms);
			(*atoms) = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));
		}

		if (strstr (lineString, "bonds"))
		{
			sscanf (lineString, "%d \n", &datafile.nBonds);
			(*bonds) = (DATA_BONDS *) malloc (datafile.nBonds * sizeof (DATA_BONDS));
		}

		if (strstr (lineString, "angles"))
		{
			sscanf (lineString, "%d \n", &datafile.nAngles);
			(*angles) = (DATA_ANGLES *) malloc (datafile.nAngles * sizeof (DATA_ANGLES));
		}

		if (strstr (lineString, "dihedrals"))
		{
			sscanf (lineString, "%d \n", &datafile.nDihedrals);
			(*dihedrals) = (DATA_DIHEDRALS *) malloc (datafile.nDihedrals * sizeof (DATA_DIHEDRALS));
		}

		if (strstr (lineString, "impropers"))
		{
			sscanf (lineString, "%d \n", &datafile.nImpropers);
			(*impropers) = (DATA_IMPROPERS *) malloc (datafile.nImpropers * sizeof (DATA_IMPROPERS));
		}

		if (strstr (lineString, "atom types"))
			sscanf (lineString, "%d \n", &datafile.nAtomTypes);

		if (strstr (lineString, "bond types"))
			sscanf (lineString, "%d \n", &datafile.nBondTypes);

		if (strstr (lineString, "angle types"))
			sscanf (lineString, "%d \n", &datafile.nAngleTypes);

		if (strstr (lineString, "dihedral types"))
			sscanf (lineString, "%d \n", &datafile.nDihedralTypes);

		if (strstr (lineString, "improper types"))
			sscanf (lineString, "%d \n", &datafile.nImproperTypes);

		if ((datafile.nAtoms >= 0) && (datafile.nBonds >= 0) && (datafile.nAngles >= 0) && (datafile.nDihedrals >= 0) && (datafile.nImpropers >= 0) && (printHeaderInfo))
			printHeaderInfo = 0;

		if (strstr (lineString, "Atoms"))
		{
			isAtomLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Bonds"))
		{
			isBondLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Angles"))
		{
			isAngleLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Dihedrals"))
		{
			isDihedralLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Impropers"))
		{
			isImproperLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (isAtomLine)
		{
			sscanf (lineString, "%d %d %d %f %f %f %d %d %d\n", 
				&(*atoms)[nAtomLine].id, 
				&(*atoms)[nAtomLine].molType, 
				&(*atoms)[nAtomLine].atomType, 
				&(*atoms)[nAtomLine].x, 
				&(*atoms)[nAtomLine].y, 
				&(*atoms)[nAtomLine].z,
				&(*atoms)[nAtomLine].ix,
				&(*atoms)[nAtomLine].iy,
				&(*atoms)[nAtomLine].iz);
			nAtomLine++;
			if (nAtomLine == datafile.nAtoms)
				isAtomLine = 0;
		}

		if (isBondLine)
		{
			sscanf (lineString, "%d %d %d %d\n", 
				&(*bonds)[nBondLine].id, 
				&(*bonds)[nBondLine].bondType, 
				&(*bonds)[nBondLine].atom1, 
				&(*bonds)[nBondLine].atom2);
			nBondLine++;
			if (nBondLine == datafile.nBonds)
				isBondLine = 0;
		}

		if (isAngleLine)
		{
			sscanf (lineString, "%d %d %d %d %d\n", 
				&(*angles)[nAngleLine].id, 
				&(*angles)[nAngleLine].angleType, 
				&(*angles)[nAngleLine].atom1, 
				&(*angles)[nAngleLine].atom2, 
				&(*angles)[nAngleLine].atom3);
			nAngleLine++;
			if (nAngleLine == datafile.nAngles)
				isAngleLine = 0;
		}

		if (isDihedralLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*dihedrals)[nDihedralLine].id, 
				&(*dihedrals)[nDihedralLine].dihedralType, 
				&(*dihedrals)[nDihedralLine].atom1, 
				&(*dihedrals)[nDihedralLine].atom2, 
				&(*dihedrals)[nDihedralLine].atom3, 
				&(*dihedrals)[nDihedralLine].atom4);
			nDihedralLine++;
			if (nDihedralLine == datafile.nDihedrals)
				isDihedralLine = 0;
		}

		if (isImproperLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*impropers)[nImproperLine].id, 
				&(*impropers)[nImproperLine].improperType, 
				&(*impropers)[nImproperLine].atom1, 
				&(*impropers)[nImproperLine].atom2, 
				&(*impropers)[nImproperLine].atom3, 
				&(*impropers)[nImproperLine].atom4);
			nImproperLine++;
			if (nImproperLine == datafile.nImpropers)
				isImproperLine = 0;
		}
	}

	printf("\nFrom input data file:\n\n nAtoms: %d\n nBonds: %d\n nAngles: %d\n nDihedrals: %d\n nImpropers: %d\n\n", datafile.nAtoms, datafile.nBonds, datafile.nAngles, datafile.nDihedrals, datafile.nImpropers);

	return datafile;
}

int main(int argc, char const *argv[])
{
	if (argc != 3)
	{
		printf("REQUIRED ARGUMENTS:\n\n {~} argv[0] = program\n {~} argv[1] = input cluster file\n {~} argv[2] = input data file\n\n");
		exit (1);
	}

	FILE *inputData, *inputCluster;
	char *pipeString, lineString[3000];
	pipeString = (char *) malloc (1000 * sizeof (char));

	inputData = fopen (argv[2], "r");

	if (strstr (argv[1], ".xz")) {
		snprintf (pipeString, 1000, "xzcat %s", argv[1]);
		inputCluster = popen (pipeString, "r"); }
	else {
		inputCluster = fopen (argv[1], "r"); }

	int nAtoms = countNAtoms (nAtoms, argv[1]), nTimesteps_cluster = countNTimesteps (nTimesteps_cluster, argv[1]);

	DATA_ATOMS *dataAtoms;
	DATA_BONDS *dataBonds;
	DATA_ANGLES *dataAngles;
	DATA_DIHEDRALS *dataDihedrals;
	DATA_IMPROPERS *dataImpropers;
	DATAFILE_INFO datafile;

	datafile = readData (argv[2], &dataAtoms, &dataBonds, &dataAngles, &dataDihedrals, &dataImpropers);

	int *cluster, *cluster_previous, *ghostParticleIDs, *ghostParticleClusterIDs, nGhostParticles = 0;
	cluster = (int *) malloc (nAtoms * sizeof (int));
	cluster_previous = (int *) malloc (nAtoms * sizeof (int));

	// Reading n-1 timesteps,
	// 'n' being the total number of timesteps.

	for (int i = 0; i < datafile.nBonds; ++i) {
		dataBonds[i].isBridge = -1; dataBonds[i].isDangle = -1; dataBonds[i].isFree = -1; dataBonds[i].isLoop = -1; }

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		if (dataAtoms[i].atomType == 2) {
			nGhostParticles++; }
	}

	ghostParticleIDs = (int *) malloc (nGhostParticles * sizeof (int));
	ghostParticleClusterIDs = (int *) malloc (nGhostParticles * sizeof (int));
	int currentGhostParticle = 0;

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		if (dataAtoms[i].atomType == 2) {
			ghostParticleIDs[currentGhostParticle] = dataAtoms[i].id;
			currentGhostParticle++; }
	}

	for (int i = 0; i < (nTimesteps_cluster - 1); ++i)
	{
		for (int i = 0; i < 9; ++i)
		{
			fgets (lineString, 2000, inputCluster);
		}

		for (int i = 0; i < nAtoms; ++i)
		{
			fgets (lineString, 2000, inputCluster);
			sscanf (lineString, "%d\n", &cluster[i]);
		}

		for (int i = 0; i < nGhostParticles; ++i)
		{
			ghostParticleClusterIDs[i] = cluster[ghostParticleIDs[i] - 1];
		}

		for (int i = 0; i < datafile.nBonds; ++i)
		{
			printf("%d (%d) %d (%d)\n", dataBonds[i].atom1, cluster[dataBonds[i].atom1 - 1], dataBonds[i].atom2, cluster[dataBonds[i].atom2 - 1]);

			// check if the cluster IDs belong to ghost particles
			// then classify the states
			// compare with the previous states to check transitions
			usleep (100000);
		}

		cluster_previous = cluster;
	}

	fclose (inputData);
	fclose (inputCluster);
	return 0;
}