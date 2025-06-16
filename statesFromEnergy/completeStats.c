#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

// These global variables are no longer needed
// I am taking these values from data file
// #define NPARTICLES 100
// #define NPOLYMERS 4000
// #define NBEADS 2
// #define COORDINATION_NUMBER 20

// This is important. Set this to high number
// Too big number requires high memory
// #define N_TIMEFRAMES_TO_CONSIDER 20000

typedef struct blocks
{
	int size;
	float average, covariance, stdev, stderr;
	float squaredAverage, averageSquare;
} BLOCKS;

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
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2;
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

typedef struct simulationBoundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
	float xLength, yLength, zLength;
} SIMULATION_BOUNDARY;

typedef struct bondStatus
{
	int id1, id2, id;
	bool isItLoop, isItBridge, isItDangle, isItFree;
	bool loop_corrected, bridge_corrected;
} BOND_STATUS;

typedef struct boundStatus
{
	bool isItBound;
	int bondNumber;
} BOUND_STATUS;

typedef struct dumpEnergy
{
	int atom1, atom2;
	float distance, energy;
} DUMP_ENERGY;

typedef struct states
{
	int nLoops, nBridges, nFree, nDangles;
} STATES;

void printStateNumber (bool isItFree, bool isItDangle, bool isItLoop, bool isItBridge)
{
	if (isItFree == 1)
	{
		fprintf(stdout, "0 ");
	}
	if (isItDangle == 1)
	{
		fprintf(stdout, "1 ");
	}
	if (isItLoop == 1)
	{
		fprintf(stdout, "2 ");
	}
	if (isItBridge == 1)
	{
		fprintf(stdout, "3 ");
	}

	fflush (stdout);
}

/*int *storeParticleIDs (int *particleIDs)
{
	for (int i = 0; i < NPARTICLES; ++i)
	{
		particleIDs[i] = (i + 1) + (i + 1) * COORDINATION_NUMBER;
	}

	return particleIDs;
}
*/

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
			sscanf (lineString, "%d %d %d %f %f %f\n", 
				&(*atoms)[nAtomLine].id, 
				&(*atoms)[nAtomLine].molType, 
				&(*atoms)[nAtomLine].atomType, 
				&(*atoms)[nAtomLine].x, 
				&(*atoms)[nAtomLine].y, 
				&(*atoms)[nAtomLine].z);
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

int countNTimeframes (const char *filename)
{
	int nTimeframes = 0;
	char *pipeString, lineString[500];
	pipeString = (char *) malloc (1000 * sizeof (char));

	if (strstr (filename, ".gz"))
	{
		snprintf (pipeString, 1000, "gzcat %s | grep \"ITEM: TIMESTEP\" | wc -l", filename);
	}
	else
	{
		snprintf (pipeString, 1000, "cat %s | grep \"ITEM: TIMESTEP\" | wc -l", filename);
	}

	FILE *lineCount;
	lineCount = popen (pipeString, "r");

	fgets (lineString, 500, lineCount);
	sscanf (lineString, "%d", &nTimeframes);
	printf("Number of lines: %s\n", lineString);
	printf("Number of timeframes: %d\n", nTimeframes);

	return nTimeframes;
}

DATA_BONDS *sortBonds (DATA_BONDS *sortedBonds, DATA_BONDS *dataBonds, int nBonds, int nAtoms)
{
	int currentIndex = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		for (int j = 0; j < nBonds; ++j)
		{
			if (dataBonds[j].atom1 == i + 1)
			{
				sortedBonds[currentIndex].id = dataBonds[j].id;
				sortedBonds[currentIndex].bondType = dataBonds[j].bondType;
				sortedBonds[currentIndex].atom1 = dataBonds[j].atom1;
				sortedBonds[currentIndex].atom2 = dataBonds[j].atom2;
				currentIndex++;

				break;
			}
		}
	}

	return sortedBonds;
}

int countDumpEntries (FILE *inputDump, int nDumpEntries)
{
	rewind (inputDump);

	char lineString[3000];

	for (int i = 0; i < 4; ++i)
	{
		fgets (lineString, 3000, inputDump);
	}

	sscanf (lineString, "%d", &nDumpEntries);

	return nDumpEntries;
}

DUMP_ENERGY *saveDumpEnergyEntries (FILE *inputDump, DUMP_ENERGY *energyEntries, int nDumpEntries, int energyColumn, int *status)
{
	char lineString[3000];
	int nDumpEntries_current = 0;
	int debugg = 0;
	(*status) = 0;

	for (int i = 0; i < 9; ++i)
	{
		fgets (lineString, 3000, inputDump);

		if (i == 3) {
			sscanf (lineString, "%d\n", &nDumpEntries_current); }
	}

	// fprintf(stdout, "nDumpEntries_current: %d\n", nDumpEntries_current);
	// fflush (stdout);

	for (int i = 0; i < nDumpEntries_current; ++i)
	{
		// fgets (lineString, 3000, inputDump);

		if (fgets (lineString, 3000, inputDump) != NULL)
		{
			(*status) = 1;
		}
		else
		{
			(*status) = 0;
			return energyEntries;
		}

		if (energyColumn == 4)
		{
			sscanf (lineString, "%d %d %*f %f\n", &energyEntries[i].atom1, &energyEntries[i].atom2, &energyEntries[i].energy);

			if (debugg == 1)
			{
				fprintf(stdout, "%d %d %f [%d]\n", energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].energy, (*status));
				usleep (100000);
			}
		}
		else if (energyColumn == 3)
		{
			sscanf (lineString, "%d %d %f\n", &energyEntries[i].atom1, &energyEntries[i].atom2, &energyEntries[i].energy);

			if (debugg == 1)
			{
				fprintf(stdout, "%d %d %f [%d]\n", energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].energy, (*status));
				usleep (100000);
			}
		}

/*		if (i < 5 || i > (nDumpEntries_current - 5))
		{
			printf("%d %d %f %f\n", energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].distance, energyEntries[i].energy);
			fflush (stdout);
		}
*/
	}

	// sleep (1);

	return energyEntries;
}

int *initBoundParticleIDs (int *boundParticleID, int nBonds)
{
	for (int i = 0; i < nBonds; ++i)
	{
		boundParticleID[i] = 0;
	}

	return boundParticleID;
}

BOND_STATUS **initPolymerBondStatus (BOND_STATUS **polymerBondStatus, int currentTimeframe, DATAFILE_INFO datafile)
{
	for (int i = 0; i < datafile.nBonds; ++i)
	{
		polymerBondStatus[i][currentTimeframe].isItFree = true;
		polymerBondStatus[i][currentTimeframe].isItDangle = false;
		polymerBondStatus[i][currentTimeframe].isItBridge = false;
		polymerBondStatus[i][currentTimeframe].isItLoop = false;
		polymerBondStatus[i][currentTimeframe].id1 = 0;
		polymerBondStatus[i][currentTimeframe].id2 = 0;
	}

	return polymerBondStatus;
}

BOUND_STATUS **initBeadBoundStatus (BOUND_STATUS **beadBoundStatus, int currentTimeframe, DATAFILE_INFO datafile)
{
	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		beadBoundStatus[i][currentTimeframe].isItBound = 0;
	}

	return beadBoundStatus;
}

BOUND_STATUS **findBoundStates (BOUND_STATUS **beadBoundStatus, DUMP_ENERGY *energyEntries, int nDumpEntries, DATAFILE_INFO datafile, int nTimeframes, int currentTimeframe, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0;
	beadBoundStatus = initBeadBoundStatus (beadBoundStatus, currentTimeframe, datafile);

	for (int i = 0; i < nDumpEntries; ++i)
	{
		if (debugg == 1) {
			printf("\n%d => %d %d %f => ", currentTimeframe, energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].energy); }

		if (currentTimeframe == 0)
		{
			if (energyEntries[i].energy < -1 && energyEntries[i].atom1 > 0 && energyEntries[i].atom2 > 0)
			{
				if (sortedAtoms[energyEntries[i].atom1 - 1].atomType == 1)
				{
					beadBoundStatus[energyEntries[i].atom1 - 1][currentTimeframe].isItBound = 1;

					if (debugg == 1) {
						printf("%d ", beadBoundStatus[energyEntries[i].atom1 - 1][currentTimeframe].isItBound); }
				}
				else if (sortedAtoms[energyEntries[i].atom2 - 1].atomType == 1)
				{
					beadBoundStatus[energyEntries[i].atom2 - 1][currentTimeframe].isItBound = 1;

					if (debugg == 1) {
						printf("%d ", beadBoundStatus[energyEntries[i].atom2 - 1][currentTimeframe].isItBound); }
				}
			}
		}
		else
		{
			if (energyEntries[i].energy < -1 && energyEntries[i].atom1 > 0 && energyEntries[i].atom2 > 0)
			{
				if (sortedAtoms[energyEntries[i].atom1 - 1].atomType == 1)
				{
					beadBoundStatus[energyEntries[i].atom1 - 1][currentTimeframe].isItBound = 1;

					if (debugg == 1) {
						printf("%d ", beadBoundStatus[energyEntries[i].atom1 - 1][currentTimeframe].isItBound); }
				}
				else if (sortedAtoms[energyEntries[i].atom2 - 1].atomType == 1)
				{
					beadBoundStatus[energyEntries[i].atom2 - 1][currentTimeframe].isItBound = 1;

					if (debugg == 1) {
						printf("%d ", beadBoundStatus[energyEntries[i].atom2 - 1][currentTimeframe].isItBound); }
				}
			}
			else if (energyEntries[i].energy < 0 && energyEntries[i].energy > -1 && energyEntries[i].atom1 > 0 && energyEntries[i].atom2 > 0)
			{
				if (sortedAtoms[energyEntries[i].atom1 - 1].atomType == 1)
				{
					if (beadBoundStatus[energyEntries[i].atom1 - 1][currentTimeframe - 1].isItBound == 1)
					{
						beadBoundStatus[energyEntries[i].atom1 - 1][currentTimeframe].isItBound = 1;

						if (debugg == 1) {
							printf("%d (prev: %f)", beadBoundStatus[energyEntries[i].atom1 - 1][currentTimeframe].isItBound, energyEntries[i - 1].energy); }
					}
				}
				else if (sortedAtoms[energyEntries[i].atom2 - 1].atomType == 1)
				{
					if (beadBoundStatus[energyEntries[i].atom2 - 1][currentTimeframe - 1].isItBound == 1)
					{
						beadBoundStatus[energyEntries[i].atom2 - 1][currentTimeframe].isItBound = 1;

						if (debugg == 1) {
							printf("%d (prev: %f)", beadBoundStatus[energyEntries[i].atom2 - 1][currentTimeframe].isItBound, energyEntries[i - 1].energy); }
					}
				}
			}
		}

		if (debugg == 1)
		{
			usleep (100000);
		}
	}

	return beadBoundStatus;
}

BOND_STATUS **checkBondStatus (BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, int nDumpEntries, DATAFILE_INFO datafile, int nTimeframes, int currentTimeframe, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0;
	int bondNumber = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;

	int *boundParticleID;
	boundParticleID = (int *) malloc (datafile.nBonds * sizeof (int));

	boundParticleID = initBoundParticleIDs (boundParticleID, datafile.nBonds);
	polymerBondStatus = initPolymerBondStatus (polymerBondStatus, currentTimeframe, datafile);

	for (int i = 0; i < nDumpEntries; ++i)
	{
		if (debugg == 1)
		{
			usleep (100000);
			fprintf(stdout, "%d [%d] %d [%d] %f ", energyEntries[i].atom1, sortedAtoms[energyEntries[i].atom1 - 1].atomType, energyEntries[i].atom2, sortedAtoms[energyEntries[i].atom2 - 1].atomType, energyEntries[i].energy);
		}

		if (energyEntries[i].energy < -1 && energyEntries[i].atom1 > 0 && energyEntries[i].atom2 > 0)
		{
			// check if one of the two atoms belongs to ghost particle and the another one belongs to polymer bead
			// if the above is true, then check if the energy is less than -1
			// if it is less than -1, then make that a dangle and store the particle ID corresponding to the bondNumber.
			// if the same particle ID is already present for the corresponding bondNumber, then make it a loop
			// if a different particle ID is present fot the bondNumber, then make it a bridge
			// if only one particle ID is associated with the bondNumber, then make it a dangle
			// if a bondNumber does not contain any particle ID, then it is a free chain

			if (sortedAtoms[energyEntries[i].atom1 - 1].atomType == 2 && sortedAtoms[energyEntries[i].atom2 - 1].atomType == 1)
			{
				// bondNumber = (energyEntries[i].atom2 - (energyEntries[i].atom2 / coordinationNumber) - 1) / NBEADS;
				bondNumber = floor ((energyEntries[i].atom2 - (int)floor (energyEntries[i].atom2 / (coordinationNumber + 1)) - 1) / 2);

				if (debugg == 1)
				{
					printf("(%d) ", bondNumber);
					printf("(%d) ", boundParticleID[bondNumber]);
					printStateNumber (polymerBondStatus[bondNumber][currentTimeframe].isItFree, polymerBondStatus[bondNumber][currentTimeframe].isItDangle, polymerBondStatus[bondNumber][currentTimeframe].isItLoop, polymerBondStatus[bondNumber][currentTimeframe].isItBridge);
				}

				if (boundParticleID[bondNumber] == 0)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItFree = false;
					boundParticleID[bondNumber] = energyEntries[i].atom1;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = true;
					polymerBondStatus[bondNumber][currentTimeframe].id1 = energyEntries[i].atom1;
					polymerBondStatus[bondNumber][currentTimeframe].id2 = 0;
				}
				else if (boundParticleID[bondNumber] == energyEntries[i].atom1)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItLoop = true;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = false;
					polymerBondStatus[bondNumber][currentTimeframe].id2 = energyEntries[i].atom1;
				}
				else if (boundParticleID[bondNumber] != energyEntries[i].atom1)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItBridge = true;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = false;
					polymerBondStatus[bondNumber][currentTimeframe].id2 = energyEntries[i].atom1;
				}

				if (debugg == 1)
				{
					printf("--> %d", boundParticleID[bondNumber]);
					printStateNumber (polymerBondStatus[bondNumber][currentTimeframe].isItFree, polymerBondStatus[bondNumber][currentTimeframe].isItDangle, polymerBondStatus[bondNumber][currentTimeframe].isItLoop, polymerBondStatus[bondNumber][currentTimeframe].isItBridge);
					printf("\n");
				}
			}
			else if (sortedAtoms[energyEntries[i].atom2 - 1].atomType == 2 && sortedAtoms[energyEntries[i].atom1 - 1].atomType == 1)
			{
				// bondNumber = (energyEntries[i].atom1 - (energyEntries[i].atom1 / coordinationNumber) - 1) / NBEADS;
				bondNumber = floor ((energyEntries[i].atom1 - (int)floor (energyEntries[i].atom1 / (coordinationNumber + 1)) - 1) / 2);

				if (debugg == 1)
				{
					printf("(%d) ", bondNumber);
					printf("(%d) ", boundParticleID[bondNumber]);
					printStateNumber (polymerBondStatus[bondNumber][currentTimeframe].isItFree, polymerBondStatus[bondNumber][currentTimeframe].isItDangle, polymerBondStatus[bondNumber][currentTimeframe].isItLoop, polymerBondStatus[bondNumber][currentTimeframe].isItBridge);
				}

				if (boundParticleID[bondNumber] == 0)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItFree = false;
					boundParticleID[bondNumber] = energyEntries[i].atom2;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = true;
					polymerBondStatus[bondNumber][currentTimeframe].id1 = energyEntries[i].atom2;
					polymerBondStatus[bondNumber][currentTimeframe].id2 = 0;
				}
				else if (boundParticleID[bondNumber] == energyEntries[i].atom2)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItLoop = true;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = false;
					polymerBondStatus[bondNumber][currentTimeframe].id2 = energyEntries[i].atom2;
				}
				else if (boundParticleID[bondNumber] != energyEntries[i].atom2)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItBridge = true;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = false;
					polymerBondStatus[bondNumber][currentTimeframe].id2 = energyEntries[i].atom2;
				}

				if (debugg == 1)
				{
					printf("--> %d", boundParticleID[bondNumber]);
					printStateNumber (polymerBondStatus[bondNumber][currentTimeframe].isItFree, polymerBondStatus[bondNumber][currentTimeframe].isItDangle, polymerBondStatus[bondNumber][currentTimeframe].isItLoop, polymerBondStatus[bondNumber][currentTimeframe].isItBridge);
					printf("\n");
				}
			}
		}
	}

	// Assigning 'free' for bonds which are not dangles/bridges/loops
	for (int i = 0; i < datafile.nBonds; ++i)
	{
		if (polymerBondStatus[i][currentTimeframe].isItFree + polymerBondStatus[i][currentTimeframe].isItDangle + polymerBondStatus[i][currentTimeframe].isItLoop + polymerBondStatus[i][currentTimeframe].isItBridge == 0)
		{
			polymerBondStatus[i][currentTimeframe].isItFree = 1;
		}
	}

/*	for (int i = 0; i < datafile.nBonds; ++i)
	{
		// printf("%d -> ", i);
		printStateNumber (polymerBondStatus[i][currentTimeframe].isItFree, polymerBondStatus[i][currentTimeframe].isItDangle, polymerBondStatus[i][currentTimeframe].isItLoop, polymerBondStatus[i][currentTimeframe].isItBridge);
	}
	printf("\n");
	sleep (100);
*/
	return polymerBondStatus;
}

BOND_STATUS **initBondStatus (BOND_STATUS **polymerBondStatus, int nBonds, int nTimeframes)
{
	for (int i = 0; i < nBonds; ++i)
	{
		for (int j = 0; j < nTimeframes; ++j)
		{
			polymerBondStatus[i][j].id = i;
			polymerBondStatus[i][j].isItBridge = false;
			polymerBondStatus[i][j].isItLoop = false;
			polymerBondStatus[i][j].isItDangle = false;
			polymerBondStatus[i][j].isItFree = false;
			polymerBondStatus[i][j].loop_corrected = false;
			polymerBondStatus[i][j].bridge_corrected = false;
		}
	}

	return polymerBondStatus;
}

DATA_ATOMS *sortAtoms (DATA_ATOMS *sortedAtoms, DATA_ATOMS *dataAtoms, int nAtoms)
{
	int currentIndex = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		for (int j = 0; j < nAtoms; ++j)
		{
			if (dataAtoms[j].id == (i + 1))
			{
				sortedAtoms[currentIndex].id = dataAtoms[j].id;
				sortedAtoms[currentIndex].atomType = dataAtoms[j].atomType;
				sortedAtoms[currentIndex].molType = dataAtoms[j].molType;

				currentIndex++;

				break;
			}
		}
	}
	return sortedAtoms;
}

DUMP_ENERGY *initEnergyEntries (DUMP_ENERGY *energyEntries, int nDumpEntries)
{
	for (int i = 0; i < nDumpEntries; ++i)
	{
		energyEntries[i].atom1 = 0;
		energyEntries[i].atom2 = 0;
		energyEntries[i].distance = 0;
		energyEntries[i].energy = 0;
	}

	return energyEntries;
}

STATES countStates (STATES polymerStates, BOND_STATUS **polymerBondStatus, DATAFILE_INFO datafile, int currentTimeframe)
{
	polymerStates.nLoops = 0;
	polymerStates.nBridges = 0;
	polymerStates.nDangles = 0;
	polymerStates.nBridges = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		if (polymerBondStatus[i][currentTimeframe].isItLoop) {
			polymerStates.nLoops++; }
		else if (polymerBondStatus[i][currentTimeframe].isItBridge) {
			polymerStates.nBridges++; }
		else if (polymerBondStatus[i][currentTimeframe].isItDangle) {
			polymerStates.nDangles++; }
	}

	return polymerStates;
}

void printStates (STATES polymerStates, FILE *outputStates, DATAFILE_INFO datafile)
{
	fprintf(outputStates, "%d %d %d %d\n", polymerStates.nBridges, polymerStates.nLoops, polymerStates.nDangles, datafile.nBonds - polymerStates.nBridges - polymerStates.nLoops - polymerStates.nDangles);
	fflush (outputStates);
}

float *initializeBlockValues (float *blocks, int nLines)
{
	for (int i = 0; i < nLines; ++i)
	{
		blocks[i] = -9999;
	}

	return blocks;
}

BLOCKS *computeBlockAverages (BLOCKS *blockAverages, int nLines, float *inputData)
{
	int currentEntry = 0, nEntries = 0;
	float *blocks;
	blocks = (float *) malloc (nLines * sizeof (blocks));

	for (int i = 0; i < nLines; ++i)
	{
		if (inputData[i] > -9999) {
			nEntries++; }
	}

	// 'i' is to iterate through various block sizes
	// 'i + 1' is the block size
	for (int i = 0; i < floor (nLines / 2); ++i)
	{
		blocks = initializeBlockValues (blocks, nLines);
		currentEntry = 0;

		// going through the blocks
		// summing the items in the blocks
		// number of blocks = floor (nLines / (i + 1))
		for (int j = 0; j < floor (nLines / (i + 1)); ++j)
		{
			// do averaging within the block
			// (i + 1) is the block size
			blocks[j] = 0;
			for (int k = 0; k < (i + 1); ++k)
			{
				blocks[j] += inputData[currentEntry];
				currentEntry++;
			}
		}

		// going through the blocks
		// finding local averages of these blocks
		// number of blocks = floor (nLines / (i + 1))
		for (int j = 0; j < floor (nLines / (i + 1)); ++j)
		{
			blocks[j] /= (i + 1);
		}

		blockAverages[i].average = 0;
		// find the global average across blocks
		for (int j = 0; j < floor (nLines / (i + 1)); ++j)
		{
			if (blocks[j] == -9999)
			{
				printf("ERROR: Blocks in default value !\n");
				exit (1);
			}

			blockAverages[i].average += blocks[j];
		}

		blockAverages[i].average /= (float)(floor (nLines / (i + 1)));

		// find the global stdev/covar across blocks
		for (int j = 0; j < floor (nLines / (i + 1)); ++j)
		{
			blockAverages[i].covariance += ((blocks[j] - blockAverages[i].average) * (blocks[j] - blockAverages[i].average));
		}

		blockAverages[i].covariance /= (float)(floor (nLines / (i + 1)) - 1);
		blockAverages[i].stdev = sqrt (blockAverages[i].covariance);
		blockAverages[i].stderr = blockAverages[i].stdev / sqrt ((float)((floor (nLines / (i + 1)))));
	}

	return blockAverages;
}

void saveInputData (float **inputData_nBridges, float **inputData_nLoops, float **inputData_nDangles, float **inputData_nFree, int nLines, FILE *file_data)
{
	char lineString[2000];
	rewind (file_data);

	for (int i = 0; i < nLines; ++i)
	{
		fgets (lineString, 2000, file_data);
		sscanf (lineString, "%f %f %f %f\n", &(*inputData_nBridges)[i], &(*inputData_nLoops)[i], &(*inputData_nDangles)[i], &(*inputData_nFree)[i]);
	}
}

BLOCKS *initializeBlocks (BLOCKS *blockAverages, int nLines)
{
	for (int i = 0; i < nLines; ++i) {
		blockAverages[i].size = i + 1; 
		blockAverages[i].average = 0;
		blockAverages[i].covariance = 0;
		blockAverages[i].stdev = 0;
		blockAverages[i].stderr = 0; }

	return blockAverages;
}

void printBlockAverageStats (FILE *file_output, BLOCKS *blockAverages, int nLines)
{
	for (int i = 0; i < floor (nLines / 2); ++i)
	{
		fprintf(file_output, "%d %f %f %f %f\n", i + 1, blockAverages[i].average, blockAverages[i].stdev, blockAverages[i].covariance*100000, blockAverages[i].stderr);
	}
}

BOND_STATUS checkingBackwards (BOND_STATUS backwardStatus, BOND_STATUS **polymerBondStatus, int currentBondIndex, int currentTimeframe, DATAFILE_INFO datafile)
{
	// printf("checking backwardStatus\n");

	for (int i = currentTimeframe; i > 0; --i)
	{
		if (polymerBondStatus[currentBondIndex][i].isItLoop == true)
		{
			backwardStatus.id = i;
			backwardStatus.isItLoop = polymerBondStatus[currentBondIndex][i].isItLoop; 
			backwardStatus.isItBridge = polymerBondStatus[currentBondIndex][i].isItBridge; 
			backwardStatus.isItDangle = polymerBondStatus[currentBondIndex][i].isItDangle; 
			backwardStatus.isItFree = polymerBondStatus[currentBondIndex][i].isItFree; 
			backwardStatus.loop_corrected = polymerBondStatus[currentBondIndex][i].loop_corrected; 
			backwardStatus.bridge_corrected = polymerBondStatus[currentBondIndex][i].bridge_corrected;

			break;
		}

		if (polymerBondStatus[currentBondIndex][i].isItBridge == true)
		{
			backwardStatus.id = i;
			backwardStatus.isItLoop = polymerBondStatus[currentBondIndex][i].isItLoop; 
			backwardStatus.isItBridge = polymerBondStatus[currentBondIndex][i].isItBridge; 
			backwardStatus.isItDangle = polymerBondStatus[currentBondIndex][i].isItDangle; 
			backwardStatus.isItFree = polymerBondStatus[currentBondIndex][i].isItFree; 
			backwardStatus.loop_corrected = polymerBondStatus[currentBondIndex][i].loop_corrected; 
			backwardStatus.bridge_corrected = polymerBondStatus[currentBondIndex][i].bridge_corrected;

			break;
		}
	}

	return backwardStatus;
}

BOND_STATUS checkingForwards (BOND_STATUS forwardStatus, BOND_STATUS **polymerBondStatus, int currentBondIndex, int currentTimeframe, DATAFILE_INFO datafile, int maxTimeframes)
{
	// printf("checking forwardStatus\n");

	for (int i = currentTimeframe; i < maxTimeframes; ++i)
	{
		if (polymerBondStatus[currentBondIndex][i].isItLoop == true)
		{
			forwardStatus.id = i;
			forwardStatus.isItLoop = polymerBondStatus[currentBondIndex][i].isItLoop; 
			forwardStatus.isItBridge = polymerBondStatus[currentBondIndex][i].isItBridge; 
			forwardStatus.isItDangle = polymerBondStatus[currentBondIndex][i].isItDangle; 
			forwardStatus.isItFree = polymerBondStatus[currentBondIndex][i].isItFree; 
			forwardStatus.loop_corrected = polymerBondStatus[currentBondIndex][i].loop_corrected; 
			forwardStatus.bridge_corrected = polymerBondStatus[currentBondIndex][i].bridge_corrected;

			break;
		}

		if (polymerBondStatus[currentBondIndex][i].isItBridge == true)
		{
			forwardStatus.id = i;
			forwardStatus.isItLoop = polymerBondStatus[currentBondIndex][i].isItLoop; 
			forwardStatus.isItBridge = polymerBondStatus[currentBondIndex][i].isItBridge; 
			forwardStatus.isItDangle = polymerBondStatus[currentBondIndex][i].isItDangle; 
			forwardStatus.isItFree = polymerBondStatus[currentBondIndex][i].isItFree; 
			forwardStatus.loop_corrected = polymerBondStatus[currentBondIndex][i].loop_corrected; 
			forwardStatus.bridge_corrected = polymerBondStatus[currentBondIndex][i].bridge_corrected;

			break;
		}
	}

	return forwardStatus;
}

BOND_STATUS **modifyDangles (BOND_STATUS **polymerBondStatus, BOND_STATUS forwardStatus, BOND_STATUS backwardStatus, DATAFILE_INFO datafile, int nTimeframes, int i)
{
	// If both the forward state and backward state are same
	if (forwardStatus.isItLoop == backwardStatus.isItLoop && forwardStatus.isItBridge == backwardStatus.isItBridge && forwardStatus.isItDangle == backwardStatus.isItDangle && forwardStatus.isItFree == backwardStatus.isItFree)
	{
		for (int j = backwardStatus.id; j < forwardStatus.id; ++j)
		{
			polymerBondStatus[i][j].isItLoop = forwardStatus.isItLoop;
			polymerBondStatus[i][j].isItBridge = forwardStatus.isItBridge;
			polymerBondStatus[i][j].isItDangle = forwardStatus.isItDangle;
			polymerBondStatus[i][j].isItFree = forwardStatus.isItFree;
			polymerBondStatus[i][j].loop_corrected = forwardStatus.loop_corrected;
			polymerBondStatus[i][j].bridge_corrected = forwardStatus.bridge_corrected;
		}
	}

	// If both the forward state and backward state are different
	if (forwardStatus.isItLoop != backwardStatus.isItLoop || forwardStatus.isItBridge != backwardStatus.isItBridge || forwardStatus.isItDangle != backwardStatus.isItDangle || forwardStatus.isItFree != backwardStatus.isItFree)
	{
		// if the forward state is zero (meaning, it is the last timeframe)
		// We don't have to worry about this condition, because these timeframes
		// will not count towards transition times

		// if the backward state is zero (meaning, it is the first timeframe)
		// Similar to the previous condition, this condition is not useful
		// to consider for transitions

		// if one state is bridge, while the other state is loop
		if ((forwardStatus.isItLoop == 1 && backwardStatus.isItBridge == 1) || 
			(forwardStatus.isItBridge == 1 && backwardStatus.isItLoop == 1))
		{
			for (int j = backwardStatus.id; j < forwardStatus.id; ++j)
			{
				polymerBondStatus[i][j].isItLoop = backwardStatus.isItLoop;
				polymerBondStatus[i][j].isItBridge = backwardStatus.isItBridge;
				polymerBondStatus[i][j].isItFree = backwardStatus.isItFree;
				polymerBondStatus[i][j].isItDangle = backwardStatus.isItDangle;
				polymerBondStatus[i][j].loop_corrected = backwardStatus.loop_corrected;
				polymerBondStatus[i][j].bridge_corrected = backwardStatus.bridge_corrected;
			}
		}
	}

	return polymerBondStatus;
}

BOND_STATUS **correctingDangles (BOND_STATUS **polymerBondStatus, DATAFILE_INFO datafile, int nTimeframes)
{
	BOND_STATUS forwardStatus, backwardStatus;

	forwardStatus.id = 0; forwardStatus.isItLoop = 0; forwardStatus.isItBridge = 0; forwardStatus.isItDangle = 0; forwardStatus.isItFree = 0; forwardStatus.loop_corrected = 0; forwardStatus.bridge_corrected = 0;

	backwardStatus.id = 0; backwardStatus.isItLoop = 0; backwardStatus.isItBridge = 0; backwardStatus.isItDangle = 0; backwardStatus.isItFree = 0; backwardStatus.loop_corrected = 0; backwardStatus.bridge_corrected = 0;

	printf("\nCorrecting dangles...\n");

/*	for (int i = 0; i < datafile.nBonds; ++i)
	{
		printStateNumber (polymerBondStatus[i][0].isItFree, polymerBondStatus[i][0].isItDangle, polymerBondStatus[i][0].isItLoop, polymerBondStatus[i][0].isItBridge);
		usleep (10000);
	}
	sleep (1000);
*/
	for (int i = 0; i < datafile.nBonds; ++i)
	{
		printf("%d/%d                                        \r", i, datafile.nBonds);
		fflush (stdout);

		for (int j = 0; j < nTimeframes; ++j)
		{
			// scroll through the time for every bond
			// once a dangle is found, go backward and go forward in time
			// while going backwards and forwards: store the index of the first instance of a bridge/loop

			// printf("%d => %d %d %d %d\n", j, polymerBondStatus[i][j].isItLoop, polymerBondStatus[i][j].isItBridge, polymerBondStatus[i][j].isItDangle, polymerBondStatus[i][j].isItFree);

			// printStateNumber (polymerBondStatus[i][j].isItFree, polymerBondStatus[i][j].isItDangle, polymerBondStatus[i][j].isItLoop, polymerBondStatus[i][j].isItBridge);

			if (polymerBondStatus[i][j].isItDangle == true)
			{
				// reset the status
				// in case a forward/backward loop/bridge is not found, then the status will be '0'
				// if the status is '0', then don't switch the dangles to any other state
				forwardStatus.id = 0; forwardStatus.isItLoop = 0; forwardStatus.isItBridge = 0; forwardStatus.isItDangle = 0; forwardStatus.isItFree = 0; forwardStatus.loop_corrected = 0; forwardStatus.bridge_corrected = 0;

				backwardStatus.id = 0; backwardStatus.isItLoop = 0; backwardStatus.isItBridge = 0; backwardStatus.isItDangle = 0; backwardStatus.isItFree = 0; backwardStatus.loop_corrected = 0; backwardStatus.bridge_corrected = 0;

				// printf("dangle found...\n");
				// go backward in time
				backwardStatus = checkingBackwards (backwardStatus, polymerBondStatus, i, j, datafile);

				// go forward in time
				forwardStatus = checkingForwards (forwardStatus, polymerBondStatus, i, j, datafile, nTimeframes);

				// printf("status => from %d (%d %d %d %d) to %d (%d %d %d %d)\n", backwardStatus.id, backwardStatus.isItLoop, backwardStatus.isItBridge, backwardStatus.isItDangle, backwardStatus.isItFree, forwardStatus.id, forwardStatus.isItLoop, forwardStatus.isItBridge, forwardStatus.isItDangle, forwardStatus.isItBridge);

				polymerBondStatus = modifyDangles (polymerBondStatus, forwardStatus, backwardStatus, datafile, nTimeframes, i);
			}

			// printStateNumber (polymerBondStatus[i][j].isItFree, polymerBondStatus[i][j].isItDangle, polymerBondStatus[i][j].isItLoop, polymerBondStatus[i][j].isItBridge);
			// printf("\n");
			// usleep (100000);

			// printf("%d => %d %d %d %d\n", j, polymerBondStatus[i][j].isItLoop, polymerBondStatus[i][j].isItBridge, polymerBondStatus[i][j].isItDangle, polymerBondStatus[i][j].isItFree);
			// usleep (100000);

			// if the state of the bond in forward/backward instances are the same, then convert dangles
			// if the states are different, then keep the dangles.
		}
	}

	return polymerBondStatus;
}

int countLBtransitions (int nLB, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	nLB = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (polymerBondStatus[i][j - 1].isItLoop == 1 && polymerBondStatus[i][j].isItBridge == 1)
			{
				nLB++;
			}
		}
	}

	return nLB;
}

int countBBtransitions (int nBB, BOND_STATUS **polymerBondStatus, BOUND_STATUS **beadBoundStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	nBB = 0;

	bool newDangle = 0;
	int boundParticleID1 = 0, boundParticleID2 = 0;

	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;

	int currentlyBridge = 0, dangleFromBridge = 0;

	int atomNumber1, atomNumber2;
	int particleID1 = 0, particleID2 = 0;

	int debugg = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		newDangle = 0;
		currentlyBridge = 0;
		dangleFromBridge = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				fprintf(stdout, "\n %d (%d %d) + %d (%d %d) => %d  [%d/%d]", 
					atomNumber1, 
					beadBoundStatus[atomNumber1 - 1][j].isItBound, 
					polymerBondStatus[i][j].id1, 
					atomNumber2, 
					beadBoundStatus[atomNumber2 - 1][j].isItBound, 
					polymerBondStatus[i][j].id2, 
					nBB,
					currentlyBridge,
					dangleFromBridge);
				fflush (stdout);

				if (nBB%50 == 0 && nBB > 0)
				{
					usleep (5000);
				}
			}

			// if it is already in a dangle state,
			// this checks if the dangle becomes a bridge or a loop
			if (dangleFromBridge == 1)
			{
				// dangle becomes a bridge
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 != polymerBondStatus[i][j].id2))
				{
					// checks if the newly formed bridge is different from the previous bridge
					if (((particleID1 != polymerBondStatus[i][j].id1 || particleID2 != polymerBondStatus[i][j].id2) &&
						(particleID2 != polymerBondStatus[i][j].id1 || particleID1 != polymerBondStatus[i][j].id2)) &&
						(polymerBondStatus[i][j].id1 != 0) &&
						(polymerBondStatus[i][j].id2 != 0))
					{
						nBB++;

						if (debugg == 1)
						{
							fprintf(stdout, ">>>Transition complete...");
							fprintf(stdout, "%d %d, %d %d\n", particleID1, polymerBondStatus[i][j].id1, particleID2, polymerBondStatus[i][j].id2);
							fflush (stdout);
							sleep (10);
						}

						currentlyBridge = 0;
						dangleFromBridge = 0;
					}
					else
					{
						currentlyBridge = 1;
						dangleFromBridge = 0;

						if (debugg == 1)
						{
							fprintf(stdout, ">Bridge formed with same particles...");
							fprintf(stdout, "%d %d, %d %d\n", particleID1, polymerBondStatus[i][j].id1, particleID2, polymerBondStatus[i][j].id2);
							fflush (stdout);
							sleep (5);
						}
					}
				}

				// dangle becomes a loop
				else if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 == polymerBondStatus[i][j].id2))
				{
					currentlyBridge = 0;
					dangleFromBridge = 0;

					if (debugg == 1)
					{
						fprintf(stdout, ">Loop formed from dangle...");
						fflush (stdout);
						// sleep (2);
					}
				}
			}

			// If it is already in a bridge state,
			// then this 'if' condition checks if it transition to a dangle state
			if ((currentlyBridge == 1) &&
				(dangleFromBridge == 0) &&
				((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) ||
				(beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyBridge = 0;
				dangleFromBridge = 1;

				if (debugg == 1)
				{
					fprintf(stdout, ">Dangle formed from bridge...");
					fflush (stdout);
					// sleep (2);
				}
			}

			// This loop checks if the dumbbell is in a bridge state
			if ((polymerBondStatus[i][j].isItBridge == 1) && (currentlyBridge == 0) && (dangleFromBridge == 0))
			{
				currentlyBridge = 1;
				dangleFromBridge = 0;

				particleID1 = polymerBondStatus[i][j].id1;
				particleID2 = polymerBondStatus[i][j].id2;

				if (debugg == 1)
				{
					fprintf(stdout, ">Bridge detected...");
					fflush (stdout);
					// sleep (2);
				}
			}


/*			// bridge becomes a dangle
			if (polymerBondStatus[i][j - 1].isItBridge == 1 && polymerBondStatus[i][j].isItDangle == 1) {
				boundParticleID1 = polymerBondStatus[i][j - 1].id1;
				boundParticleID2 = polymerBondStatus[i][j - 1].id2;
				newDangle = 1; }
			else if (polymerBondStatus[i][j - 1].isItBridge == true && polymerBondStatus[i][j].isItLoop == true) {
				newDangle = 0; }
			else if (polymerBondStatus[i][j - 1].isItBridge == true && polymerBondStatus[i][j].isItFree == true) {
				newDangle = 0; }

			// dangle becomes a bridge again
			// it is counted as a bridge to bridge transition only if the 
			// bridge is formed between two different set of micelles
			if (polymerBondStatus[i][j].isItBridge == 1 && newDangle == 1) 
			{
				newDangle = 0;

				if ((boundParticleID1 != polymerBondStatus[i][j].id1) || (boundParticleID2 != polymerBondStatus[i][j].id2))
				{
					nBB++; 
				}
			}

			// dangle becomes a loop or free chain
			if ((polymerBondStatus[i][j].isItFree == true || polymerBondStatus[i][j].isItLoop == true) && newDangle == 1) {
				newDangle = 0; }*/
		}
	}

	fprintf(stdout, "nBB: %d\n", nBB);
	fflush (stdout);

	return nBB;
}

float *initTauBB (float *tauBB, int nBB)
{
	for (int i = 0; i < nBB; ++i)
	{
		tauBB[i] = 0;
	}

	return tauBB;
}

float *countTauBB (float *tauBB, BOND_STATUS **polymerBondStatus, BOUND_STATUS **beadBoundStatus, int nTimeframes, DATAFILE_INFO datafile, int nBB)
{
	bool newDangle = 0;

	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;

	int currentlyBridge = 0, dangleFromBridge = 0;
	int atomNumber1, atomNumber2;

	int counter = 0;
	int debugg = 0;

	// initialize tauBB
	tauBB = initTauBB (tauBB, nBB);
	int currentIndex = 0;
	int particleID1 = 0, particleID2 = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		newDangle = 0;
		counter = 0;
		currentlyBridge = 0;
		dangleFromBridge = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				fprintf(stdout, "\n %d (%d %d) + %d (%d %d) => %d  [%d/%d]", 
					atomNumber1, 
					beadBoundStatus[atomNumber1 - 1][j].isItBound, 
					polymerBondStatus[i][j].id1, 
					atomNumber2, 
					beadBoundStatus[atomNumber2 - 1][j].isItBound, 
					polymerBondStatus[i][j].id2, 
					counter,
					currentlyBridge,
					dangleFromBridge);
				fflush (stdout);

				if (counter%50 == 0 && counter > 0)
				{
					usleep (5000);
				}
			}

			if (dangleFromBridge == 1)
			{
				// dangle becomes a bridge
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 != polymerBondStatus[i][j].id2) &&
					(polymerBondStatus[i][j].id1 > 0) &&
					(polymerBondStatus[i][j].id2 > 0))
				{
					// checks if the newly formed bridge is different from the previous bridge
					if (((particleID1 != polymerBondStatus[i][j].id1 || particleID2 != polymerBondStatus[i][j].id2) &&
						(particleID2 != polymerBondStatus[i][j].id1 || particleID1 != polymerBondStatus[i][j].id2)) &&
						(polymerBondStatus[i][j].id1 != 0) &&
						(polymerBondStatus[i][j].id2 != 0))
					{
						tauBB[currentIndex] = counter;

						currentIndex++;

						if (debugg == 1)
						{
							fprintf(stdout, ">>>Transition complete...");
							fprintf(stdout, "%d %d, %d %d\n", particleID1, polymerBondStatus[i][j].id1, particleID2, polymerBondStatus[i][j].id2);
							fflush (stdout);
							sleep (10);
						}

						currentlyBridge = 0;
						dangleFromBridge = 0;
						counter = 0;
					}
					// if it goes back to the same bridge
					else
					{
						currentlyBridge = 1;
						dangleFromBridge = 0;
						counter = 0;

						if (debugg == 1)
						{
							fprintf(stdout, ">Bridge formed with same particles...");
							fprintf(stdout, "%d %d, %d %d\n", particleID1, polymerBondStatus[i][j].id1, particleID2, polymerBondStatus[i][j].id2);
							fflush (stdout);
							sleep (5);
						}
					}
				}

				// dangle becomes a loop
				else if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 == polymerBondStatus[i][j].id2))
				{
					currentlyBridge = 0;
					dangleFromBridge = 0;
					counter = 0;

					if (debugg == 1)
					{
						fprintf(stdout, ">Loop formed from dangle...");
						fflush (stdout);
						// sleep (2);
					}
				}
			}

			if ((currentlyBridge == 1) &&
				(dangleFromBridge == 0) &&
				((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) ||
				(beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyBridge = 0;
				dangleFromBridge = 1;

				if (debugg == 1)
				{
					fprintf(stdout, ">Dangle formed from bridge...");
					fflush (stdout);
					// sleep (2);
				}
			}

			if ((polymerBondStatus[i][j].isItBridge == 1) && (currentlyBridge == 0) && (dangleFromBridge == 0))
			{
				currentlyBridge = 1;
				dangleFromBridge = 0;

				particleID1 = polymerBondStatus[i][j].id1;
				particleID2 = polymerBondStatus[i][j].id2;

				if (debugg == 1)
				{
					fprintf(stdout, ">Bridge detected...");
					fflush (stdout);
					// sleep (2);
				}
			}

			if ((currentlyBridge == 1) || (dangleFromBridge == 1))
			{
				counter++;

				if (debugg == 1)
				{
					if (currentlyBridge == 1)
					{
						fprintf(stdout, "  [bridge]");
						fflush (stdout);
					}
					else if (dangleFromBridge == 1)
					{
						fprintf(stdout, "  [dangle]");
						fflush (stdout);
					}
				}
			}

/*			// Bridge becomes a dangle, I am storing the IDs of bound core/ghost particles
			// And the variable 'newDangle' is activated
			if (polymerBondStatus[i][j - 1].isItBridge == 1 && polymerBondStatus[i][j].isItDangle == 1) {
				boundParticleID1 = polymerBondStatus[i][j - 1].id1;
				boundParticleID2 = polymerBondStatus[i][j - 1].id2;
				newDangle = 1; }
			else if (polymerBondStatus[i][j - 1].isItBridge == true && polymerBondStatus[i][j].isItLoop == true) {
				newDangle = 0; }
			else if (polymerBondStatus[i][j - 1].isItBridge == true && polymerBondStatus[i][j].isItFree == true) {
				newDangle = 0; }
			else if (polymerBondStatus[i][j].isItLoop == 1) {
				newDangle = 0; }
			else if (polymerBondStatus[i][j].isItFree == 1) {
				newDangle = 0; }

			if (newDangle == 1 || polymerBondStatus[i][j].isItBridge == 1) {
				tauBB[currentBridge]++; }

			// dangle becomes a bridge again
			if (polymerBondStatus[i][j].isItBridge == 1 && newDangle == 1) 
			{
				newDangle = 0; 

				if ((boundParticleID1 != polymerBondStatus[i][j].id1) || (boundParticleID2 != polymerBondStatus[i][j].id2)) {
					currentBridge++; }
				else {
					tauBB[currentBridge] = 0; }
			}

			// dangle becomes a loop or free chain
			if ((polymerBondStatus[i][j].isItFree == true || polymerBondStatus[i][j].isItLoop == true) && newDangle == 1) {
				tauBB[currentBridge] = 0; newDangle = 0; }*/
		}
	}

	return tauBB;
}

int countLLtransitions (int nLL, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, BOUND_STATUS **beadBoundStatus, DATA_ATOMS *sortedAtoms)
{
	nLL = 0;
	int counter = 0, currentIndex = 0;
	int bondNumber;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	int atomNumber1, atomNumber2;

	int debugg = 0;

	printf("\nCounting loop to loop transitions...\n");

	for (int i = 0; i < datafile.nBonds; ++i)
	{		
		counter = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (beadBoundStatus[atomNumber1 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j - 1].isItBound == 1 && polymerBondStatus[i][j - 1].isItLoop == 1)
			{
				if (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 || beadBoundStatus[atomNumber2 - 1][j].isItBound == 0)
				{
					counter = 1; // here, counter is only to identify the state. In the next connected function, counter will be used to measure the transition time.
				}
			}

			if (debugg)
			{
				printf("%d (%d: %d/%d) + %d (%d: %d/%d) => %d, loop: %d, dangle: %d, timeframe: %d  ", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, polymerBondStatus[i][j].id2, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id1, polymerBondStatus[i][j].id2, counter, polymerBondStatus[i][j].isItLoop, polymerBondStatus[i][j].isItDangle, j);

				if (counter == 1)
				{
					usleep (100000);
				}
			}

			if ((counter == 1) && 
				(polymerBondStatus[i][j].isItLoop == 1) && 
				((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) || 
				(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				counter = 0;
				nLL++;

				if (debugg)
				{
					fprintf(stdout, ">counter stopped<");
					fflush (stdout);
					sleep (10);
				}
			}

			if (debugg)
			{
				printf("\n");
			}
		}
	}

	return nLL;
}

float *countTauLL (float *tauLL, int nLL, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, BOUND_STATUS **beadBoundStatus, DATA_ATOMS *sortedAtoms)
{
	int counter = 0, currentIndex = 0;
	int bondNumber;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	int atomNumber1, atomNumber2;

	int transitionToDangle = 0;

	int debugg = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;
		transitionToDangle = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1 && polymerBondStatus[i][j].isItLoop == 1)
			{
				counter++;

/*				if ((beadBoundStatus[atomNumber1 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber1 - 1][j].isItBound == 0) || 
					(beadBoundStatus[atomNumber2 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0))
				{
					counter++; // here, counter is only to identify the state. In the next connected function, counter will be used to measure the transition time.
				}*/
			}

			if (counter > 0)
			{
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) || 
					(beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0))
				{
					transitionToDangle = 1;
					counter++;
				}
			}

			if (debugg == 1)
			{
				fprintf(stdout, "%d (%d) + %d (%d) => %d\n", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, counter);
				usleep (100000);
			}

			if ((transitionToDangle == 1) && 
				(polymerBondStatus[i][j].isItLoop == 1) && 
				((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) || 
				(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				tauLL[currentIndex] = counter;
				currentIndex++;

				if (debugg == 1)
				{
					printf("tauLL: %d\n", counter);
					usleep (100000);
				}

				counter = 0;
				transitionToDangle = 0;
			}
		}
	}

	return tauLL;
}

int countBLtransitions (int nBL, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	nBL = 0;
	printf("\nCounting bridge to loop transitions...\n");

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (polymerBondStatus[i][j - 1].isItBridge == 1 && polymerBondStatus[i][j].isItLoop == 1)
			{
				nBL++;
			}
		}
	}

	return nBL;
}

float *initTauBL (int nBL, float *tauBL)
{
	for (int i = 0; i < nBL; ++i)
	{
		tauBL[i] = 0;
	}

	return tauBL;
}

float *countTauBL2 (float *tauBL, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	int currentTransition = 0, currentTau = 0;
	int boundParticleID1 = 0, boundParticleID2 = 0;
	int bbcorrection = 1;
	int debugg = 1;
	int intermediateBridge = 0;
	int tempPrint = 0;

	printf("\nCalculating transition time from bridge to loop (version 2.0)...\n");

	if (bbcorrection == 1)
	{
		fprintf(stdout, "In tauBL2.0, only the bridge to loop transitions with a bridge as an intermediate are counted. The rest of the B -> L transitions are discarded...\n");
	}

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentTau = 0;
		intermediateBridge = 0;

		if (debugg == 1)
		{
			printf("\nNext bond\n\n");
			sleep (3);
		}

		boundParticleID1 = 0;
		boundParticleID2 = 0;

		for (int j = 0; j < nTimeframes; ++j)
		{
			if (debugg == 1)
			{
				printf("\n");
				printStateNumber (polymerBondStatus[i][j].isItFree, polymerBondStatus[i][j].isItDangle, polymerBondStatus[i][j].isItLoop, polymerBondStatus[i][j].isItBridge);
				printf("  %d %d  ", polymerBondStatus[i][j].id1, polymerBondStatus[i][j].id2);
				usleep (100000);
			}

			if (j > 0)
			{
				// if it becomes a loop from an intermediate bridge
				if (currentTau > 0 && polymerBondStatus[i][j].isItLoop == 1 && intermediateBridge == 1)
				{
					tauBL[currentTransition] = (float)currentTau;
					currentTransition++;
					currentTau = 0;
					intermediateBridge = 0;

					if (debugg == 1)
					{
						printf("--> %d\n", currentTransition);
						sleep (5);
					}
				}
				// if it becomes a loop without an intermediate bridge state, (i.e) direct bridge to loop
				else if (currentTau > 0 && polymerBondStatus[i][j].isItLoop == 1 && intermediateBridge == 0)
				{
					currentTau = 0;
					intermediateBridge = 0;

					if (debugg == 1)
					{
						printf("timer discarded\n");
					}
				}
			}

			// counter starts if the dumbbell is in a bridge state
			// counter continues if it a bridge
			if (polymerBondStatus[i][j].isItBridge == 1)
			{
				currentTau++;

				if (debugg == 1)
				{
					printf("(%d) \n", currentTau);
					usleep (10000);
				}

				// turns on the 'intermediateBridge' if it detects a bridge to bridge transition
				if ((boundParticleID1 != 0) && (boundParticleID2 != 0) && (polymerBondStatus[i][j].id1 != 0) && (polymerBondStatus[i][j].id2 != 0) && (bbcorrection == 1) && (currentTau > 1))
				{
					if ((boundParticleID1 != polymerBondStatus[i][j].id1) || (boundParticleID2 != polymerBondStatus[i][j].id2))
					{
						if ((boundParticleID1 != polymerBondStatus[i][j].id2) && (boundParticleID2 != polymerBondStatus[i][j].id1))
						{
							intermediateBridge = 1;

							if (debugg == 1)
							{
								printf("\n[b2b] %d -> %d; %d -> %d\n", boundParticleID1, polymerBondStatus[i][j].id1, boundParticleID2, polymerBondStatus[i][j].id2);
							}
						}
					}
				}

				if (polymerBondStatus[i][j].id1 != 0 && polymerBondStatus[i][j].id2 != 0)
				{
					if ((polymerBondStatus[i][j].id1 != boundParticleID1) && (debugg == 1))
					{
						printf("\nboundParticleID1: %d -> ", boundParticleID1);
						tempPrint = 1;
					}
					
					boundParticleID1 = polymerBondStatus[i][j].id1;

					if (tempPrint == 1)
					{
						printf(" %d [state: %d/%d]\n", boundParticleID1, polymerBondStatus[i][j].isItBridge, polymerBondStatus[i][j].isItLoop);
						tempPrint = 0;
					}
				}

				if (polymerBondStatus[i][j].id1 != 0 && polymerBondStatus[i][j].id2 != 0)
				{
					if ((polymerBondStatus[i][j].id2 != boundParticleID2) && (debugg == 1))
					{
						printf("\nboundParticleID2: %d -> ", boundParticleID2);
						tempPrint = 1;
					}

					boundParticleID2 = polymerBondStatus[i][j].id2;

					if (tempPrint == 1)
					{
						printf(" %d [state: %d/%d]\n", boundParticleID2, polymerBondStatus[i][j].isItBridge, polymerBondStatus[i][j].isItLoop);
						tempPrint = 0;
					}
				}
			}

			else if (currentTau > 0 && polymerBondStatus[i][j].isItLoop == 0 && polymerBondStatus[i][j].isItFree == 0)
			{
				currentTau++;
			}
		}

		if (debugg == 1)
		{
			printf("\n");
		}
	}

	return tauBL;
}

float *countTauBL (float *tauBL, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, int bbcorrection)
{
	int currentTransition = 0, currentTau = 0;
	int boundParticleID1 = 0, boundParticleID2 = 0;
	int tempPrint = 0;
	int debugg = 0;

	printf("\nCalculating transition time from bridge to loop...\n");

	if (bbcorrection == 1)
	{
		fprintf(stdout, "Correcting the bridge to loop transitions based on bridge to bridge transitions...\n");
	}
	else if (bbcorrection == 0)
	{
		fprintf(stdout, "Not correcting for bridge to bridge transitions in tauBL...\n");
	}
	else
	{
		fprintf(stdout, "Undefined behaviour for bridge to bridge correction...\n");
	}

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentTau = 0;
		// printf("\nNext bond\n\n");
		// sleep (5);

		boundParticleID1 = 0;
		boundParticleID2 = 0;

		for (int j = 0; j < nTimeframes; ++j)
		{

			// printStateNumber (polymerBondStatus[i][j].isItFree, polymerBondStatus[i][j].isItDangle, polymerBondStatus[i][j].isItLoop, polymerBondStatus[i][j].isItBridge);
			// usleep (100000);

			if (polymerBondStatus[i][j].isItBridge == 1)
			{
				currentTau++;

				if (debugg == 1)
				{
					printf("(%d) \n", currentTau);
					usleep (10000);
				}

				// turns on the 'intermediateBridge' if it detects a bridge to bridge transition
				if ((boundParticleID1 != 0) && (boundParticleID2 != 0) && (polymerBondStatus[i][j].id1 != 0) && (polymerBondStatus[i][j].id2 != 0) && (bbcorrection == 1) && (currentTau > 1))
				{
					if ((boundParticleID1 != polymerBondStatus[i][j].id1) || (boundParticleID2 != polymerBondStatus[i][j].id2))
					{
						if ((boundParticleID1 != polymerBondStatus[i][j].id2) && (boundParticleID2 != polymerBondStatus[i][j].id1))
						{
							currentTau = 0;

							if (debugg == 1)
							{
								printf("\n[b2b] %d -> %d; %d -> %d\n", boundParticleID1, polymerBondStatus[i][j].id1, boundParticleID2, polymerBondStatus[i][j].id2);
							}
						}
					}
				}

				if (polymerBondStatus[i][j].id1 != 0 && polymerBondStatus[i][j].id2 != 0)
				{
					if ((polymerBondStatus[i][j].id1 != boundParticleID1) && (debugg == 1))
					{
						printf("\nboundParticleID1: %d -> ", boundParticleID1);
						tempPrint = 1;
					}
					
					boundParticleID1 = polymerBondStatus[i][j].id1;

					if (tempPrint == 1)
					{
						printf(" %d [state: %d/%d]\n", boundParticleID1, polymerBondStatus[i][j].isItBridge, polymerBondStatus[i][j].isItLoop);
						tempPrint = 0;
					}
				}

				if (polymerBondStatus[i][j].id1 != 0 && polymerBondStatus[i][j].id2 != 0)
				{
					if ((polymerBondStatus[i][j].id2 != boundParticleID2) && (debugg == 1))
					{
						printf("\nboundParticleID2: %d -> ", boundParticleID2);
						tempPrint = 1;
					}

					boundParticleID2 = polymerBondStatus[i][j].id2;

					if (tempPrint == 1)
					{
						printf(" %d [state: %d/%d]\n", boundParticleID2, polymerBondStatus[i][j].isItBridge, polymerBondStatus[i][j].isItLoop);
						tempPrint = 0;
					}
				}

				// Old code block
				/*
				currentTau++;
				// printf("(%d) ", currentTau);
				// usleep (10000);

				if ((boundParticleID1 != 0) && (boundParticleID2 != 0) && (polymerBondStatus[i][j].id1 != 0) && (polymerBondStatus[i][j].id2 != 0) && (bbcorrection == 1))
				{
					if ((boundParticleID1 != polymerBondStatus[i][j].id1) || (boundParticleID2 != polymerBondStatus[i][j].id2))
					{
						currentTau = 0;
					}
				}

				if (polymerBondStatus[i][j].id1 != 0)
				{
					boundParticleID1 = polymerBondStatus[i][j].id1;
				}

				if (polymerBondStatus[i][j].id2 != 0)
				{
					boundParticleID2 = polymerBondStatus[i][j].id2;
				}
				*/
			}
			else if (currentTau > 0 && polymerBondStatus[i][j].isItLoop == 0 && polymerBondStatus[i][j].isItFree == 0)
			{
				currentTau++;
			}
			else if (j > 0)
			{
				if (currentTau > 0 && polymerBondStatus[i][j].isItLoop == 1)
				{
					tauBL[currentTransition] = (float)currentTau;
					currentTransition++;
					currentTau = 0;
					// printf("--> %d\n", currentTransition);
					// sleep (1);
				}
			}
		}

		// printf("\n");
	}

	return tauBL;
}

float *countTauLB (float *tauLB, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	int currentTransition = 0, currentTau = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentTau = 0;

		for (int j = 0; j < nTimeframes; ++j)
		{
			if (polymerBondStatus[i][j].isItLoop == 1)
			{
				currentTau++;
			}
			else if (j > 0 && polymerBondStatus[i][j].isItBridge == 1)
			{
				if (polymerBondStatus[i][j - 1].isItLoop == 1)
				{
					tauLB[currentTransition] = (float)currentTau;
					currentTransition++;
					currentTau = 0;
				}
			}
		}
	}

	return tauLB;
}

void printTauLL (float *tauLL, int nLL, const char *folderName)
{
	char *tauLL_file_filename;
	tauLL_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauLL_file_filename, 500, "%s/tauLL.output", folderName);
	fprintf(stdout, "Printing %s\n", tauLL_file_filename);

	FILE *tauLL_file;
	tauLL_file = fopen (tauLL_file_filename, "w");

	for (int i = 0; i < nLL; ++i)
	{
		fprintf(tauLL_file, "%f\n", tauLL[i]);
	}

	fclose (tauLL_file);
}

void printTauBL (float *tauBL, int nBL, const char *folderName)
{
	char *tauBL_file_filename;
	tauBL_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauBL_file_filename, 500, "%s/tauBL.output", folderName);
	fprintf(stdout, "Printing %s\n", tauBL_file_filename);

	FILE *tauBL_file;
	tauBL_file = fopen (tauBL_file_filename, "w");

	for (int i = 0; i < nBL; ++i)
	{
		fprintf(tauBL_file, "%f\n", tauBL[i]);
	}

	fclose (tauBL_file);
}

void printTauBL2 (float *tauBL, int nBL, const char *folderName)
{
	char *tauBL_file_filename;
	tauBL_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauBL_file_filename, 500, "%s/tauBL2.output", folderName);
	fprintf(stdout, "Printing %s\n", tauBL_file_filename);

	FILE *tauBL_file;
	tauBL_file = fopen (tauBL_file_filename, "w");

	for (int i = 0; i < nBL; ++i)
	{
		fprintf(tauBL_file, "%f\n", tauBL[i]);
	}

	fclose (tauBL_file);
}

void printTauLB (float *tauLB, int nLB, const char *folderName)
{
	char *tauLB_file_filename;
	tauLB_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauLB_file_filename, 500, "%s/tauLB.output", folderName);
	fprintf(stdout, "Printing %s\n", tauLB_file_filename);

	FILE *tauLB_file;
	tauLB_file = fopen (tauLB_file_filename, "w");

	for (int i = 0; i < nLB; ++i)
	{
		fprintf(tauLB_file, "%f\n", tauLB[i]);
	}

	fclose (tauLB_file);
}

void printTauBB (float *tauBB, int nBB, const char *folderName)
{
	char *tauBB_file_filename;
	tauBB_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauBB_file_filename, 500, "%s/tauBB.output", folderName);
	fprintf (stdout, "Printing %s\n", tauBB_file_filename);
	fflush (stdout);

	FILE *tauBB_file;
	tauBB_file = fopen (tauBB_file_filename, "w");

	for (int i = 0; i < nBB; ++i)
	{
		fprintf(tauBB_file, "%f\n", tauBB[i]);
	}

	fclose (tauBB_file);
}

void printTauE (float *tauLB, int nLB, const char *folderName)
{
	char *tauLB_file_filename;
	tauLB_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauLB_file_filename, 500, "%s/tauE.output", folderName);
	fprintf(stdout, "Printing %s\n", tauLB_file_filename);

	FILE *tauLB_file;
	tauLB_file = fopen (tauLB_file_filename, "w");

	for (int i = 0; i < nLB; ++i)
	{
		fprintf(tauLB_file, "%f\n", tauLB[i]);
	}

	fclose (tauLB_file);
}

void printTauEm (float *tauLB, int nLB, const char *folderName)
{
	char *tauLB_file_filename;
	tauLB_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauLB_file_filename, 500, "%s/tauEm.output", folderName);
	fprintf(stdout, "Printing %s\n", tauLB_file_filename);

	FILE *tauLB_file;
	tauLB_file = fopen (tauLB_file_filename, "w");

	for (int i = 0; i < nLB; ++i)
	{
		fprintf(tauLB_file, "%f\n", tauLB[i]);
	}

	fclose (tauLB_file);
}

void printTauS (float *tauLB, int nLB, const char *folderName)
{
	char *tauLB_file_filename;
	tauLB_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauLB_file_filename, 500, "%s/tauS.output", folderName);
	fprintf(stdout, "Printing %s\n", tauLB_file_filename);

	FILE *tauLB_file;
	tauLB_file = fopen (tauLB_file_filename, "w");

	for (int i = 0; i < nLB; ++i)
	{
		fprintf(tauLB_file, "%f\n", tauLB[i]);
	}

	fclose (tauLB_file);
}

void printTauSm (float *tauLB, int nLB, const char *folderName)
{
	char *tauLB_file_filename;
	tauLB_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauLB_file_filename, 500, "%s/tauSm.output", folderName);
	fprintf(stdout, "Printing %s\n", tauLB_file_filename);

	FILE *tauLB_file;
	tauLB_file = fopen (tauLB_file_filename, "w");

	for (int i = 0; i < nLB; ++i)
	{
		fprintf(tauLB_file, "%f\n", tauLB[i]);
	}

	fclose (tauLB_file);
}

void computeStats (float *average, float *stdev, float *stderr, float *inputData, int nLines)
{
	(*average) = 0;
	(*stderr) = 0;
	(*stdev) = 0;

	for (int i = 0; i < nLines; ++i)
	{
		(*average) += inputData[i];
	}

	(*average) /= nLines;

	for (int i = 0; i < nLines; ++i)
	{
		(*stdev) += ((inputData[i] - (*average)) * (inputData[i] - (*average)));
	}

	(*stdev) /= nLines;
	(*stdev) = sqrt ((*stdev));
	(*stdev) /= sqrt (nLines);
}

// tauE is the time it takes for a bridge to eject
// tauEm is the time it takes for a bridge/loop to eject
// tauS is the time the bead spends in the solvent as a dangle before returning to bridge
// tauSm is the time the bead spends in the solvent as a dangle before returning to bridge/loop

float *initTauE (float *tauE, int nE)
{
	for (int i = 0; i < nE; ++i)
	{
		tauE[i] = 0;
	}

	return tauE;
}

float *countTauE (float *tauE, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, int nE)
{
	int currentIndex = 0, counter = 0;
	int justEjected = 0;
	int boundParticleID1, boundParticleID2;

	tauE = initTauE (tauE, nE);

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;
		justEjected = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			// this loop checks if the bead has just recently ejected from a micelle core
			// in addition, the first 'if' loop checks if the transtion is from bridge to dangle
			if (polymerBondStatus[i][j - 1].isItBridge == 1 && polymerBondStatus[i][j].isItDangle == 1)
			{
				if ((polymerBondStatus[i][j - 1].id1 > polymerBondStatus[i][j].id1) || (polymerBondStatus[i][j - 1].id2 > polymerBondStatus[i][j].id2))
				{
					if ((polymerBondStatus[i][j - 1].id1*polymerBondStatus[i][j].id1 == 0) || (polymerBondStatus[i][j - 1].id2*polymerBondStatus[i][j].id2 == 0))
					{
						justEjected = 1;
						boundParticleID1 = polymerBondStatus[i][j].id1;
						boundParticleID2 = polymerBondStatus[i][j].id2;
					}
				}
			}

			if (justEjected == 1)
			{
				counter++;
			}

			// this loop checks if the bead gets re-attaches to some other micelle core
			// the 'if' loop checks if the transition is from dangle to bridge
			if (justEjected == 1 && polymerBondStatus[i][j - 1].isItDangle && polymerBondStatus[i][j].isItBridge)
			{
				if ((polymerBondStatus[i][j - 1].id1 < polymerBondStatus[i][j].id1) || (polymerBondStatus[i][j - 1].id2 < polymerBondStatus[i][j].id2))
				{
					if ((polymerBondStatus[i][j - 1].id1*polymerBondStatus[i][j].id1 == 0) || (polymerBondStatus[i][j - 1].id2*polymerBondStatus[i][j].id2 == 0))
					{
						justEjected = 0;
						
						// if the bead is getting attached to a different micelle core, then the counter is stored
						// if it re-attaches to the same micelle core, then the counter is rejected
						if (boundParticleID1 != polymerBondStatus[i][j].id1 || boundParticleID2 != polymerBondStatus[i][j].id2)
						{
							tauE[currentIndex] = counter;
							counter = 0;
							currentIndex++;
						}
					}
				}
			}
		}
	}

	return tauE;
}

int countEtransitions (int nE, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	nE = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (polymerBondStatus[i][j - 1].isItBridge == 1 && polymerBondStatus[i][j].isItDangle == 1)
			{
				nE++;
			}
		}
	}

	return nE;
}

float *initTauEm (float *tauEm, int nEm)
{
	for (int i = 0; i < nEm; ++i)
	{
		tauEm[i] = 0;
	}

	return tauEm;
}

float *countTauEm (float *tauEm, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, int nEm)
{
	int currentIndex = 0, counter = 0;
	int justEjected = 0;
	int boundParticleID1, boundParticleID2;

	tauEm = initTauEm (tauEm, nEm);

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;
		justEjected = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			// this loop checks if the bead has just recently ejected from a micelle core
			if ((polymerBondStatus[i][j - 1].id1 > polymerBondStatus[i][j].id1) || (polymerBondStatus[i][j - 1].id2 > polymerBondStatus[i][j].id2))
			{
				if ((polymerBondStatus[i][j - 1].id1*polymerBondStatus[i][j].id1 == 0) || (polymerBondStatus[i][j - 1].id2*polymerBondStatus[i][j].id2 == 0))
				{
					justEjected = 1;
					boundParticleID1 = polymerBondStatus[i][j].id1;
					boundParticleID2 = polymerBondStatus[i][j].id2;
				}
			}

			if (justEjected == 1)
			{
				counter++;
			}

			// this loop checks if the bead gets re-attaches to some other micelle core
			if (justEjected == 1)
			{
				if ((polymerBondStatus[i][j - 1].id1 < polymerBondStatus[i][j].id1) || (polymerBondStatus[i][j - 1].id2 < polymerBondStatus[i][j].id2))
				{
					if ((polymerBondStatus[i][j - 1].id1*polymerBondStatus[i][j].id1 == 0) || (polymerBondStatus[i][j - 1].id2*polymerBondStatus[i][j].id2 == 0))
					{
						justEjected = 0;
						
						// if the bead is getting attached to a different micelle core, then the counter is stored
						// if it re-attaches to the same micelle core, then the counter is rejected
						if (boundParticleID1 != polymerBondStatus[i][j].id1 || boundParticleID2 != polymerBondStatus[i][j].id2)
						{
							tauEm[currentIndex] = counter;
							counter = 0;
							currentIndex++;
						}
					}
				}
			}
		}
	}

	return tauEm;
}

int countEmtransitions (int nEm, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	nEm = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if ((polymerBondStatus[i][j - 1].id1 != polymerBondStatus[i][j].id1) || (polymerBondStatus[i][j - 1].id2 != polymerBondStatus[i][j].id2))
			{
				if (polymerBondStatus[i][j].id1 == 0 || polymerBondStatus[i][j].id2 == 0)
				{
					nEm++;
				}
			}
		}
	}

	return nEm;
}

// countTauS (tauS, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, nS);
float *countTauS (float *tauS, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, int nS)
{
	int currentIndex = 0, counter = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			if (polymerBondStatus[i][j].isItDangle == 1) {
				counter++; }

			if (polymerBondStatus[i][j].isItBridge == 1 && polymerBondStatus[i][j - 1].isItDangle == 1) {
					tauS[currentIndex] = counter;
					counter = 0;
					currentIndex++; }

			if (polymerBondStatus[i][j].isItLoop == 1 && polymerBondStatus[i][j - 1].isItDangle == 1) {
					counter = 0; }
		}
	}

	return tauS;
}

int countStransitions (int nS, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	nS = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (polymerBondStatus[i][j - 1].isItDangle == 1 && polymerBondStatus[i][j].isItBridge == 1)
			{
				nS++;
			}
		}
	}

	return nS;
}

float *countTauSm (float *tauSm, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, int nSm)
{
	int currentIndex = 0, counter = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			if (polymerBondStatus[i][j].isItDangle == 1) {
				counter++; }

			if (polymerBondStatus[i][j].isItLoop == 1 || polymerBondStatus[i][j].isItBridge == 1) {
				if (polymerBondStatus[i][j - 1].isItDangle == 1) {
					tauSm[currentIndex] = counter;
					counter = 0;
					currentIndex++; }
			}
		}
	}

	return tauSm;
}

int countSmtransitions (int nSm, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	nSm = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (polymerBondStatus[i][j - 1].isItDangle == 1 && polymerBondStatus[i][j].isItLoop == 1)
			{
				nSm++;
			}
			else if (polymerBondStatus[i][j - 1].isItDangle == 1 && polymerBondStatus[i][j].isItBridge == 1)
			{
				nSm++;
			}
		}
	}

	return nSm;
}

int countE_at_tansitions (int nE_at, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes, int nDumpEntries, BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, DATA_ATOMS *sortedAtoms)
{
	nE_at = 0;

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (beadBoundStatus[i][j - 1].isItBound == 1 && beadBoundStatus[i][j].isItBound == 0)
			{
				nE_at++;
			}
		}
	}

	return nE_at;
}

int countE_at_bridge_transitions (int nE_at_bridge, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes, BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, DATA_ATOMS *sortedAtoms, int nDumpEntries)
{
	nE_at_bridge = 0;
	int bondNumber;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	int atomNumber1, atomNumber2;

	int debugg = 0, previousCounter = 0;
	int currentlyBridge = 0;

	int particleID1 = 0, particleID2 = 0, counter = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentlyBridge = 0;
		counter = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			// If the bridge is stable
			// this check prevents bridge to bridge transitions at the beginning
			// after a stable bridge was found, the counter starts 
			if ((polymerBondStatus[i][j].isItBridge == 1) && 
				(polymerBondStatus[i][j - 1].id1 == polymerBondStatus[i][j].id1) && 
				(polymerBondStatus[i][j - 1].id2 == polymerBondStatus[i][j].id2) &&
				(counter == 0))
			{
				currentlyBridge = 1;

				particleID1 = polymerBondStatus[i][j].id1;
				particleID2 = polymerBondStatus[i][j].id2;
			}

			// Sometimes, a bead is still bound (with pair interaction energy of greater than -1 kBT, but less than 0 kBT).
			// During these situations, it is not considered as a bridge and thus, it will not have a bound particle ID.
			// The 'if' condition below is checking for that situation.
			// If this is the case, then the counter must continue and not break.
			// The counter only breaks if both the beads are 'unbound', or
			// it can break when both the beads are bound, but one of them have a different 
			// bound particle ID from the previous timestep and it must not be 'zero'.
			// EDIT: The first 'if' condition now checks for this condition as well.
			// else if (counter > 0 && 
			// 		((polymerBondStatus[i][j].id1 == 0) || (polymerBondStatus[i][j].id2 == 0)) &&
			// 		(beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) &&
			// 		(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1))
			// {
			// 	counter++;
			// }

			// If the dumbbell is currently in a bridge state,
			if (currentlyBridge == 1)
			{
				// The counter stops when one of the two beads eject or,
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 0) || 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 0))
				{
					nE_at_bridge++;

					if (debugg == 1)
					{
						printf(">Bridge ejection detected...<\n");
					}
				}
				// when both the beads stay bounded, but to a different particles, indicating bridge to bridge transition.
				else if (((particleID1 != polymerBondStatus[i][j].id1) || 
						(particleID2 != polymerBondStatus[i][j].id2)) &&
						(polymerBondStatus[i][j].id1 != 0) &&
						(polymerBondStatus[i][j].id2 != 0))
				{
					nE_at_bridge++;

					if (debugg == 1)
					{
						printf(">Bridge ejected and became another bridge very quickly...<\n");
					}
				}
				// when both of the beads are still bounded, but they are now a loop.
				else if ((polymerBondStatus[i][j].id1 == polymerBondStatus[i][j].id2) &&
						(polymerBondStatus[i][j].id1 != 0) &&
						(polymerBondStatus[i][j].id2 != 0))
				{
					nE_at_bridge++;

					if (debugg == 1)
					{
						printf(">Bridge ejected and became another loop very quickly...<\n");
					}
				}
			}

			if (debugg == 1)
			{
				printf("%d (%d %d) + %d (%d %d) [bridge: %d] => %d\n", 
					atomNumber1, 
					beadBoundStatus[atomNumber1 - 1][j].isItBound, 
					polymerBondStatus[i][j].id1, 
					atomNumber2, 
					beadBoundStatus[atomNumber2 - 1][j].isItBound, 
					polymerBondStatus[i][j].id2,
					polymerBondStatus[i][j].isItBridge, 
					nE_at_bridge);

				if (nE_at_bridge != previousCounter)
				{
					sleep (2);
				}

				previousCounter = nE_at_bridge;
			}
		}
	}

	// The following code uses a different approach for calculation, not used in the paper
	/*for (int i = 0; i < datafile.nAtoms; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (sortedAtoms[i].atomType == 1)
			{
				bondNumber = floor ((sortedAtoms[i].id - (int)floor (sortedAtoms[i].id / (coordinationNumber + 1)) - 1) / 2);

				// if (beadBoundStatus[sortedAtoms[i].id - 1][j - 1].isItBound == 1 && beadBoundStatus[sortedAtoms[i].id - 1][j].isItBound == 0 && polymerBondStatus[bondNumber][j - 1].isItBridge == 1 && polymerBondStatus[bondNumber][j].isItDangle == 1)
				// {
				// 	nE_at_bridge++;
				// }

				if (polymerBondStatus[bondNumber][j - 1].isItBridge == 1)
				{
					if ((polymerBondStatus[bondNumber][j - 1].id1 != polymerBondStatus[bondNumber][j].id1) || (polymerBondStatus[bondNumber][j - 1].id2 != polymerBondStatus[bondNumber][j].id2))
					{
						nE_at_bridge++;
					}
				}
			}
		}
	}*/

	return nE_at_bridge;
}

float *countTauE_at_bridge (float *tauE_at_bridge, int nE_at_bridge, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes, BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, DATA_ATOMS *sortedAtoms, int nDumpEntries)
{
	int counter = 0, currentIndex = 0;

	int bondNumber;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;

	int atomNumber1, atomNumber2;
	int particleID1 = 0, particleID2 = 0;

	int debugg = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			// Once the counter starts, the counter continues as long as both
			// the beads are bound to the same two particles
			if (counter > 0 && 
				(beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) &&
				(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
				((particleID1 == polymerBondStatus[i][j].id1 && particleID2 == polymerBondStatus[i][j].id2) ||
				(polymerBondStatus[i][j].id1 == 0 && particleID2 == polymerBondStatus[i][j].id2) ||
				(particleID1 == polymerBondStatus[i][j].id1 && polymerBondStatus[i][j].id2 == 0) ||
				(polymerBondStatus[i][j].id1 == 0 && polymerBondStatus[i][j].id2 == 0)))
			{
				counter++;
			}

			// After the counter starts,
			if (counter > 0)
			{
				// The counter stops when one of the two beads eject or,
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 0) || 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 0))
				{
					tauE_at_bridge[currentIndex] = counter;
					counter = 0;
					currentIndex++;

					if (debugg == 1)
					{
						printf(">counter ends: 1<\n");
					}
				}
				// when both the beads stay bounded, but to a different particles, indicating bridge to bridge transition.
				else if (((particleID1 != polymerBondStatus[i][j].id1) || 
						(particleID2 != polymerBondStatus[i][j].id2)) &&
						(polymerBondStatus[i][j].id1 != 0) &&
						(polymerBondStatus[i][j].id2 != 0))
				{
					tauE_at_bridge[currentIndex] = counter;
					counter = 0;
					currentIndex++;

					if (debugg == 1)
					{
						printf(">counter ends: 2<\n");
					}
				}
				// when both of the beads are still bounded, but they are now a loop.
				else if ((polymerBondStatus[i][j].id1 == polymerBondStatus[i][j].id2) &&
						(polymerBondStatus[i][j].id1 != 0) &&
						(polymerBondStatus[i][j].id2 != 0))
				{
					tauE_at_bridge[currentIndex] = counter;
					counter = 0;
					currentIndex++;

					if (debugg == 1)
					{
						printf(">counter ends: 3<\n");
					}
				}
			}

			// If the bridge is stable
			// this check prevents bridge to bridge transitions
			// after a stable bridge was found, the counter starts 
			if ((polymerBondStatus[i][j].isItBridge == 1) && 
				(polymerBondStatus[i][j - 1].id1 == polymerBondStatus[i][j].id1) && 
				(polymerBondStatus[i][j - 1].id2 == polymerBondStatus[i][j].id2) &&
				(counter == 0))
			{
				counter++;

				particleID1 = polymerBondStatus[i][j].id1;
				particleID2 = polymerBondStatus[i][j].id2;
			}

			// Sometimes, a bead is still bound (with pair interaction energy of greater than -1 kBT, but less than 0 kBT).
			// During these situations, it is not considered as a bridge and thus, it will not have a bound particle ID.
			// The 'if' condition below is checking for that situation.
			// If this is the case, then the counter must continue and not break.
			// The counter only breaks if both the beads are 'unbound', or
			// it can break when both the beads are bound, but one of them have a different 
			// bound particle ID from the previous timestep and it must not be 'zero'.
			// EDIT: The first 'if' condition now checks for this condition as well.
			// else if (counter > 0 && 
			// 		((polymerBondStatus[i][j].id1 == 0) || (polymerBondStatus[i][j].id2 == 0)) &&
			// 		(beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) &&
			// 		(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1))
			// {
			// 	counter++;
			// }

			if (debugg == 1)
			{
				printf("at1: %d (%d %d) + at2: %d (%d %d) [%d] => count: %d\n", 
					atomNumber1, 
					beadBoundStatus[atomNumber1 - 1][j].isItBound, 
					polymerBondStatus[i][j].id1, 
					atomNumber2, 
					beadBoundStatus[atomNumber2 - 1][j].isItBound, 
					polymerBondStatus[i][j].id2, 
					polymerBondStatus[i][j].isItBridge,
					counter);

				if (polymerBondStatus[i][j].isItBridge == 1 && (counter%100 == 0))
				{
					usleep (1000000);
				}
			}
		}
	}

	// The following code was not used in the paper
	// The counter works for every polymer beads individually.
	// Previous loop works for every polymer dumbbell, 
	// which is compatible with other functions, thus error free.

	/*for (int i = 0; i < datafile.nAtoms; ++i)
	{
		counter = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			if (sortedAtoms[i].atomType == 1)
			{
				bondNumber = floor ((sortedAtoms[i].id - (int)floor (sortedAtoms[i].id / (coordinationNumber + 1)) - 1) / 2);

				// if (beadBoundStatus[sortedAtoms[i].id - 1][j - 1].isItBound == 1 && beadBoundStatus[sortedAtoms[i].id - 1][j].isItBound == 0 && polymerBondStatus[bondNumber][j - 1].isItBridge == 1 && polymerBondStatus[bondNumber][j].isItDangle == 1)
				// {
				// 	nE_at_bridge++;
				// }

				if (polymerBondStatus[bondNumber][j].isItBridge == 1)
				{
					counter++;
				}

				if (polymerBondStatus[bondNumber][j - 1].isItBridge == 1)
				{
					if ((polymerBondStatus[bondNumber][j - 1].id1 != polymerBondStatus[bondNumber][j].id1) || (polymerBondStatus[bondNumber][j - 1].id2 != polymerBondStatus[bondNumber][j].id2))
					{
						tauE_at_bridge[currentIndex] = counter;
						counter = 0;
						currentIndex++;
					}
				}
			}
		}
	}*/

	return tauE_at_bridge;
}

int countE_at_dangle_transitions (int nE_at_dangle, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes, BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, DATA_ATOMS *sortedAtoms, int nDumpEntries)
{
	nE_at_dangle = 0;
	int bondNumber;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;

	// for (int i = 0; i < datafile.nBonds; ++i)
	// {
	// 	for (int j = 1; j < nTimeframes; ++j)
	// 	{
	// 		atomNumber1 = ((i + 1) * 2) - 1 + m.floor ((i * 2) / (coordinationNumber));
	// 		atomNumber2 = ((i + 1) * 2) + m.floor ((i * 2) / (coordinationNumber));

	// 		if ((beadBoundStatus[atomNumber1 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber1 - 1][j].isItBound == 0) || 
	// 			(beadBoundStatus[atomNumber2 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0))
	// 		{
	// 			if (polymerBondStatus[i][j - 1].isItBridge == 1)
	// 			{
	// 				nE_at_bridge++;
	// 			}
	// 		}
	// 		else if ((polymerBondStatus[i][j - 1].id1 != polymerBondStatus[i][j].id1) || 
	// 				(polymerBondStatus[i][j - 1].id2 != polymerBondStatus[i][j].id2))
	// 		{
	// 			nE_at_bridge++;
	// 		}
	// 	}
	// }

	// The following code was not used in the paper
	// The counter works for every polymer beads individually.
	// Previous loop works for every polymer dumbbell, 
	// which is compatible with other functions, thus error free.

	/*for (int i = 0; i < datafile.nAtoms; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (sortedAtoms[i].atomType == 1)
			{
				bondNumber = floor ((sortedAtoms[i].id - (int)floor (sortedAtoms[i].id / (coordinationNumber + 1)) - 1) / 2);

				// if (beadBoundStatus[i][j - 1].isItBound == 1 && beadBoundStatus[i][j].isItBound == 0 && polymerBondStatus[][j - 1].isItDangle == 1 && polymerBondStatus[][j].isItFree == 1)
				// {
				// 	nE_at_dangle++;
				// }

				if (polymerBondStatus[bondNumber][j - 1].isItDangle == 1 && polymerBondStatus[bondNumber][j].isItDangle == 1)
				{
					if ((polymerBondStatus[bondNumber][j - 1].id1 != polymerBondStatus[bondNumber][j].id1) || (polymerBondStatus[bondNumber][j - 1].id2 != polymerBondStatus[bondNumber][j].id2))
					{
						nE_at_dangle++;
					}
				}
				else if (polymerBondStatus[bondNumber][j - 1].isItDangle == 1 && polymerBondStatus[bondNumber][j].isItFree == 1)
				{
					nE_at_dangle++;
				}
			}
		}
	}*/

	return nE_at_dangle;
}

float *countTauE_at_dangle (float *tauE_at_dangle, int nE_at_dangle, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes, BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, DATA_ATOMS *sortedAtoms, int nDumpEntries)
{
	int counter = 0, currentIndex = 0;

	int bondNumber;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;

	// The following code was not used in the paper
	// The counter works for every polymer beads individually.
	// Previous loop works for every polymer dumbbell, 
	// which is compatible with other functions, thus error free.

	/*for (int i = 0; i < datafile.nAtoms; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (sortedAtoms[i].atomType == 1)
			{
				bondNumber = floor ((sortedAtoms[i].id - (int)floor (sortedAtoms[i].id / (coordinationNumber + 1)) - 1) / 2);

				// if (beadBoundStatus[i][j - 1].isItBound == 1 && beadBoundStatus[i][j].isItBound == 0 && polymerBondStatus[][j - 1].isItDangle == 1 && polymerBondStatus[][j].isItFree == 1)
				// {
				// 	nE_at_dangle++;
				// }

				if (polymerBondStatus[bondNumber][j].isItDangle == 1)
				{
					counter++;
				}

				if (polymerBondStatus[bondNumber][j - 1].isItDangle == 1 && polymerBondStatus[bondNumber][j].isItDangle == 1)
				{
					if ((polymerBondStatus[bondNumber][j - 1].id1 != polymerBondStatus[bondNumber][j].id1) || (polymerBondStatus[bondNumber][j - 1].id2 != polymerBondStatus[bondNumber][j].id2))
					{
						tauE_at_dangle[currentIndex] = counter;
						counter = 0;
						currentIndex++;
					}
				}
				else if (polymerBondStatus[bondNumber][j - 1].isItDangle == 1 && polymerBondStatus[bondNumber][j].isItFree == 1)
				{
					tauE_at_dangle[currentIndex] = counter;
					counter = 0;
					currentIndex++;
				}
			}
		}
	}*/

	return tauE_at_dangle;

}

int countE_at_loop_transitions (int nE_at_loop, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes, BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, DATA_ATOMS *sortedAtoms, int nDumpEntries)
{
	nE_at_loop = 0;

	int bondNumber;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;

	int atomNumber1, atomNumber2;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if ((beadBoundStatus[atomNumber1 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber1 - 1][j].isItBound == 0) || 
				(beadBoundStatus[atomNumber2 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0))
			{
				if (polymerBondStatus[i][j - 1].isItLoop == 1)
				{
					nE_at_loop++;
				}
			}
		}
	}

	// The following code was not used in the paper
	// The counter works for every polymer beads individually.
	// Previous loop works for every polymer dumbbell, 
	// which is compatible with other functions, thus error free.

	/*for (int i = 0; i < datafile.nAtoms; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (sortedAtoms[i].atomType == 1)
			{
				bondNumber = floor ((sortedAtoms[i].id - (int)floor (sortedAtoms[i].id / (coordinationNumber + 1)) - 1) / 2);

				// if (beadBoundStatus[i][j - 1].isItBound == 1 && beadBoundStatus[i][j].isItBound == 0 && polymerBondStatus[][j - 1].isItLoop == 1 && polymerBondStatus[][j].isItDangle == 1)
				// {
				// 	nE_at_loop++;
				// }

				if (polymerBondStatus[bondNumber][j - 1].isItLoop == 1 && (polymerBondStatus[bondNumber][j].isItDangle == 1 || polymerBondStatus[bondNumber][j].isItFree == 1))
				{
					nE_at_loop++;
				}
			}
		}
	}*/

	return nE_at_loop;
}

float *countTauE_at_loop (float *tauE_at_loop, int nE_at_loop, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes, BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, DATA_ATOMS *sortedAtoms, int nDumpEntries)
{
	int counter = 0, currentIndex = 0;

	int bondNumber;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	int atomNumber1, atomNumber2;

	int debugg = 0;
	int particleID1 = 0, particleID2 = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (counter > 0 && 
				(beadBoundStatus[atomNumber1 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) &&
				(beadBoundStatus[atomNumber2 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1))
			{
				counter++;
			}

			if ((polymerBondStatus[i][j].isItLoop == 1) && 
				(polymerBondStatus[i][j - 1].isItLoop == 1) &&
				(counter == 0))
			{
				counter++;
				particleID1 = polymerBondStatus[i][j].id1;
				particleID2 = polymerBondStatus[i][j].id2;
			}

			if (debugg == 1)
			{
				printf("at1: %d (%d %d); at2: %d (%d %d) => count: %d\n", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, counter);

				if (counter%100 == 0 && counter > 0)
				{
					sleep (1);
				}
			}

			if (counter > 0)
			{
				if ((beadBoundStatus[atomNumber1 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber1 - 1][j].isItBound == 0) || 
					(beadBoundStatus[atomNumber2 - 1][j - 1].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0))
				{
					tauE_at_loop[currentIndex] = counter;
					counter = 0;
					currentIndex++;

					if (debugg == 1)
					{
						printf(">counter ends 1<\n");
						sleep (1);
					}
				}
				else if (((particleID1 != polymerBondStatus[i][j].id1) || 
						(particleID2 != polymerBondStatus[i][j].id2)) &&
						(polymerBondStatus[i][j].id1 != 0) &&
						(polymerBondStatus[i][j].id2 != 0) &&
						(polymerBondStatus[i][j].isItLoop == 0))
				{
					tauE_at_loop[currentIndex] = counter;
					counter = 0;
					currentIndex++;

					if (debugg == 1)
					{
						printf(">counter ends 2<\n");
						sleep (1);
					}
				}
			}
		}
	}

	// The following code was not used in the paper
	// The counter works for every polymer beads individually.
	// Previous loop works for every polymer dumbbell, 
	// which is compatible with other functions, thus error free.

	/*for (int i = 0; i < datafile.nAtoms; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (sortedAtoms[i].atomType == 1)
			{
				bondNumber = floor ((sortedAtoms[i].id - (int)floor (sortedAtoms[i].id / (coordinationNumber + 1)) - 1) / 2);

				// if (beadBoundStatus[i][j - 1].isItBound == 1 && beadBoundStatus[i][j].isItBound == 0 && polymerBondStatus[][j - 1].isItLoop == 1 && polymerBondStatus[][j].isItDangle == 1)
				// {
				// 	nE_at_loop++;
				// }

				if (polymerBondStatus[bondNumber][j].isItLoop == 1)
				{
					counter++;
				}

				if (polymerBondStatus[bondNumber][j - 1].isItLoop == 1 && (polymerBondStatus[bondNumber][j].isItDangle == 1 || polymerBondStatus[bondNumber][j].isItFree == 1))
				{
					tauE_at_loop[currentIndex] = counter;
					counter = 0;
					currentIndex++;
				}
			}
		}
	}*/

	return tauE_at_loop;
}

int countS_at_transitions (int nS_at, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes, int nDumpEntries, BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0;
	nS_at = 0;

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (debugg == 1)
			{
				printf("%d\n", beadBoundStatus[i][j - 1].isItBound);
				usleep (100000);
			}

			if (beadBoundStatus[i][j - 1].isItBound == 0 && beadBoundStatus[i][j].isItBound == 1)
			{
				nS_at++;
			}
		}
	}

	return nS_at;
}

float *countTauE_at (float *tauE_at, int nE_at, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes)
{
	int counter = 0, currentIndex = 0;

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (beadBoundStatus[i][j].isItBound == 1)
			{
				counter++;
			}

			if (beadBoundStatus[i][j - 1].isItBound == 1 && beadBoundStatus[i][j].isItBound == 0 && counter > 0)
			{
				tauE_at[currentIndex] = counter;
				counter = 0;
				currentIndex++;
			}
		}
	}

	return tauE_at;
}

float *countTauS_at (float *tauS_at, int nS_at, BOUND_STATUS **beadBoundStatus, DATAFILE_INFO datafile, int nTimeframes)
{
	int counter = 0, currentIndex = 0;

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		for (int j = 1; j < nTimeframes; ++j)
		{
			if (beadBoundStatus[i][j].isItBound == 0)
			{
				counter++;
			}

			if (beadBoundStatus[i][j - 1].isItBound == 0 && beadBoundStatus[i][j].isItBound == 1 && counter > 0)
			{
				tauS_at[currentIndex] = counter;
				counter = 0;
				currentIndex++;
			}
		}
	}

	return tauS_at;
}

int count_nBDBs (int nBDBs, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromBridge = 0, currentlyBridge = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	nBDBs = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentlyBridge = 0;
		dangleFromBridge = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				fprintf(stdout, "\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, nBDBs);

				if (dangleFromBridge == 1 || currentlyBridge == 1)
				{
					// usleep (100000);
				}
			}

			if ((currentlyBridge == 1) && 
				((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyBridge = 0;
				dangleFromBridge = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from bridge<");
					fflush (stdout);
					sleep (3);
				}
			}

			if ((polymerBondStatus[i][j - 1].isItBridge == 1) && 
				(currentlyBridge == 0) && 
				(dangleFromBridge == 0))
			{
				currentlyBridge = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >bridge detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromBridge == 1)
			{
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 != polymerBondStatus[i][j].id2))
				{
					dangleFromBridge = 0;

					if ((particleID1 == polymerBondStatus[i][j].id1) && 
						(particleID2 == polymerBondStatus[i][j].id2))
					{
						nBDBs++;
					}

					if (debugg == 1)
					{
						fprintf(stdout, " >bridge forms from dangle<\n");
						fflush (stdout);
						sleep (3);
					}
				}
			}
		}
	}

	return nBDBs;
}

float *count_tau_BDBs (float *tau_BDBs, int nBDBs, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromBridge = 0, currentlyBridge = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	int counter = 0, currentIndex = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;
		currentlyBridge = 0;
		dangleFromBridge = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				printf("\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, counter);

				if (dangleFromBridge == 1 || currentlyBridge == 1)
				{
					// usleep (100000);
				}
			}

			if (currentlyBridge == 1 && ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyBridge = 0;
				dangleFromBridge = 1;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from bridge<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromBridge == 1)
			{
				counter++;
			}

			if (polymerBondStatus[i][j - 1].isItBridge == 1 && currentlyBridge == 0 && dangleFromBridge == 0)
			{
				currentlyBridge = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >bridge detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			// Possible transitions
			// dangle -> dangle
			// dangle -> bridge (same)
			// dangle -> bridge (different)
			// dangle -> loop
			// dangle -> free
			if (dangleFromBridge == 1)
			{
				// dangle -> bridge
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 != polymerBondStatus[i][j].id2))
				{
					// dangle -> bridge (same)
					if ((particleID1 == polymerBondStatus[i][j].id1) && 
						(particleID2 == polymerBondStatus[i][j].id2))
					{
						currentIndex++;
						tau_BDBs[currentIndex] = counter;
					}

					// counter is zeroed even for 'dangle -> bridge (different)'
					dangleFromBridge = 0;
					counter = 0;
					currentlyBridge = 0;

					if (debugg == 1)
					{
						fprintf(stdout, " >bridge forms from dangle<\n");
						fflush (stdout);
						sleep (3);
					}
				}
				// dangle -> loop
				else if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 == polymerBondStatus[i][j].id2))
				{
					counter = 0;
					dangleFromBridge = 0;
					currentlyBridge = 0;
				}
				// dangle -> free
				else if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 0) &&
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1))
				{
					counter = 0;
					dangleFromBridge = 0;
					currentlyBridge = 0;
				}
				// dangle -> dangle
				else
				{
					dangleFromBridge = 1;
				}
			}
		}
	}

	return tau_BDBs;
}

int count_nBDBd (int nBDBd, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromBridge = 0, currentlyBridge = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	nBDBd = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentlyBridge = 0;
		dangleFromBridge = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				printf("\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, nBDBd);

				if (dangleFromBridge == 1 || currentlyBridge == 1)
				{
					// usleep (100000);
				}
			}

			if ((currentlyBridge == 1) && 
				((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyBridge = 0;
				dangleFromBridge = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from bridge<");
					fflush (stdout);
					sleep (3);
				}
			}

			if ((polymerBondStatus[i][j - 1].isItBridge == 1) && 
				(currentlyBridge == 0) && 
				(dangleFromBridge == 0))
			{
				currentlyBridge = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >bridge detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromBridge == 1)
			{
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 != polymerBondStatus[i][j].id2))
				{
					dangleFromBridge = 0;

					if ((particleID1 != polymerBondStatus[i][j].id1) && 
						(particleID2 != polymerBondStatus[i][j].id2))
					{
						nBDBd++;

						if (debugg == 1)
						{
							fprintf(stdout, " >bridge forms from dangle<\n");
							fflush (stdout);
							sleep (3);
						}
					}
				}
			}
		}
	}

	return nBDBd;
}

float *count_tau_BDBd (float *tau_BDBd, int nBDBd, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromBridge = 0, currentlyBridge = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	int counter = 0, currentIndex = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;
		currentlyBridge = 0;
		dangleFromBridge = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				printf("\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, counter);

				if (dangleFromBridge == 1 || currentlyBridge == 1)
				{
					// usleep (100000);
				}
			}

			if (currentlyBridge == 1 && ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyBridge = 0;
				dangleFromBridge = 1;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from bridge<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromBridge == 1)
			{
				counter++;
			}

			if (polymerBondStatus[i][j - 1].isItBridge == 1 && currentlyBridge == 0 && dangleFromBridge == 0)
			{
				currentlyBridge = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >bridge detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			// POSSIBLE TRANSITIONS FROM A DANGLE:
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// dangle -> free
			// dangle -> loop
			// dangle -> bridge (same)
			// dangle -> bridge (different)

			if (dangleFromBridge == 1)
			{
				// if the dangle becomes a bridge
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 != polymerBondStatus[i][j].id2))
				{
					// if the dangle becomes a different bridge
					if ((particleID1 != polymerBondStatus[i][j].id1) && 
						(particleID2 != polymerBondStatus[i][j].id2))
					{
						currentIndex++;
						tau_BDBd[currentIndex] = counter;

						if (debugg == 1)
						{
							fprintf(stdout, " >bridge forms from dangle<\n");
							fflush (stdout);
							sleep (3);
						}
					}
					// if the dangle becomes a bridge with the same particles
					else
					{
						counter = 0;
					}

					dangleFromBridge = 0;
					counter = 0;
				}
				// if it stays as a dangle
				else if (((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || 
					(beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
				{
					dangleFromBridge = 1; // no change in status
				}
				// if it is not a dangle and not a bridge
				else
				{
					counter = 0;
				}
			}
		}
	}

	return tau_BDBd;
}

int count_nBDL (int nBDL, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromBridge = 0, currentlyBridge = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	nBDL = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentlyBridge = 0;
		dangleFromBridge = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				printf("\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, nBDL);

				if (dangleFromBridge == 1 || currentlyBridge == 1)
				{
					// usleep (100000);
				}
			}

			if ((currentlyBridge == 1) && 
				((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyBridge = 0;
				dangleFromBridge = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from bridge<");
					fflush (stdout);
					sleep (3);
				}
			}

			if ((polymerBondStatus[i][j - 1].isItBridge == 1) && 
				(currentlyBridge == 0) && 
				(dangleFromBridge == 0))
			{
				currentlyBridge = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >bridge detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromBridge == 1)
			{
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 == polymerBondStatus[i][j].id2))
				{
					dangleFromBridge = 0;
					nBDL++;

					if (debugg == 1)
					{
						fprintf(stdout, " >loop forms from dangle<\n");
						fflush (stdout);
						sleep (3);
					}
				}
			}
		}
	}

	return nBDL;
}

float *count_tau_BDL (float *tau_BDL, int nBDL, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromBridge = 0, currentlyBridge = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	int counter = 0, currentIndex = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;
		currentlyBridge = 0;
		dangleFromBridge = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				printf("\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, counter);

				if (dangleFromBridge == 1 || currentlyBridge == 1)
				{
					// usleep (100000);
				}
			}

			if (currentlyBridge == 1 && ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyBridge = 0;
				dangleFromBridge = 1;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from bridge<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromBridge == 1)
			{
				counter++;
			}

			if (polymerBondStatus[i][j - 1].isItBridge == 1 && currentlyBridge == 0 && dangleFromBridge == 0)
			{
				currentlyBridge = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >bridge detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			// POSSIBLE TRANSITIONS FROM A DANGLE:
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// dangle -> free
			// dangle -> loop
			// dangle -> bridge (same)
			// dangle -> bridge (different)

			if (dangleFromBridge == 1)
			{
				// to a loop
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 == polymerBondStatus[i][j].id2) &&
					(polymerBondStatus[i][j].id1 != 0) &&
					(polymerBondStatus[i][j].id2 != 0))
				{
					currentIndex++;
					tau_BDL[currentIndex] = counter;

					dangleFromBridge = 0;
					counter = 0;

					if (debugg == 1)
					{
						fprintf(stdout, " >loop forms from dangle<\n");
						fflush (stdout);
						sleep (3);
					}
				}

				// to a dangle
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || 
					(beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1))
				{
					dangleFromBridge = 1;
				}
				// if it is not a dangle and
				// not a loop (counter was not zeroed if the dumbbell did not form a loop)
				else if (counter > 0)
				{
					counter = 0;
					dangleFromBridge = 0;
					currentlyBridge = 0;
				}
			}
		}
	}

	return tau_BDL;
}

int count_nLDL (int nLDL, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromLoop = 0, currentlyLoop = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	nLDL = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentlyLoop = 0;
		dangleFromLoop = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				printf("\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, nLDL);

				if (dangleFromLoop == 1 || currentlyLoop == 1)
				{
					// usleep (100000);
				}
			}

			if ((currentlyLoop == 1) && 
				((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyLoop = 0;
				dangleFromLoop = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from loop<");
					fflush (stdout);
					sleep (3);
				}
			}

			if ((polymerBondStatus[i][j - 1].isItLoop == 1) && 
				(currentlyLoop == 0) && 
				(dangleFromLoop == 0))
			{
				currentlyLoop = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >loop detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromLoop == 1)
			{
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 == polymerBondStatus[i][j].id2))
				{
					dangleFromLoop = 0;
					nLDL++;

					if (debugg == 1)
					{
						fprintf(stdout, " >loop forms from dangle<\n");
						fflush (stdout);
						sleep (3);
					}
				}
			}
		}
	}

	return nLDL;
}

float *count_tau_LDL (float *tau_LDL, int nLDL, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)

{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromLoop = 0, currentlyLoop = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	int counter = 0, currentIndex = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;
		dangleFromLoop = 0;
		currentlyLoop = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				printf("\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, counter);

				if (dangleFromLoop == 1 || currentlyLoop == 1)
				{
					// usleep (100000);
				}
			}

			if (currentlyLoop == 1 && ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyLoop = 0;
				dangleFromLoop = 1;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from loop<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromLoop == 1)
			{
				counter++;
			}

			if (polymerBondStatus[i][j - 1].isItLoop == 1 && currentlyLoop == 0 && dangleFromLoop == 0)
			{
				currentlyLoop = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >loop detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			// Possible transitions
			// dangle -> loop
			// dangle -> bridge
			// dangle -> dangle
			// dangle -> free
			if (dangleFromLoop == 1)
			{
				// dangle -> loop
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 == polymerBondStatus[i][j].id2))
				{
					currentIndex++;
					tau_LDL[currentIndex] = counter;

					dangleFromLoop = 0;
					counter = 0;
					currentlyLoop = 0;

					if (debugg == 1)
					{
						fprintf(stdout, " >loop forms from dangle<\n");
						fflush (stdout);
						sleep (3);
					}
				}
				// dangle -> bridge
				else if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) &&
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 != polymerBondStatus[i][j].id2))
				{
					counter = 0;
					dangleFromLoop = 0;
					currentlyLoop = 0;
				}
				// dangle -> free
				else if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 0) &&
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 0))
				{
					counter = 0;
					dangleFromLoop = 0;
					currentlyLoop = 0;
				}
				// dangle -> dangle
				else
				{
					dangleFromLoop = 1;
				}
			}
		}
	}

	return tau_LDL;
}

int count_nLDB (int nLDB, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromLoop = 0, currentlyLoop = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	nLDB = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentlyLoop = 0;
		dangleFromLoop = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				printf("\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, nLDB);

				if (dangleFromLoop == 1 || currentlyLoop == 1)
				{
					// usleep (100000);
				}
			}

			if ((currentlyLoop == 1) && 
				((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyLoop = 0;
				dangleFromLoop = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from loop<");
					fflush (stdout);
					sleep (3);
				}
			}

			if ((polymerBondStatus[i][j - 1].isItLoop == 1) && 
				(currentlyLoop == 0) && 
				(dangleFromLoop == 0))
			{
				currentlyLoop = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >loop detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromLoop == 1)
			{
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 != polymerBondStatus[i][j].id2))
				{
					dangleFromLoop = 0;
					nLDB++;

					if (debugg == 1)
					{
						fprintf(stdout, " >bridge forms from dangle<\n");
						fflush (stdout);
						sleep (3);
					}
				}
			}
		}
	}

	return nLDB;
}

float *count_tau_LDB (float *tau_LDB, int nLDB, BOUND_STATUS **beadBoundStatus, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile, DATA_ATOMS *sortedAtoms)
{
	int debugg = 0, particleID1 = 0, particleID2 = 0, dangleFromLoop = 0, currentlyLoop = 0, atomNumber1 = 0, atomNumber2 = 0;
	int nPolymers = datafile.nBonds, nParticles = datafile.nAtoms - (datafile.nBonds * 2);
	int coordinationNumber = (nPolymers * 2) / nParticles;
	int counter = 0, currentIndex = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		counter = 0;
		dangleFromLoop = 0;
		currentlyLoop = 0;

		for (int j = 1; j < nTimeframes; ++j)
		{
			atomNumber1 = ((i + 1) * 2) - 1 + floor ((i * 2) / (coordinationNumber));
			atomNumber2 = ((i + 1) * 2) + floor ((i * 2) / (coordinationNumber));

			if (debugg == 1)
			{
				printf("\n%d (%d %d) + %d (%d %d) => %d", atomNumber1, beadBoundStatus[atomNumber1 - 1][j].isItBound, polymerBondStatus[i][j].id1, atomNumber2, beadBoundStatus[atomNumber2 - 1][j].isItBound, polymerBondStatus[i][j].id2, counter);

				if (dangleFromLoop == 1 || currentlyLoop == 1)
				{
					// usleep (100000);
				}
			}

			if (currentlyLoop == 1 && ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || (beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1)))
			{
				currentlyLoop = 0;
				dangleFromLoop = 1;

				if (debugg == 1)
				{
					fprintf(stdout, " >dangle formed from loop<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromLoop == 1)
			{
				counter++;
			}

			if (polymerBondStatus[i][j - 1].isItLoop == 1 && currentlyLoop == 0 && dangleFromLoop == 0)
			{
				currentlyLoop = 1;
				particleID1 = polymerBondStatus[i][j - 1].id1;
				particleID2 = polymerBondStatus[i][j - 1].id2;

				if (debugg == 1)
				{
					fprintf(stdout, " >loop detected<");
					fflush (stdout);
					sleep (3);
				}
			}

			if (dangleFromLoop == 1)
			{
				// to a bridge
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1) && 
					(beadBoundStatus[atomNumber2 - 1][j].isItBound == 1) &&
					(polymerBondStatus[i][j].id1 != polymerBondStatus[i][j].id2))
				{
					tau_LDB[currentIndex] = counter;
					currentIndex++;

					dangleFromLoop = 0;
					counter = 0;

					if (debugg == 1)
					{
						fprintf(stdout, " >bridge forms from dangle<\n");
						fflush (stdout);
						sleep (3);
					}
				}

				// to dangle
				if ((beadBoundStatus[atomNumber1 - 1][j].isItBound == 1 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 0) || 
					(beadBoundStatus[atomNumber1 - 1][j].isItBound == 0 && beadBoundStatus[atomNumber2 - 1][j].isItBound == 1))
				{
					dangleFromLoop = 1;
				}
				// if it is not a dangle and
				// not a bridge (counter was not zeroed if the dumbbell did not form a bridge)
				// So, the dumbbell must be in a loop or free state. Both requires us to zero the counter
				else if (counter > 0)
				{
					counter = 0;
					dangleFromLoop = 0;
					currentlyLoop = 0;
				}
			}
		}
	}

	return tau_LDB;
}

void printTauE_at (float *tauE_at, int nE_at, const char *folderName)
{
	char *tauE_at_file_filename;
	tauE_at_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauE_at_file_filename, 500, "%s/tauE_at.output", folderName);
	fprintf(stdout, "Printing %s\n", tauE_at_file_filename);

	FILE *tauE_at_file;
	tauE_at_file = fopen (tauE_at_file_filename, "w");

	for (int i = 0; i < nE_at; ++i)
	{
		fprintf(tauE_at_file, "%f\n", tauE_at[i]);
	}

	fclose (tauE_at_file);
}

void printTauE_at_bridge (float *tauE_at_bridge, int nE_at_bridge, const char *folderName)
{
	char *tauE_at_file_filename;
	tauE_at_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauE_at_file_filename, 500, "%s/tauE_at_bridge.output", folderName);

	FILE *tauE_at_bridge_file;
	tauE_at_bridge_file = fopen (tauE_at_file_filename, "w");

	for (int i = 0; i < nE_at_bridge; ++i)
	{
		fprintf(tauE_at_bridge_file, "%f\n", tauE_at_bridge[i]);
	}

	fclose (tauE_at_bridge_file);
}

void printTauE_at_loop (float *tauE_at_loop, int nE_at_loop, const char *folderName)
{
	char *tauE_at_file_filename;
	tauE_at_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauE_at_file_filename, 500, "%s/tauE_at_loop.output", folderName);

	FILE *tauE_at_loop_file;
	tauE_at_loop_file = fopen (tauE_at_file_filename, "w");

	for (int i = 0; i < nE_at_loop; ++i)
	{
		fprintf(tauE_at_loop_file, "%f\n", tauE_at_loop[i]);
	}

	fclose (tauE_at_loop_file);
}

void printTauE_at_dangle (float *tauE_at_dangle, int nE_at_dangle, const char *folderName)
{
	char *tauE_at_file_filename;
	tauE_at_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauE_at_file_filename, 500, "%s/tauE_at_dangle.output", folderName);

	FILE *tauE_at_dangle_file;
	tauE_at_dangle_file = fopen (tauE_at_file_filename, "w");

	for (int i = 0; i < nE_at_dangle; ++i)
	{
		fprintf(tauE_at_dangle_file, "%f\n", tauE_at_dangle[i]);
	}

	fclose (tauE_at_dangle_file);
}

void printTauS_at (float *tauS_at, int nS_at, const char *folderName)
{
	char *tauS_at_file_filename;
	tauS_at_file_filename = (char *) malloc (500 * sizeof (char));
	snprintf (tauS_at_file_filename, 500, "%s/tauS_at.output", folderName);

	FILE *tauS_at_file;
	tauS_at_file = fopen (tauS_at_file_filename, "w");

	for (int i = 0; i < nS_at; ++i)
	{
		fprintf(tauS_at_file, "%f\n", tauS_at[i]);
	}

	fclose (tauS_at_file);
}

void printTauBDBs (float *tau_BDBs, int nBDBs, const char *folderName)
{
	char *filename;
	filename = (char *) malloc (1000 * sizeof (char));
	snprintf (filename, 1000, "%s/tau_BDBs.output", folderName);
	fprintf(stdout, "Printing %s\n", filename);

	FILE *tau_file_output;
	tau_file_output = fopen (filename, "w");

	for (int i = 0; i < nBDBs; ++i)
	{
		fprintf(tau_file_output, "%f\n", tau_BDBs[i]);
	}

	fclose (tau_file_output);
}

void printTauBDBd (float *tau_BDBd, int nBDBd, const char *folderName)
{
	char *filename;
	filename = (char *) malloc (1000 * sizeof (char));
	snprintf (filename, 1000, "%s/tau_BDBd.output", folderName);
	fprintf(stdout, "Printing %s\n", filename);

	FILE *tau_file_output;
	tau_file_output = fopen (filename, "w");

	for (int i = 0; i < nBDBd; ++i)
	{
		fprintf(tau_file_output, "%f\n", tau_BDBd[i]);
	}

	fclose (tau_file_output);
}

void printTauBDL (float *tau_BDL, int nBDL, const char *folderName)
{
	char *filename;
	filename = (char *) malloc (1000 * sizeof (char));
	snprintf (filename, 1000, "%s/tau_BDL.output", folderName);
	fprintf(stdout, "Printing %s\n", filename);

	FILE *tau_file_output;
	tau_file_output = fopen (filename, "w");

	for (int i = 0; i < nBDL; ++i)
	{
		fprintf(tau_file_output, "%f\n", tau_BDL[i]);
	}

	fclose (tau_file_output);
}

void printTauLDL (float *tau_LDL, int nLDL, const char *folderName)
{
	char *filename;
	filename = (char *) malloc (1000 * sizeof (char));
	snprintf (filename, 1000, "%s/tau_LDL.output", folderName);
	fprintf(stdout, "Printing %s\n", filename);

	FILE *tau_file_output;
	tau_file_output = fopen (filename, "w");

	for (int i = 0; i < nLDL; ++i)
	{
		fprintf(tau_file_output, "%f\n", tau_LDL[i]);
	}

	fclose (tau_file_output);
}

void printTauLDB (float *tau_LDB, int nLDB, const char *folderName)
{
	char *filename;
	filename = (char *) malloc (1000 * sizeof (char));
	snprintf (filename, 1000, "%s/tau_LDB.output", folderName);
	fprintf(stdout, "Printing %s\n", filename);

	FILE *tau_file_output;
	tau_file_output = fopen (filename, "w");

	for (int i = 0; i < nLDB; ++i)
	{
		fprintf(tau_file_output, "%f\n", tau_LDB[i]);
	}

	fclose (tau_file_output);
}

int main(int argc, char const *argv[])
{
	if (argc != 7)
	{
		printf("ERROR: INCORRECT ARGUMENTS PASSED.\n\n Required arguments:\n\n{~} argv[0] = ./program\n{~} argv[1] = input dump file (ascii text or *.gz)\n{~} argv[2] = input data file.\n{~} argv[3] = dt for calculations.\n{~} argv[4] = Column number for energy entries\n{~} argv[5] = Number of timeframes to consider\n{~} argv[6] = Bridge to bridge correction (0 or 1)\n\n");
		exit (1);
	}

	FILE *inputDump, *outputStates;
	char *pipeString, *folderName;
	pipeString = (char *) malloc (500 * sizeof (char));
	folderName = (char *) malloc (500 * sizeof (char));
	snprintf (pipeString, 500, "gzcat %s", argv[1]);
	snprintf (folderName, 500, "stats_%s_%s_%d_%d_%d_%d", argv[1], argv[2], atoi (argv[3]), atoi (argv[4]), atoi (argv[5]), atoi (argv[6]));

	int N_TIMEFRAMES_TO_CONSIDER = atoi (argv[5]);
	int bbcorrection = atoi (argv[6]);

	fprintf(stdout, "N_TIMEFRAMES_TO_CONSIDER: %d\ndt: %d\nbbcorrection: %d\n\n", N_TIMEFRAMES_TO_CONSIDER, atoi (argv[3]), bbcorrection);
	fflush (stdout);

	int energyColumn = atoi (argv[4]);

	char *createFolder_command;
	createFolder_command = (char *) malloc (500 * sizeof (char));
	snprintf (createFolder_command, 500, "mkdir %s", folderName);
	system (createFolder_command);

	if (strstr (argv[1], ".gz"))
	{
		inputDump = popen (pipeString, "r");
	}
	else
	{
		inputDump = fopen (argv[1], "r");
	}

	int dt = atoi (argv[3]);

	char *outputStates_filename;
	outputStates_filename = (char *) malloc (500 * sizeof (char));
	snprintf (outputStates_filename, 500, "%s/polymerStates.timeseries", folderName);
	outputStates = fopen (outputStates_filename, "w");

/*	int nChains = NPARTICLES * NPOLYMERS, *particleIDs, nAtoms = NPARTICLES + (NPARTICLES * NPOLYMERS * NBEADS);
	particleIDs = (int *) malloc (NPARTICLES * sizeof (int));
	particleIDs = storeParticleIDs (particleIDs);
*/

	DATA_ATOMS *dataAtoms, *sortedAtoms;
	DATA_BONDS *dataBonds, *sortedBonds;
	DATA_ANGLES *dataAngles;
	DATA_DIHEDRALS *dataDihedrals;
	DATA_IMPROPERS *dataImpropers;
	DATAFILE_INFO datafile;

	datafile = readData (argv[2], &dataAtoms, &dataBonds, &dataAngles, &dataDihedrals, &dataImpropers);

	sortedBonds = (DATA_BONDS *) malloc (datafile.nBonds * sizeof (DATA_BONDS));
	sortedBonds = sortBonds (sortedBonds, dataBonds, datafile.nBonds, datafile.nAtoms);

	sortedAtoms = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));
	sortedAtoms = sortAtoms (sortedAtoms, dataAtoms, datafile.nAtoms);

	int nTimeframes = 0; //countNTimeframes (argv[1]);
	int N_TIMEFRAMES_TO_CONSIDER2;
	//N_TIMEFRAMES_TO_CONSIDER2 = nTimeframes;
	N_TIMEFRAMES_TO_CONSIDER2 = N_TIMEFRAMES_TO_CONSIDER;

	// if (N_TIMEFRAMES_TO_CONSIDER <= nTimeframes)
	// {
	// 	N_TIMEFRAMES_TO_CONSIDER2 = nTimeframes;
	// }
	// else
	// {
	// 	N_TIMEFRAMES_TO_CONSIDER2 = N_TIMEFRAMES_TO_CONSIDER;
	// }

	BOND_STATUS **polymerBondStatus;
	polymerBondStatus = (BOND_STATUS **) malloc (datafile.nBonds * sizeof (BOND_STATUS *));

	BOUND_STATUS **beadBoundStatus;
	beadBoundStatus = (BOUND_STATUS **) malloc (datafile.nAtoms * sizeof (BOUND_STATUS *));

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		polymerBondStatus[i] = (BOND_STATUS *) malloc ((N_TIMEFRAMES_TO_CONSIDER2 + 1) * sizeof (BOND_STATUS));
	}

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		beadBoundStatus[i]= (BOUND_STATUS *) malloc ((N_TIMEFRAMES_TO_CONSIDER2 + 1) * sizeof (BOUND_STATUS));
	}

	polymerBondStatus = initBondStatus (polymerBondStatus, datafile.nBonds, N_TIMEFRAMES_TO_CONSIDER2);

	int file_status = fgetc (inputDump);
	char lineString[3000];

	DUMP_ENERGY *energyEntries;
	int nDumpEntries = 0;
	nDumpEntries = countDumpEntries (inputDump, nDumpEntries);

	if (strstr (argv[1], ".gz"))
	{
		pclose (inputDump);
		FILE *inputDump;
		inputDump = popen (pipeString, "r");
	}
	else
	{
		rewind (inputDump);
	}

	energyEntries = (DUMP_ENERGY *) malloc (nDumpEntries * 2 * sizeof (DUMP_ENERGY));
	printf("Allocating %d mem (type DUMP_ENERGY) for variable energyEntries\n", nDumpEntries * 2);
	int currentTimeframe = 0;

	STATES polymerStates;

	int timeframesToSkip = 0;

	printf("Skipping %d timeframes initially...\n", timeframesToSkip);
	printf("To consider: %d\n", N_TIMEFRAMES_TO_CONSIDER2);

	timeframesToSkip = 0;
	int status = 0;

	int effectiveCurrentTimeframe = 0;

	while (file_status > 0)
	{
		printf("Scanning timeframe: %d (%d) / %d                  \r", currentTimeframe + 1, file_status, N_TIMEFRAMES_TO_CONSIDER2);
		fflush (stdout);

		energyEntries = initEnergyEntries (energyEntries, nDumpEntries);
		energyEntries = saveDumpEnergyEntries (inputDump, energyEntries, nDumpEntries, energyColumn, &status);

		if (status == 0)
		{
			goto leaveThisLoop;
		}

		if ((currentTimeframe % dt) == 0)
		{
			polymerBondStatus = checkBondStatus (polymerBondStatus, energyEntries, nDumpEntries, datafile, N_TIMEFRAMES_TO_CONSIDER2, effectiveCurrentTimeframe, sortedAtoms);

			beadBoundStatus = findBoundStates (beadBoundStatus, energyEntries, nDumpEntries, datafile, N_TIMEFRAMES_TO_CONSIDER2, effectiveCurrentTimeframe, sortedAtoms);

			polymerStates = countStates (polymerStates, polymerBondStatus, datafile, effectiveCurrentTimeframe);

			printStates (polymerStates, outputStates, datafile);

			effectiveCurrentTimeframe += 1;
		}

		if (currentTimeframe > N_TIMEFRAMES_TO_CONSIDER2) {
			goto leaveThisLoop; }

		file_status = fgetc (inputDump);
		currentTimeframe++;
	}

	leaveThisLoop:;

	free (energyEntries);

	if (strstr (argv[1], ".gz")) {
		pclose (inputDump); }
	else {
		fclose (inputDump); }

	fclose (outputStates);

	int nBB = countBBtransitions (nBB, polymerBondStatus, beadBoundStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile);
	fprintf(stdout, "nBB: %d\n", nBB);
	fflush (stdout);
	float *tauBB;

	tauBB = (float *) malloc (nBB * sizeof (float));
	tauBB = countTauBB (tauBB, polymerBondStatus, beadBoundStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, nBB);
	printTauBB (tauBB, nBB, folderName);

	int nE = countEtransitions (nE, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile), nEm = countEmtransitions (nEm, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile), nS = countStransitions (nS, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile), nSm = countSmtransitions (nSm, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile);

	float *tauE, *tauEm, *tauS, *tauSm;

	tauE = (float *) malloc (nE * sizeof (float));
	tauE = countTauE (tauE, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, nE);
	printTauE (tauE, nE, folderName);

	tauEm = (float *) malloc (nEm * sizeof (float));
	tauEm = countTauEm (tauEm, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, nEm);
	printTauEm (tauEm, nEm, folderName);

	tauS = (float *) malloc (nS * sizeof (float));
	tauS = countTauS (tauS, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, nS);
	printTauS (tauS, nS, folderName);

	tauSm = (float *) malloc (nSm * sizeof (float));
	tauSm = countTauSm (tauSm, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, nSm);
	printTauSm (tauSm, nSm, folderName);

	// Calculating tauE and tauS using a third route.
	// Here, the ejection is calculated as mentioned in Alyssa Travitz' paper. 
	// If the energy becomes less than -1, then it is considered as bound,
	// but the bead is considered unbound only if the energy is greater than 0.
	// This is another way of preventing quick transitions because BD simulations can lead
	// to very quick movements, which can be mis-read as very short ejection time.

	int nE_at, nS_at, nE_at_bridge, nE_at_loop, nE_at_dangle;
	nE_at = countE_at_tansitions (nE_at, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2, nDumpEntries, polymerBondStatus, energyEntries, sortedAtoms);
	nS_at = countS_at_transitions (nS_at, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2, nDumpEntries, polymerBondStatus, energyEntries, sortedAtoms);
	nE_at_bridge = countE_at_bridge_transitions (nE_at_bridge, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2, polymerBondStatus, energyEntries, sortedAtoms, nDumpEntries);
	nE_at_loop = countE_at_loop_transitions (nE_at_loop, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2, polymerBondStatus, energyEntries, sortedAtoms, nDumpEntries);
	nE_at_dangle = countE_at_dangle_transitions (nE_at_dangle, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2, polymerBondStatus, energyEntries, sortedAtoms, nDumpEntries);

	printf("Number of nE_at: %d; nS_at: %d\n", nE_at, nS_at);
	printf("nE_at_bridge: %d\nnE_at_loop: %d\nnE_at_dangle: %d\n", nE_at_bridge, nE_at_loop, nE_at_dangle);

	float *tauE_at, *tauS_at, *tauE_at_bridge, *tauE_at_loop, *tauE_at_dangle;
	tauE_at = (float *) malloc (nE_at * sizeof (float));
	tauE_at = countTauE_at (tauE_at, nE_at, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2);

	tauE_at_bridge = (float *) malloc (nE_at_bridge * sizeof (float));
	tauE_at_bridge = countTauE_at_bridge (tauE_at_bridge, nE_at_bridge, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2, polymerBondStatus, energyEntries, sortedAtoms, nDumpEntries);

	tauE_at_loop = (float *) malloc (nE_at_loop * sizeof (float));
	tauE_at_loop = countTauE_at_loop (tauE_at_loop, nE_at_loop, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2, polymerBondStatus, energyEntries, sortedAtoms, nDumpEntries);

	tauE_at_dangle = (float *) malloc (nE_at_dangle * sizeof (float));
	tauE_at_dangle = countTauE_at_dangle (tauE_at_dangle, nE_at_dangle, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2, polymerBondStatus, energyEntries, sortedAtoms, nDumpEntries);

	// Here, I am calculating the following transitions,
	// bridge -> dangle -> bridge (same)
	// bridge -> dangle -> bridge (different)
	// bridge -> dangle -> loop
	// loop -> dangle -> loop
	// loop -> dangle -> bridge
	int nBDBs = 0, nBDBd = 0, nBDL = 0, nLDL = 0, nLDB = 0;
	nBDBs = count_nBDBs (nBDBs, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);
	nBDBd = count_nBDBd (nBDBd, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);
	nBDL = count_nBDL (nBDL, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);
	nLDL = count_nLDL (nLDL, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);
	nLDB = count_nLDB (nLDB, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);

	fprintf(stdout, "nBDBs: %d\nnBDBd: %d\nnBDL: %d\nnLDL: %d\nnLDB: %d\n", nBDBs, nBDBd, nBDL, nLDL, nLDB);
	fflush (stdout);

	float *tau_BDBs, *tau_BDBd, *tau_BDL, *tau_LDL, *tau_LDB;
	tau_BDBs = (float *) malloc (nBDBs * sizeof (float));
	tau_BDBd = (float *) malloc (nBDBd * sizeof (float));
	tau_BDL = (float *) malloc (nBDL * sizeof (float));
	tau_LDL = (float *) malloc (nLDL * sizeof (float));
	tau_LDB = (float *) malloc (nLDB * sizeof (float));

	tau_BDBs = count_tau_BDBs (tau_BDBs, nBDBs, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);
	tau_BDBd = count_tau_BDBd (tau_BDBd, nBDBd, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);
	tau_BDL = count_tau_BDL (tau_BDL, nBDL, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);
	tau_LDL = count_tau_LDL (tau_LDL, nLDL, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);
	tau_LDB = count_tau_LDB (tau_LDB, nLDB, beadBoundStatus, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, sortedAtoms);

	int nLL;
	nLL = countLLtransitions (nLL, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, beadBoundStatus, sortedAtoms);

	fprintf(stdout, "nLL: %d\n", nLL);
	fflush (stdout);

	float *tauLL;
	tauLL = (float *) malloc (nLL * sizeof (float));
	tauLL = countTauLL (tauLL, nLL, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, beadBoundStatus, sortedAtoms);
	printTauLL (tauLL, nLL, folderName);

	printTauE_at (tauE_at, nE_at, folderName);
	printTauE_at_bridge (tauE_at_bridge, nE_at_bridge, folderName);
	printTauE_at_loop (tauE_at_loop, nE_at_loop, folderName);
	printTauE_at_dangle (tauE_at_dangle, nE_at_dangle, folderName);

	printTauBDBs (tau_BDBs, nBDBs, folderName);
	printTauBDBd (tau_BDBd, nBDBd, folderName);
	printTauBDL (tau_BDL, nBDL, folderName);
	printTauLDL (tau_LDL, nLDL, folderName);
	printTauLDB (tau_LDB, nLDB, folderName);
	fflush (stdout);

	tauS_at = (float *) malloc (nS_at * sizeof (float));
	tauS_at = countTauS_at (tauS_at, nS_at, beadBoundStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2);
	printTauS_at (tauS_at, nS_at, folderName);

	// computing transitions
	polymerBondStatus = correctingDangles (polymerBondStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2);

	float *tauBL, *tauBL2, *tauLB;
	int nBL = countBLtransitions (nBL, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile), nLB = countLBtransitions (nLB, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile);

	tauBL = (float *) malloc (nBL * sizeof (float));
	tauBL2 = (float *) malloc (nBL * sizeof (float));

	tauBL = initTauBL (nBL, tauBL);
	tauBL2 = initTauBL (nBL, tauBL2);

	tauBL = countTauBL (tauBL, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile, bbcorrection);
	tauBL2 = countTauBL2 (tauBL2, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile);

	printTauBL (tauBL, nBL, folderName);
	printTauBL2 (tauBL2, nBL, folderName);
	
	tauLB = (float *) malloc (nLB * sizeof (float));
	tauLB = countTauLB (tauLB, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile);
	printTauLB (tauLB, nLB, folderName);

	float averageLB, stdevLB, stderrLB;
	float averageBL, stdevBL, stderrBL;
	float averageBB, stdevBB, stderrBB;

	computeStats (&averageLB, &stdevLB, &stderrLB, tauLB, nLB);
	computeStats (&averageBL, &stdevBL, &stderrBL, tauBL, nBL);
	computeStats (&averageBB, &stdevBB, &stderrBB, tauBB, nBB);

	char *transitionStats_filename;
	transitionStats_filename = (char *) malloc (500 * sizeof (char));
	snprintf (transitionStats_filename, 500, "%s/transitions.stats", folderName);

	FILE *transitionStats;
	transitionStats = fopen (transitionStats_filename, "w");
	fprintf(transitionStats, "Bridge to loop:\n\naverage: %f\nstdev: %f\nstderr: %f\n\nLoop to bridge:\n\naverage: %f\nstdev: %f\nstderr: %f\n", averageBL, stdevBL, stderrBL, averageLB, stdevLB, stderrLB);
	fclose (transitionStats);

	// block averaging of transitions
	// BLOCKS *blockAverages_tauBL, *blockAverages_tauLB;
	// blockAverages_tauBL = (BLOCKS *) malloc (nBL * sizeof (BLOCKS));
	// blockAverages_tauLB = (BLOCKS *) malloc (nLB * sizeof (BLOCKS));

	// blockAverages_tauBL = initializeBlocks (blockAverages_tauBL, nBL);
	// blockAverages_tauLB = initializeBlocks (blockAverages_tauLB, nLB);

	// printf("\nComputing block averages on transitions...\n");
	// blockAverages_tauBL = computeBlockAverages (blockAverages_tauBL, nBL, tauBL);
	// blockAverages_tauLB = computeBlockAverages (blockAverages_tauLB, nLB, tauLB);

	// printing block averages of transitions
	// FILE *file_block_BL, *file_block_LB;
	// file_block_BL = fopen ("bridgeToLoop.block", "w");
	// file_block_LB = fopen ("loopToBridge.block", "w");

	// printf("Printing block averages of transitions...\n");
	// printBlockAverageStats (file_block_BL, blockAverages_tauBL, nBL);
	// printBlockAverageStats (file_block_LB, blockAverages_tauLB, nLB);

	// fclose (file_block_BL);
	// fclose (file_block_LB);
	// free (blockAverages_tauBL);
	// free (blockAverages_tauLB);
	free (tauBL);
	free (tauLB);
	free (polymerBondStatus);

	// block averaging of microstates
	FILE *file_data, *file_block_nBridges, *file_block_nLoops, *file_block_nDangles, *file_block_nFree;

	char *file_data_filename;
	file_data_filename = (char *) malloc (500 * sizeof (char));
	snprintf (file_data_filename, 500, "%s/polymerStates.timeseries", folderName);
	file_data = fopen (file_data_filename, "r");

	char *file_block_nBridges_filename;
	file_block_nBridges_filename = (char *) malloc (500 * sizeof (char));
	snprintf (file_block_nBridges_filename, 500, "%s/bridges.average.block", folderName);
	file_block_nBridges = fopen (file_block_nBridges_filename, "w");

	char *file_block_nLoops_filename;
	file_block_nLoops_filename = (char *) malloc (500 * sizeof (char));
	snprintf (file_block_nLoops_filename, 500, "%s/loops.average.block", folderName);
	file_block_nLoops = fopen (file_block_nLoops_filename, "w");

	char *file_block_nDangles_filename;
	file_block_nDangles_filename = (char *) malloc (500 * sizeof (char));
	snprintf (file_block_nDangles_filename, 500, "%s/dangles.average.block", folderName);
	file_block_nDangles = fopen (file_block_nDangles_filename, "w");

	char *file_block_nFree_filename;
	file_block_nFree_filename = (char *) malloc (500 * sizeof (char));
	snprintf (file_block_nFree_filename, 500, "%s/free.average.block", folderName);
	file_block_nFree = fopen (file_block_nFree_filename, "w");

	int nLines = N_TIMEFRAMES_TO_CONSIDER2;
	float *inputData_nBridges, *inputData_nLoops, *inputData_nDangles, *inputData_nFree;

	inputData_nBridges = (float *) malloc (nLines * sizeof (float));
	inputData_nLoops = (float *) malloc (nLines * sizeof (float));
	inputData_nDangles = (float *) malloc (nLines * sizeof (float));
	inputData_nFree = (float *) malloc (nLines * sizeof (float));

	printf("Reading bridges/loops/dangles/free from saved file...\n");
	saveInputData (&inputData_nBridges, &inputData_nLoops, &inputData_nDangles, &inputData_nFree, nLines, file_data);

	/*	
	BLOCKS *blockAverages_nBridges, *blockAverages_nLoops, *blockAverages_nDangles, *blockAverages_nFree;
	blockAverages_nBridges = (BLOCKS *) malloc (nLines * sizeof (BLOCKS));
	blockAverages_nLoops = (BLOCKS *) malloc (nLines * sizeof (BLOCKS));
	blockAverages_nDangles = (BLOCKS *) malloc (nLines * sizeof (BLOCKS));
	blockAverages_nFree = (BLOCKS *) malloc (nLines * sizeof (BLOCKS));

	blockAverages_nBridges = initializeBlocks (blockAverages_nBridges, nLines);
	blockAverages_nLoops = initializeBlocks (blockAverages_nLoops, nLines);
	blockAverages_nDangles = initializeBlocks (blockAverages_nDangles, nLines);
	blockAverages_nFree = initializeBlocks (blockAverages_nFree, nLines);

	printf("Computing bridges...\n");
	blockAverages_nBridges = computeBlockAverages (blockAverages_nBridges, nLines, inputData_nBridges);
	printf("Computing loops...\n");
	blockAverages_nLoops = computeBlockAverages (blockAverages_nLoops, nLines, inputData_nLoops);
	printf("Computing dangles...\n");
	blockAverages_nDangles = computeBlockAverages (blockAverages_nDangles, nLines, inputData_nDangles);
	printf("Computing free...\n");
	blockAverages_nFree = computeBlockAverages (blockAverages_nFree, nLines, inputData_nFree);

	// printing the results from block averaging of microstates
	printf("Printing block averages of bridges/loops/dangles/free...\n");
	printBlockAverageStats (file_block_nBridges, blockAverages_nBridges, nLines);
	printBlockAverageStats (file_block_nLoops, blockAverages_nLoops, nLines);
	printBlockAverageStats (file_block_nDangles, blockAverages_nDangles, nLines);
	printBlockAverageStats (file_block_nFree, blockAverages_nFree, nLines);

	free (blockAverages_nBridges);
	free (blockAverages_nLoops);
	free (blockAverages_nDangles);
	free (blockAverages_nFree);
	*/

	fclose (file_data);
	fclose (file_block_nBridges);
	fclose (file_block_nLoops);
	fclose (file_block_nDangles);
	fclose (file_block_nFree);

	/*
	To do:
		calculate coordination number
		calculate lifetime of loops/bridges/free/dangles
		calculate bridge to bridge transition time
	*/

	return 0;
}

/*
for dirname, dirpath, files in os.walk ("."):
    for file in files:
            if ("stat" in file and "dump" in file):
                    os.chdir (dirname)
                    os.system ("cp {}/completeStats.sh .".format (pd))
                    os.chdir (pd)
                    print (dirname)

FOR LATEX PARTICLE SIMULATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BL:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tauBL.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tauBL.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tauBL.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tauBL.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tauBL.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tauBL.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tauBL.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tauBL.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tauBL.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tauBL.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tauBL.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tauBL.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tauBL.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tauBL.output

BL2:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tauBL2.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tauBL2.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tauBL2.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tauBL2.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tauBL2.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tauBL2.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tauBL2.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tauBL2.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tauBL2.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tauBL2.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tauBL2.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tauBL2.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tauBL2.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tauBL2.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tauBL2.output

BB:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tauBB.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tauBB.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tauBB.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tauBB.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tauBB.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tauBB.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tauBB.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tauBB.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tauBB.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tauBB.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tauBB.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tauBB.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tauBB.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tauBB.output

LB:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tauLB.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tauLB.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tauLB.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tauLB.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tauLB.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tauLB.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tauLB.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tauLB.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tauLB.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tauLB.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tauLB.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tauLB.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tauLB.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tauLB.output

LL:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tauLL.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tauLL.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tauLL.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tauLL.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tauLL.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tauLL.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tauLL.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tauLL.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tauLL.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tauLL.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tauLL.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tauLL.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tauLL.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tauLL.output

tauE_bridge:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tauE_at_bridge.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tauE_at_bridge.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tauE_at_bridge.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tauE_at_bridge.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tauE_at_bridge.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tauE_at_bridge.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tauE_at_bridge.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tauE_at_bridge.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tauE_at_bridge.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tauE_at_bridge.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tauE_at_bridge.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tauE_at_bridge.output

BDBd:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tau_BDBd.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tau_BDBd.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tau_BDBd.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tau_BDBd.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tau_BDBd.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tau_BDBd.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tau_BDBd.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tau_BDBd.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tau_BDBd.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tau_BDBd.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tau_BDBd.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tau_BDBd.output

BDBs:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tau_BDBs.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tau_BDBs.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tau_BDBs.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tau_BDBs.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tau_BDBs.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tau_BDBs.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tau_BDBs.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tau_BDBs.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tau_BDBs.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tau_BDBs.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tau_BDBs.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tau_BDBs.output

BDL:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tau_BDL.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tau_BDL.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tau_BDL.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tau_BDL.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tau_BDL.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tau_BDL.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tau_BDL.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tau_BDL.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tau_BDL.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tau_BDL.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tau_BDL.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tau_BDL.output

LDB:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tau_LDB.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tau_LDB.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tau_LDB.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tau_LDB.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tau_LDB.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tau_LDB.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tau_LDB.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tau_LDB.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tau_LDB.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tau_LDB.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tau_LDB.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tau_LDB.output

LDL:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_0/tau_LDL.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_0/tau_LDL.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_0/tau_LDL.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_0/tau_LDL.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_0/tau_LDL.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_0/tau_LDL.output

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1;
wc -l eps5/stats_dump.short_data.poly_1_4_*_1/tau_LDL.output; wc -l eps6/stats_dump.short_data.poly_1_4_*_1/tau_LDL.output; wc -l eps7/stats_dump.short_data.poly_1_4_*_1/tau_LDL.output; wc -l eps8/stats_dump.short_data.poly_1_4_*_1/tau_LDL.output; wc -l eps9/stats_dump.short_data.poly_1_4_*_1/tau_LDL.output; wc -l eps10/stats_dump.short_data.poly_1_4_*_1/tau_LDL.output


Bridges:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1;

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1;


Loop:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1;

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1;


Dangle:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1;

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1;


Free:

head -10 eps5/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1;

head -10 eps5/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps6/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps7/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps8/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps9/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps10/stats_dump.short_data.poly_1_4_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1;


FOR GHOST PARTICLE SIMULATION:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BL:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauBL.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tauBL.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tauBL.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tauBL.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tauBL.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauBL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauBL.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tauBL.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tauBL.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tauBL.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tauBL.output

BB:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauBB.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tauBB.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tauBB.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tauBB.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tauBB.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauBB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauBB.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tauBB.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tauBB.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tauBB.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tauBB.output

LB:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauLB.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tauLB.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tauLB.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tauLB.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tauLB.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauLB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauLB.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tauLB.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tauLB.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tauLB.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tauLB.output

LL:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauLL.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tauLL.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tauLL.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tauLL.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tauLL.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauLL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauLL.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tauLL.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tauLL.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tauLL.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tauLL.output

tauE_bridge:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tauE_at_bridge.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tauE_at_bridge.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tauE_at_bridge.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tauE_at_bridge.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tauE_at_bridge.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tauE_at_bridge.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tauE_at_bridge.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tauE_at_bridge.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tauE_at_bridge.output

BDBd:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDBd.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDBd.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDBd.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDBd.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDBd.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDBd.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDBd.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDBd.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDBd.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDBd.output

BDBs:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDBs.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDBs.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDBs.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDBs.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDBs.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDBs.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDBs.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDBs.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDBs.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDBs.output

BDL:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_BDL.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDL.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDL.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDL.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tau_BDL.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_BDL.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDL.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDL.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDL.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tau_BDL.output

LDB:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_LDB.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tau_LDB.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tau_LDB.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tau_LDB.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tau_LDB.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_LDB.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tau_LDB.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tau_LDB.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tau_LDB.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tau_LDB.output

LDL:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_tau_LDL.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_0/tau_LDL.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_0/tau_LDL.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_0/tau_LDL.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_0/tau_LDL.output

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_tau_LDL.output_n999999_c1_2.block| tail -1;
wc -l eps7/stats_dump*.stat.gz_input.data_1_3_*_1/tau_LDL.output; wc -l eps8/stats_dump*.stat.gz_input.data_1_3_*_1/tau_LDL.output; wc -l eps9/stats_dump*.stat.gz_input.data_1_3_*_1/tau_LDL.output; wc -l eps10/stats_dump*.stat.gz_input.data_1_3_*_1/tau_LDL.output


Bridges:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c1_2.block | tail -1;

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c1_2.block | tail -1;


Loop:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c2_2.block | tail -1;

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c2_2.block | tail -1;


Dangle:

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c3_2.block | tail -1;

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c3_2.block | tail -1;


Free:

hhead -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_0/average_polymerStates.timeseries_n999999_c4_2.block | tail -1;

head -10 eps7/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps8/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps9/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1; head -10 eps10/stats_dump*.stat.gz_input.data_1_3_*_1/average_polymerStates.timeseries_n999999_c4_2.block | tail -1;

*/