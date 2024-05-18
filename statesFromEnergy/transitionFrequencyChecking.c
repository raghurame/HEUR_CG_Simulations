#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

#define NPARTICLES 100
#define NPOLYMERS 4000
#define NBEADS 2
#define COORDINATION_NUMBER 80
#define N_TIMEFRAMES_TO_CONSIDER 20000000

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
	int id;
	bool isItLoop, isItBridge, isItDangle, isItFree;
	bool loop_corrected, bridge_corrected;
} BOND_STATUS;

typedef struct dumpEnergy
{
	int atom1, atom2;
	float distance, energy;
} DUMP_ENERGY;

typedef struct states
{
	int nLoops, nBridges, nFree, nDangles;
} STATES;

int countNTimeframes (const char *filename)
{
	int nTimeframes = 0;
	char *pipeString, lineString[500];
	pipeString = (char *) malloc (1000 * sizeof (char));

	if (strstr (filename, ".gz"))
	{
		snprintf (pipeString, 1000, "zcat %s | grep \"ITEM: TIMESTEP\" | wc -l", filename);
	}
	else
	{
		snprintf (pipeString, 1000, "cat %s | grep \"ITEM: TIMESTEP\" | wc -l", filename);
	}

	FILE *lineCount;
	lineCount = popen (pipeString, "r");

	fgets (lineString, 500, lineCount);
	sscanf (lineString, "%d", &nTimeframes);
	printf("Number of timeframes: %d\n", nTimeframes);

	return nTimeframes;
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

DUMP_ENERGY *saveDumpEnergyEntries (FILE *inputDump, DUMP_ENERGY *energyEntries, int nDumpEntries)
{
	char lineString[3000];
	int nDumpEntries_current;

	for (int i = 0; i < 9; ++i)
	{
		fgets (lineString, 3000, inputDump);

		if (i == 3) {
			sscanf (lineString, "%d\n", &nDumpEntries_current); }
	}

	for (int i = 0; i < nDumpEntries_current; ++i)
	{
		fgets (lineString, 3000, inputDump);
		sscanf (lineString, "%d %d %f %f\n", &energyEntries[i].atom1, &energyEntries[i].atom2, &energyEntries[i].distance, &energyEntries[i].energy);
	}

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

BOND_STATUS **checkBondStatus (BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, int nDumpEntries, DATAFILE_INFO datafile, int nTimeframes, int currentTimeframe, DATA_ATOMS *sortedAtoms)
{
	int bondNumber = 0;

	int *boundParticleID;
	boundParticleID = (int *) malloc (datafile.nBonds * sizeof (int));

	boundParticleID = initBoundParticleIDs (boundParticleID, datafile.nBonds);

	for (int i = 0; i < nDumpEntries; ++i)
	{
		// usleep (100000);
		// fprintf(stdout, "%d %d %f ", energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].energy);

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
				// bondNumber = (energyEntries[i].atom2 - (energyEntries[i].atom2 / COORDINATION_NUMBER) - 1) / NBEADS;
				bondNumber = floor ((energyEntries[i].atom2 - (int)floor (energyEntries[i].atom2 / (COORDINATION_NUMBER + 1)) - 1) / 2);
				// printf("(%d) ", bondNumber);

				if (boundParticleID[bondNumber] == 0)
				{
					boundParticleID[bondNumber] = energyEntries[i].atom1;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = true;
				}
				else if (boundParticleID[bondNumber] == energyEntries[i].atom1)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItLoop = true;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = false;
				}
				else if (boundParticleID[bondNumber] != energyEntries[i].atom1)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItBridge = true;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = false;
				}

				// printf("(%d) --> ", boundParticleID[bondNumber]);

				// printStateNumber (polymerBondStatus[bondNumber][currentTimeframe].isItFree, polymerBondStatus[bondNumber][currentTimeframe].isItDangle, polymerBondStatus[bondNumber][currentTimeframe].isItLoop, polymerBondStatus[bondNumber][currentTimeframe].isItBridge);
				// printf("\n");
			}
			else if (sortedAtoms[energyEntries[i].atom2 - 1].atomType == 2 && sortedAtoms[energyEntries[i].atom1 - 1].atomType == 1)
			{
				// bondNumber = (energyEntries[i].atom1 - (energyEntries[i].atom1 / COORDINATION_NUMBER) - 1) / NBEADS;
				bondNumber = floor ((energyEntries[i].atom1 - (int)floor (energyEntries[i].atom1 / (COORDINATION_NUMBER + 1)) - 1) / 2);
				// printf("(%d) ", bondNumber);

				if (boundParticleID[bondNumber] == 0)
				{
					boundParticleID[bondNumber] = energyEntries[i].atom2;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = true;
				}
				else if (boundParticleID[bondNumber] == energyEntries[i].atom2)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItLoop = true;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = false;
				}
				else if (boundParticleID[bondNumber] != energyEntries[i].atom2)
				{
					polymerBondStatus[bondNumber][currentTimeframe].isItBridge = true;
					polymerBondStatus[bondNumber][currentTimeframe].isItDangle = false;
				}

				// printf("(%d) --> ", boundParticleID[bondNumber]);

				// printStateNumber (polymerBondStatus[bondNumber][currentTimeframe].isItFree, polymerBondStatus[bondNumber][currentTimeframe].isItDangle, polymerBondStatus[bondNumber][currentTimeframe].isItLoop, polymerBondStatus[bondNumber][currentTimeframe].isItBridge);
				// printf("\n");
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
		printf("%d -> ", i);
		printStateNumber (polymerBondStatus[i][currentTimeframe].isItFree, polymerBondStatus[i][currentTimeframe].isItDangle, polymerBondStatus[i][currentTimeframe].isItLoop, polymerBondStatus[i][currentTimeframe].isItBridge);
		printf("\n");
	}

	sleep (100);
*/
	return polymerBondStatus;
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

int main(int argc, char const *argv[])
{
	if (argc != 4)
	{
		printf("ERROR: INCORRECT ARGUMENTS PASSED.\n\n Required arguments:\n\n{~} argv[0] = ./program\n{~} argv[1] = input dump file (ascii text or *.gz)\n{~} argv[2] = input data file.\n{~} argv[3] = max atom number in input dump file.\n\n");
		exit (1);
	}

	FILE *input;
	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));
	snprintf (pipeString, 500, "zcat %s", argv[1]);

	if (strstr (argv[1], ".gz"))
	{
		input = popen (pipeString, "r");
	}
	else
	{
		input = fopen (argv[1], "r");
	}

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

	int file_status = fgetc (input);
	char lineString[3000];
	int nTimeframes = countNTimeframes (argv[1]);

/*	BOND_STATUS **polymerBondStatus;
	polymerBondStatus = (BOND_STATUS **) malloc (datafile.nBonds * sizeof (BOND_STATUS *));

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		polymerBondStatus[i] = (BOND_STATUS *) malloc ((nTimeframes + 1) * sizeof (BOND_STATUS));
	}

	polymerBondStatus = initBondStatus (polymerBondStatus, datafile.nBonds, nTimeframes);
*/
	int nDumpEntries;
	nDumpEntries = countDumpEntries (input, nDumpEntries);

	DUMP_ENERGY *energyEntries;
	energyEntries = (DUMP_ENERGY *) malloc (nDumpEntries * 2 * sizeof (DUMP_ENERGY));

	int maxAtomID = atoi (argv[3]);

	DUMP_ENERGY **beadEnergies;
	beadEnergies = (DUMP_ENERGY **) malloc (maxAtomID * sizeof (DUMP_ENERGY *));

	for (int i = 0; i < maxAtomID; ++i)
	{
		beadEnergies[i] = (DUMP_ENERGY *) malloc (nTimeframes * sizeof (DUMP_ENERGY));
	}

	int currentTimeframe = 0;

	while (file_status > 0)
	{
		printf("Scanning timeframe: %d / %d                  \r", currentTimeframe + 1, nTimeframes);
		fflush (stdout);

		energyEntries = initEnergyEntries (energyEntries, nDumpEntries);
		energyEntries = saveDumpEnergyEntries (input, energyEntries, nDumpEntries);

		// polymerBondStatus = checkBondStatus (polymerBondStatus, energyEntries, nDumpEntries, datafile, nTimeframes, currentTimeframe, sortedAtoms);

		file_status = fgetc (input);
		currentTimeframe++;
	}

	fclose (input);
	return 0;
}