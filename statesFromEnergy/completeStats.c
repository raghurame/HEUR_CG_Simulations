#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#define NPARTICLES 100
#define NPOLYMERS 4000
#define NBEADS 2
#define COORDINATION_NUMBER 80
#define N_TIMEFRAMES_TO_CONSIDER 200

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

int *storeParticleIDs (int *particleIDs)
{
	for (int i = 0; i < NPARTICLES; ++i)
	{
		particleIDs[i] = (i + 1) + (i + 1) * COORDINATION_NUMBER;
	}

	return particleIDs;
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

int countNTimeframes (const char *filename)
{
	int nTimeframes = 0;
	char *pipeString, lineString[300];
	pipeString = (char *) malloc (1000 * sizeof (char));
	snprintf (pipeString, 1000, "cat %s | grep \"ITEM: TIMESTEP\" | wc -l", filename);
	FILE *lineCount;
	lineCount = popen (pipeString, "r");

	fgets (lineString, 300, lineCount);
	sscanf (lineString, "%d", &nTimeframes);
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

float *countTauBL (float *tauBL, BOND_STATUS **polymerBondStatus, int nTimeframes, DATAFILE_INFO datafile)
{
	int currentTransition = 0, currentTau = 0;
	printf("\nCalculating transition time from bridge to loop...\n");

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		currentTau = 0;

		for (int j = 0; j < nTimeframes; ++j)
		{

			// printStateNumber (polymerBondStatus[i][j].isItFree, polymerBondStatus[i][j].isItDangle, polymerBondStatus[i][j].isItLoop, polymerBondStatus[i][j].isItBridge);

			if (polymerBondStatus[i][j].isItBridge == 1)
			{
				currentTau++;
				// printf("(%d) ", currentTau);
				// usleep (10000);
			}
			else if (j > 0)
			{
				if (polymerBondStatus[i][j - 1].isItBridge == 1 && polymerBondStatus[i][j].isItLoop == 1)
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

void printTauBL (float *tauBL, int nBL)
{
	FILE *tauBL_file;
	tauBL_file = fopen ("tauBL.output", "w");

	for (int i = 0; i < nBL; ++i)
	{
		fprintf(tauBL_file, "%f\n", tauBL[i]);
	}

	fclose (tauBL_file);
}

void printTauLB (float *tauLB, int nLB)
{
	FILE *tauLB_file;
	tauLB_file = fopen ("tauLB.output", "w");

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

int main(int argc, char const *argv[])
{
	if (argc != 3)
	{
		printf("ERROR: INCORRECT ARGUMENTS PASSED.\n\n Required arguments:\n\n{~} argv[0] = ./program\n{~} argv[1] = input dump file (ascii text or *.gz)\n{~} argv[2] = input data file.\n\n");
		exit (1);
	}

	FILE *inputDump, *outputStates;
	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));
	snprintf (pipeString, 500, "zcat %s", argv[1]);

	if (strstr (argv[1], ".gz"))
	{
		inputDump = popen (pipeString, "r");
	}
	else
	{
		inputDump = fopen (argv[1], "r");
	}

	outputStates = fopen ("polymerStates.timeseries", "w");

	int nChains = NPARTICLES * NPOLYMERS, *particleIDs, nAtoms = NPARTICLES + (NPARTICLES * NPOLYMERS * NBEADS);
	particleIDs = (int *) malloc (NPARTICLES * sizeof (int));
	particleIDs = storeParticleIDs (particleIDs);

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

	int nTimeframes = countNTimeframes (argv[1]), N_TIMEFRAMES_TO_CONSIDER2;

	if (N_TIMEFRAMES_TO_CONSIDER >= nTimeframes)
	{
		N_TIMEFRAMES_TO_CONSIDER2 = nTimeframes;
	}
	else
	{
		N_TIMEFRAMES_TO_CONSIDER2 = N_TIMEFRAMES_TO_CONSIDER;
	}

	BOND_STATUS **polymerBondStatus;
	polymerBondStatus = (BOND_STATUS **) malloc (datafile.nBonds * sizeof (BOND_STATUS *));

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		polymerBondStatus[i] = (BOND_STATUS *) malloc ((N_TIMEFRAMES_TO_CONSIDER2 + 1) * sizeof (BOND_STATUS));
	}

	polymerBondStatus = initBondStatus (polymerBondStatus, datafile.nBonds, N_TIMEFRAMES_TO_CONSIDER2);

	int file_status = fgetc (inputDump);
	char lineString[3000];

	DUMP_ENERGY *energyEntries;
	int nDumpEntries;
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
	int currentTimeframe = 0;

	STATES polymerStates;

	int timeframesToSkip = 0;

	if (nTimeframes > N_TIMEFRAMES_TO_CONSIDER2) {
		timeframesToSkip = nTimeframes - N_TIMEFRAMES_TO_CONSIDER2; }

	printf("Skipping %d timeframes initially...\n", timeframesToSkip);
	printf("To consider: %d\n", N_TIMEFRAMES_TO_CONSIDER2);

	while (file_status > 0)
	{
		printf("Scanning timeframe: %d / %d                  \r", currentTimeframe + 1, nTimeframes);
		fflush (stdout);

		energyEntries = initEnergyEntries (energyEntries, nDumpEntries);
		energyEntries = saveDumpEnergyEntries (inputDump, energyEntries, nDumpEntries);

		if ((currentTimeframe + 1) > timeframesToSkip)
		{
			polymerBondStatus = checkBondStatus (polymerBondStatus, energyEntries, nDumpEntries, datafile, N_TIMEFRAMES_TO_CONSIDER2, currentTimeframe - timeframesToSkip, sortedAtoms);

			polymerStates = countStates (polymerStates, polymerBondStatus, datafile, currentTimeframe - timeframesToSkip);

			if ((currentTimeframe + 1) < nTimeframes) {
				printStates (polymerStates, outputStates, datafile); }
		}

		if (currentTimeframe > nTimeframes) {
			goto leaveThisLoop; }

		file_status = fgetc (inputDump);
		currentTimeframe++;
	}

	leaveThisLoop:;

	free (energyEntries);

	if (strstr (argv[1], ".gz"))
	{
		pclose (inputDump);
	}
	else
	{
		fclose (inputDump);
	}

	fclose (outputStates);

	// computing transitions
	polymerBondStatus = correctingDangles (polymerBondStatus, datafile, N_TIMEFRAMES_TO_CONSIDER2);

	float *tauBL, *tauLB;
	int nBL = countBLtransitions (nBL, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile), nLB = countLBtransitions (nLB, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile);

	tauBL = (float *) malloc (nBL * sizeof (float));
	tauLB = (float *) malloc (nLB * sizeof (float));

	tauBL = countTauBL (tauBL, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile);
	tauLB = countTauLB (tauLB, polymerBondStatus, N_TIMEFRAMES_TO_CONSIDER2, datafile);

	printTauBL (tauBL, nBL);
	printTauLB (tauLB, nLB);

	float averageLB, stdevLB, stderrLB;
	float averageBL, stdevBL, stderrBL;
	computeStats (&averageLB, &stdevLB, &stderrLB, tauLB, nLB);
	computeStats (&averageBL, &stdevBL, &stderrBL, tauBL, nBL);

	FILE *transitionStats;
	transitionStats = fopen ("transitions.stats", "w");
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
	file_data = fopen ("polymerStates.timeseries", "r");
	file_block_nBridges = fopen ("bridges.average.block", "w");
	file_block_nLoops = fopen ("loops.average.block", "w");
	file_block_nDangles = fopen ("dangles.average.block", "w");
	file_block_nFree = fopen ("free.average.block", "w");

	int nLines = N_TIMEFRAMES_TO_CONSIDER2;
	float *inputData_nBridges, *inputData_nLoops, *inputData_nDangles, *inputData_nFree;

	inputData_nBridges = (float *) malloc (nLines * sizeof (float));
	inputData_nLoops = (float *) malloc (nLines * sizeof (float));
	inputData_nDangles = (float *) malloc (nLines * sizeof (float));
	inputData_nFree = (float *) malloc (nLines * sizeof (float));

	printf("Reading bridges/loops/dangles/free from saved file...\n");
	saveInputData (&inputData_nBridges, &inputData_nLoops, &inputData_nDangles, &inputData_nFree, nLines, file_data);

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

	fclose (file_data);
	fclose (file_block_nBridges);
	fclose (file_block_nLoops);
	fclose (file_block_nDangles);
	fclose (file_block_nFree);
	return 0;
}