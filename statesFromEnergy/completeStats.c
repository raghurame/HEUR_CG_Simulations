#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#define NPARTICLES 961
#define NPOLYMERS 10
#define NBEADS 2
#define COORDINATION_NUMBER 20

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

	rewind (inputDump);
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
				bondNumber = (energyEntries[i].atom2 - (energyEntries[i].atom2 / COORDINATION_NUMBER) - 1) / NBEADS;

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
			}
			else if (sortedAtoms[energyEntries[i].atom2 - 1].atomType == 2 && sortedAtoms[energyEntries[i].atom1 - 1].atomType == 1)
			{
				bondNumber = (energyEntries[i].atom1 - (energyEntries[i].atom1 / COORDINATION_NUMBER) - 1) / NBEADS;

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
			}
		}
	}

/*	for (int i = 0; i < datafile.nBonds; ++i)
	{
		if (boundParticleID[i] == 0)
		{
			polymerBondStatus[i][currentTimeframe].isItFree = true;
		}
	}
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

int main(int argc, char const *argv[])
{
	FILE *inputDump, *outputStates;
	inputDump = fopen (argv[1], "r");
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

	int nTimeframes = countNTimeframes (argv[1]);

	BOND_STATUS **polymerBondStatus;
	polymerBondStatus = (BOND_STATUS **) malloc (datafile.nBonds * sizeof (BOND_STATUS *));

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		polymerBondStatus[i] = (BOND_STATUS *) malloc ((nTimeframes + 1) * sizeof (BOND_STATUS));
	}

	polymerBondStatus = initBondStatus (polymerBondStatus, datafile.nBonds, nTimeframes);

	int file_status = fgetc (inputDump);
	char lineString[3000];

	DUMP_ENERGY *energyEntries;
	int nDumpEntries;
	nDumpEntries = countDumpEntries (inputDump, nDumpEntries);
	energyEntries = (DUMP_ENERGY *) malloc (nDumpEntries * 2 * sizeof (DUMP_ENERGY));
	int currentTimeframe = 0;

	STATES polymerStates;

	while (file_status > 0)
	{
		printf("Scanning timeframe: %d / %d                  \r", currentTimeframe + 1, nTimeframes);
		fflush (stdout);

		energyEntries = initEnergyEntries (energyEntries, nDumpEntries);
		energyEntries = saveDumpEnergyEntries (inputDump, energyEntries, nDumpEntries);

		polymerBondStatus = checkBondStatus (polymerBondStatus, energyEntries, nDumpEntries, datafile, nTimeframes, currentTimeframe, sortedAtoms);

		polymerStates = countStates (polymerStates, polymerBondStatus, datafile, currentTimeframe);

		if ((currentTimeframe + 1) < nTimeframes)
		{
			printStates (polymerStates, outputStates, datafile);
		}

		file_status = fgetc (inputDump);
		currentTimeframe++;
	}

	fclose (inputDump);
	fclose (outputStates);

	// block averaging
	FILE *file_data, *file_block_nBridges, *file_block_nLoops, *file_block_nDangles, *file_block_nFree;
	file_data = fopen ("polymerStates.timeseries", "r");
	file_block_nBridges = fopen ("bridges.average.block", "w");
	file_block_nLoops = fopen ("loops.average.block", "w");
	file_block_nDangles = fopen ("dangles.average.block", "w");
	file_block_nFree = fopen ("free.average.block", "w");

	int nLines = nTimeframes;
	float *inputData_nBridges, *inputData_nLoops, *inputData_nDangles, *inputData_nFree;

	inputData_nBridges = (float *) malloc (nLines * sizeof (float));
	inputData_nLoops = (float *) malloc (nLines * sizeof (float));
	inputData_nDangles = (float *) malloc (nLines * sizeof (float));
	inputData_nFree = (float *) malloc (nLines * sizeof (float));

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

	blockAverages_nBridges = computeBlockAverages (blockAverages_nBridges, nLines, inputData_nBridges);
	blockAverages_nLoops = computeBlockAverages (blockAverages_nLoops, nLines, inputData_nLoops);
	blockAverages_nDangles = computeBlockAverages (blockAverages_nDangles, nLines, inputData_nDangles);
	blockAverages_nFree = computeBlockAverages (blockAverages_nFree, nLines, inputData_nFree);

	printBlockAverageStats (file_block_nBridges, blockAverages_nBridges, nLines);
	printBlockAverageStats (file_block_nLoops, blockAverages_nLoops, nLines);
	printBlockAverageStats (file_block_nDangles, blockAverages_nDangles, nLines);
	printBlockAverageStats (file_block_nFree, blockAverages_nFree, nLines);

	fclose (file_data);
	fclose (file_block_nBridges);
	fclose (file_block_nLoops);
	fclose (file_block_nDangles);
	fclose (file_block_nFree);
	return 0;
}