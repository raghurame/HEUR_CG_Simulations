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
	bool isItLoop, isItBridge;
} BOND_STATUS;

typedef struct dumpEnergy
{
	int atom1, atom2;
	float distance, energy;
} DUMP_ENERGY;

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

	for (int i = 0; i < 9; ++i)
	{
		fgets (lineString, 3000, inputDump);
	}

	for (int i = 0; i < nDumpEntries; ++i)
	{
		fgets (lineString, 3000, inputDump);
		sscanf (lineString, "%d %d %f %f\n", &energyEntries[i].atom1, &energyEntries[i].atom2, &energyEntries[i].distance, &energyEntries[i].energy);
	}

	return energyEntries;
}

BOND_STATUS **checkBondStatus (BOND_STATUS **polymerBondStatus, DUMP_ENERGY *energyEntries, int nDumpEntries, DATAFILE_INFO datafile, int nTimeframes, int currentTimeframe)
{
	int bondNumber = 0;

	for (int i = 0; i < nDumpEntries; ++i)
	{
		if (energyEntries[i].energy < -1)
		{
			bondNumber = (energyEntries[i].atom1 - (energyEntries[i].atom1 / COORDINATION_NUMBER) - 1) / NBEADS;

			if (polymerBondStatus[bondNumber][currentTimeframe].isItLoop == false)
			{
				polymerBondStatus[bondNumber][currentTimeframe].isItLoop = true;
			}
			else
			{
				polymerBondStatus[bondNumber][currentTimeframe].isItBridge = true;
			}
		}
	}

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
		}
	}

	return polymerBondStatus;
}

int main(int argc, char const *argv[])
{
	FILE *inputDump;
	inputDump = fopen (argv[1], "r");

	int nChains = NPARTICLES * NPOLYMERS, *particleIDs, nAtoms = NPARTICLES + (NPARTICLES * NPOLYMERS * NBEADS);
	particleIDs = (int *) malloc (NPARTICLES * sizeof (int));
	particleIDs = storeParticleIDs (particleIDs);

	DATA_ATOMS *dataAtoms;
	DATA_BONDS *dataBonds, *sortedBonds;
	DATA_ANGLES *dataAngles;
	DATA_DIHEDRALS *dataDihedrals;
	DATA_IMPROPERS *dataImpropers;
	DATAFILE_INFO datafile;

	datafile = readData (argv[2], &dataAtoms, &dataBonds, &dataAngles, &dataDihedrals, &dataImpropers);
	sortedBonds = (DATA_BONDS *) malloc (datafile.nBonds * sizeof (DATA_BONDS));
	sortedBonds = sortBonds (sortedBonds, dataBonds, datafile.nBonds, datafile.nAtoms);

/*	for (int i = 0; i < datafile.nBonds; ++i)
	{
		printf("%d -> %d %d %d %d %d\n", i, sortedBonds[i].id, sortedBonds[i].bondType, sortedBonds[i].atom1, sortedBonds[i].atom2, (sortedBonds[i].atom2 - (sortedBonds[i].atom2 / COORDINATION_NUMBER) - 1) / NBEADS);
		usleep (100000);
	}
*/
	int nTimeframes = countNTimeframes (argv[1]);

	BOND_STATUS **polymerBondStatus;
	polymerBondStatus = (BOND_STATUS **) malloc (datafile.nBonds * sizeof (BOND_STATUS *));

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		polymerBondStatus[i] = (BOND_STATUS *) malloc (nTimeframes * sizeof (BOND_STATUS));
	}

	polymerBondStatus = initBondStatus (polymerBondStatus, datafile.nBonds, nTimeframes);

	int file_status = fgetc (inputDump);
	char lineString[3000];

	DUMP_ENERGY *energyEntries;
	int nDumpEntries = countDumpEntries (inputDump, nDumpEntries);
	energyEntries = (DUMP_ENERGY *) malloc (nDumpEntries * sizeof (DUMP_ENERGY));
	int currentTimeframe = 0;

	while (file_status)
	{
		printf("Scanning timeframe: %d                   \r", currentTimeframe + 1);
		fflush (stdout);

		energyEntries = saveDumpEnergyEntries (inputDump, energyEntries, nDumpEntries);

		polymerBondStatus = checkBondStatus (polymerBondStatus, energyEntries, nDumpEntries, datafile, nTimeframes, currentTimeframe);

		file_status = fgetc (inputDump);
		currentTimeframe++;
	}

	fclose (inputDump);
	return 0;
}