#ifndef COMPUTESTATESHEUR_H
#define COMPUTESTATESHEUR_H

void printAverageStates (const char *filename_avgStates, STATES avgStates, STATES stdevStates);
STATES computeStdevStates (STATES stdevStates, STATES avgStates, STATES *allStates, int nTimeframes);
STATES *readAllStates (STATES *allStates, int nTimeframes, const char *filename_states);
STATES computeAvgStates (STATES avgStates, int nTimeframes);
void printCurrentStates (FILE *file_printStates, STATES currentStates);
STATES sumAllStates (STATES currentStates, STATES avgStates);
STATES computeAllStates (STATES currentStates, TRAJECTORY *atoms, int nAtoms);
STATES initializeStates (STATES currentStates);

#endif