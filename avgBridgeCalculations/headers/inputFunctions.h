#ifndef INPUTFUNCTIONSHEUR_H
#define INPUTFUNCTIONSHEUR_H

BOUNDARY getBoundary (FILE *file_inputTrj, BOUNDARY simBoundary);
int countNAtoms (FILE *file_inputTrj);
TRAJECTORY *getAtoms (TRAJECTORY *atoms, int nAtoms, BOUNDARY *simBoundary, float distanceCutoff, FILE *file_inputTrj, int file_status, TRAJECTORY **micelles, int nMicelles);
int countNMicelles (FILE *file_inputTrj, int nAtoms);

#endif