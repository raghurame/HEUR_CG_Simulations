#ifndef BRIDGESBETWEENBINS_H
#define BRIDGESBETWEENBINS_H

BRIDGESBIN *assignBinBounds (BRIDGESBIN *bridgeBetweenBins, BOUNDARY simBoundary, float binWidth, float delBinDistance, int nBins);
BRIDGESBIN *countBridgesBetweenBins (TRAJECTORY **atoms, BOUNDARY simBoundary, float distanceCutoff, BRIDGESBIN *bridgeBetweenBins, int nAtoms, TRAJECTORY *micelles, int nMicelles, int nBins);

#endif
