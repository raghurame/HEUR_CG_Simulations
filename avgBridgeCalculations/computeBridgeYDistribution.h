#ifndef BRIDGEYDISTRIBUTION_H
#define BRIDGEYDISTRIBUTION_H

YDIST *assignBridgeYDistribution (float maxFeneExtension, int nBins, float binWidth, YDIST *bridgeYDistribution);
YDIST *computeBridgeDistribution (TRAJECTORY *atoms, int nAtoms, YDIST *bridgeYDistribution, int nBins, BOUNDARY simBoundary);


#endif
