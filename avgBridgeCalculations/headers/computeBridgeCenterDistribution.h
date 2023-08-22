#ifndef BRIDGECENTERDISTRIBUTION_H
#define BRIDGECENTERDISTRIBUTION_H

BRIDGESBIN *assignBridgeCenterDistribution (BRIDGESBIN *bridgeCenterDistribution, int nBins_centerDistribution, float binWidth_centerDistribution, BOUNDARY simBoundary);
BONDINFO *computeBridgeCenter (TRAJECTORY *atoms, int nAtoms, BONDINFO *allBonds, BOUNDARY simBoundary, FILE *file_bondDump);
BRIDGESBIN *computeBridgeCenterDistribution (BONDINFO *allBonds, int nBonds, BRIDGESBIN *bridgeCenterDistribution, int nBins_centerDistribution);

#endif
