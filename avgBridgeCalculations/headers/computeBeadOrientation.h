#ifndef COMPUTEBEDORIENTATIONHEUR_H
#define COMPUTEBEDORIENTATIONHEUR_H

ANGLE_DISTRIBUTION *assignBeadOrientation (ANGLE_DISTRIBUTION *beadOrientation, int nBins_orientation);
ANGLE_DISTRIBUTION *computeBridgeOrientationDistribution (ANGLE_DISTRIBUTION *bridgeOrientation, int nBins_orientation, int nBonds, BONDINFO *allBonds);
ANGLE_DISTRIBUTION *computeLoopOrientationDistribution (ANGLE_DISTRIBUTION *loopOrientation, int nBins_orientation, int nBonds, BONDINFO *allBonds);
ANGLE_DISTRIBUTION *computeDangleOrientationDistribution (ANGLE_DISTRIBUTION *dangleOrientation, int nBins_orientation, int nBonds, BONDINFO *allBonds);
ANGLE_DISTRIBUTION *computeFreeOrientationDistribution (ANGLE_DISTRIBUTION *freeOrientation, int nBins_orientation, int nBonds, BONDINFO *allBonds);

#endif