#ifndef INPUTPARAMETERSHEUR_H
#define INPUTPARAMETERSHEUR_H

#define PI 3.14159
#define MSLEEP usleep (1000000)
#define NTHREADS 1
#define INPUTFILE "dump.lammpstrj"
#define BINWIDTH_VERTICALBRIDGES 6.0
#define ADSORPTION_DISTANCE_CUTOFF 3.38
#define DELBINDISTANCE_VERTICALBRIDGES 0.1
#define MAX_FENE_EXTENSION 60.0
#define NBINS_YDISTRIBUTIONBRIDGES 20
#define BINWIDTH_BONDCENTERDISTRIBUTION 3.0
#define NFREQ_PROGRESS_SCREEN 1
#define MAXTIMESTEPS 0
#define BINWIDTH_ANGLEORIENTATION 10.0
#define OUTPUT_BRIDGESBETWEENBINS "outputs/nBridgesBetweenBins.count"
#define OUTPUT_YDISTRIBUTIONBRIDGES "outputs/bridges.distribution"
#define OUTPUT_BONDCENTERDISTRIBUTION "outputs/bridgeCenter.distribution"
#define OUTPUT_CURRENT_STATES "outputs/current.states"
#define OUTPUT_AVG_STATES "outputs/average.states"
#define OUTPUT_BRIDGEORIENTATIONDISTRIBUTION "outputs/bridges.orientationDistribution"
#define OUTPUT_LOOPORIENTATIONDISTRIBUTION "outputs/loops.orientationDistribution"
#define OUTPUT_FREEORIENTATIONDISTRIBUTION "outputs/free.orientationDistribution"
#define OUTPUT_DANGLEORIENTATIONDISTRIBUTION "outputs/dangle.orientationDistribution"

#endif