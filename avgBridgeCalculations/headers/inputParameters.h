#ifndef STRUCTDEFINITIONS_H
#define STRUCTDEFINITIONS_H

#define MSLEEP usleep (1000000)
#define NTHREADS 8
#define INPUTFILE "dump.lammpstrj"
#define BINWIDTH_VERTICALBRIDGES 6.0
#define ADSORPTION_DISTANCE_CUTOFF 3.38
#define DELBINDISTANCE_VERTICALBRIDGES 0.1
#define MAX_FENE_EXTENSION 60.0
#define NBINS_YDISTRIBUTIONBRIDGES 20
#define BINWIDTH_BONDCENTERDISTRIBUTION 3.0
#define OUTPUT_BRIDGESBETWEENBINS "outputs/nBridgesBetweenBins.count"
#define OUTPUT_YDISTRIBUTIONBRIDGES "outputs/bridges.distribution"
#define OUTPUT_BONDCENTERDISTRIBUTION "outputs/bridgeCenter.distribution"
#define NFREQ_PROGRESS_SCREEN 5
#define OUTPUT_CURRENT_STATES "outputs/current.states"
#define OUTPUT_AVG_STATES "outputs/average.states"

#endif