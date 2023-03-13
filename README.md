# HEUR_CG_Simulations

The C code calculates the following,

1. Replicates the simulation box along the Y axis, to emulate the periodic boundary. Creates two bins, with a thickness along Y-axis. The yhi of lower bin is the same as ylo of the upper bin. The number of bridges (originating from one bin and ending in another bin) are then calculated between the two bins. Then the two bins are translated along the Y-axis, to measure the bridges as a function of Y.
2. Counts the number of bridges as a function of del(y), where del(y) equals y2 - y1 (the two polymer end groups are located at x1, y1, z1, and x2, y2, z2)

# To compile
Type 'make'

#To run
Type 'make run'
