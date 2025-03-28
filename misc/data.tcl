set sel [atomselect top all]
measure minmax atomselect0
pbc set {150.0 150.0 150.0 90.0 90.0 90.0}
topo writelammpsdata data.run bond
