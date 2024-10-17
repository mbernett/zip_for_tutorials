#!/bin/bash
if [ -z "$VMDPATH" ]; then
   echo "VMDPATH env var is not set, source the env.sh in the top folder first"
   exit
fi
if [ -z "$PATHTOOLS" ]; then
   echo "PATHTOOLS env var is not set, source the env.sh in the top folder first"
   exit
fi

PDB=./structure_protlig.pdb
TRAJ=./steer_protlig.xtc
PROTSELECTION='resid 47 48 49 50 51 52 53 54 86 87 88 89 90 91 92 93 179 180 181 182 184 185 221 222 225 and name CA'
LIGSELECTION='resname ALP'

MULTIPDB="multipdb_occbeta.pdb"
MORPHED="REPARAM_FINAL_morphed.pdb"

OPTS="-dispdev text"

# Generate multipdb for rmsd_mc.sh
echo "Tagging occupancy and b-factor ...."
$VMDPATH $OPTS -e $PATHTOOLS/occ_beta.tcl -args -pdb $PDB -traj $TRAJ \
  -protsel "${PROTSELECTION// /__}"   -ligsel "${LIGSELECTION// /__}"  \
  -outpdb $MULTIPDB >/dev/null
echo "Done!"
echo ""

# Extract optimal frames: round 1
echo "Extracting optimal frames ...."
$PATHTOOLS/rmsd_mc.sh -f $MULTIPDB -r 80 -s 1 -t 1 -n 10000 -c 1000 > out_rmsdmc.dat
echo "Done!"
echo ""

echo "Interpolating ...."
$VMDPATH $OPTS -e $PATHTOOLS/morph_frameset.tcl -args -pdb REPARAM_FINAL.pdb \
  -morph 5 -outpdb $MORPHED >/dev/null
echo "Done!"
echo ""

# Extract optimal frames: round 2
echo "Extracting optimal frames ...."
$PATHTOOLS/rmsd_mc.sh -f $MORPHED -r 80 -s 1 -t 1 -n 10000 -c 1000 > out_rmsdmc.dat
echo "Done!"


