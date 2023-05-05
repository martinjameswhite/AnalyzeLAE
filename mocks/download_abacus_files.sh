#!/bin/bash
#
# Downloads the Abacus halo and particle data from NERSC HPSS.
# Should be run from within an AbacusSummit_high_c000_ph100
# subdirectory on e.g. $SCRATCH.
#
db=/nersc/projects/desi/cosmosim/Abacus/AbacusSummit_high_c000_ph100/
fn=Abacus_AbacusSummit_high_c000_ph100_halos.tar
htar -xvf ${db}${fn} './halos/z3.000'
#
