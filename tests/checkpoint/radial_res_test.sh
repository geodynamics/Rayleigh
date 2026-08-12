#!/usr/bin/env bash

# This test ensures that Rayleigh handles changes in radial resolution
# as expected and ensures compatibility with checkpoint formats from 
# version 1.3 and prior.

export RA_ROOT=../../..
export rayleigh_v13=$RA_ROOT/Rayleigh_1.3/bin/rayleigh.dbg
export rayleigh_dbg=$RA_ROOT/bin/rayleigh.dbg

cp ../../post_processing/rayleigh_diagnostics.py .

# Generate a checkpoint in the old format using Rayleigh v1.3
cd old_format
echo ""
echo "Now in: "$PWD
echo ""
mpirun -np 4 $rayleigh_v13 -niter 100
cd ..

# Restart from the old checkpoint format at increased resolution using Rayleigh 1.3.
lastdir=old_format
newdir=old_format_up

mkdir $newdir
mkdir $newdir/Checkpoints
cp $lastdir/main_input $newdir/.
cp -r $lastdir/Checkpoints/00000050 $newdir/Checkpoints/.

cd $newdir
echo ""
echo "Now in: "$PWD
echo ""
sed -i 's/n_r=48/n_r=64/g' main_input
sed -i 's/init_type=1/init_type=-1/g' main_input
sed -i 's/restart_iter=-1/restart_iter=50/g' main_input
mpirun -np 4 $rayleigh_v13 -niter 100
cd ..


# Restart from the old checkpoint format at same resolution using current Rayleigh.
lastdir=old_format
newdir=new_format

mkdir $newdir
mkdir $newdir/Checkpoints
cp $lastdir/main_input $newdir/.
cp -r $lastdir/Checkpoints/00000050 $newdir/Checkpoints/.

cd $newdir
echo ""
echo "Now in: "$PWD
echo ""
sed -i 's/init_type=1/init_type=-1/g' main_input
sed -i 's/restart_iter=-1/restart_iter=50/g' main_input
mpirun -np 4 $rayleigh_dbg -niter 100
cd ..


# Restart from the new checkpoint format at same resolution.
lastdir=new_format
newdir=new_format_same_res
mkdir $newdir
mkdir $newdir/Checkpoints
cp $lastdir/main_input $newdir/.
cp -r $lastdir/Checkpoints/00000100 $newdir/Checkpoints/.

cd $newdir
echo ""
echo "Now in: "$PWD
echo ""
sed -i 's/restart_iter=50/restart_iter=100/g' main_input
mpirun -np 4 $rayleigh_dbg -niter 100
cd ..

# Restart from the new checkpoint format at increased resolution.
lastdir=new_format_same_res
newdir=new_format_up_res
mkdir $newdir
mkdir $newdir/Checkpoints
cp $lastdir/main_input $newdir/.
cp -r $lastdir/Checkpoints/00000150 $newdir/Checkpoints/.

cd $newdir
echo ""
echo "Now in: "$PWD
echo ""
sed -i 's/n_r=48/n_r=64/g' main_input
sed -i 's/restart_iter=100/restart_iter=150/g' main_input
mpirun -np 4 $rayleigh_dbg -niter 100
cd ..

# Restart from the new checkpoint format and degrade resolution.
lastdir=new_format_up_res
newdir=new_format_down_res
mkdir $newdir
mkdir $newdir/Checkpoints
cp $lastdir/main_input $newdir/.
cp -r $lastdir/Checkpoints/00000200 $newdir/Checkpoints/.

cd $newdir
echo ""
echo "Now in: "$PWD
echo ""
sed -i 's/n_r=64/n_r=48/g' main_input
sed -i 's/restart_iter=150/restart_iter=200/g' main_input
mpirun -np 4 $rayleigh_dbg -niter 100
cd ..

python3 compare_check.py
