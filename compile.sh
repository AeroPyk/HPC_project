#!/bin/bash

cmake CMakeLists.txt
make

# srun -n 4 --ntasks-per-node=1 ./main > res.txt # example