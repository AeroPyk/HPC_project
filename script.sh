#!/bin/bash

cmake CMakeLists.txt
make

srun -n 16 ./main