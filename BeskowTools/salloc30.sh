#!/bin/bash

# Easy way to perform salloc for 10 min
# 1 node with Hashwell processor

salloc --nodes=1 -t 00:30:00 -A edu20.DD2356 -C Haswell

# Then run with srun ./yourprogram [-n amountOfTask]

