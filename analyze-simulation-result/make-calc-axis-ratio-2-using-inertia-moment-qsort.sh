#!/bin/csh

mpicxx -O3  -c ./mysorcecodes-for-FDPS/SPH.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -c ./mysorcecodes-for-FDPS/space.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -c ./mysorcecodes-for-FDPS/calc-pressure-force.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -c ./mysorcecodes-for-FDPS/calc-self-gravity.c -I ./mysorcecodes-for-FDPS
	
mpicxx -O3  -c ./mysorcecodes-for-FDPS/input.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -c ./mysorcecodes-for-FDPS/mirror.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -c ./mysorcecodes-for-FDPS/time-development.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -c ./mysorcecodes-for-FDPS/tillotson.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -c ./mysorcecodes-for-FDPS/make-initial-condition-file.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -c ./mysorcecodes-for-FDPS/fracture-and-porosity-model.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -c calc-axis-ratio-2-using-inertia-moment-qsort.c -I ./mysorcecodes-for-FDPS

mpicxx -O3  -o calc-axis-ratio-2-using-inertia-moment-qsort calc-axis-ratio-2-using-inertia-moment-qsort.o ./SPH.o ./space.o ./calc-pressure-force.o ./calc-self-gravity.o ./input.o ./mirror.o ./time-development.o ./tillotson.o ./make-initial-condition-file.o ./fracture-and-porosity-model.o

rm *.o
