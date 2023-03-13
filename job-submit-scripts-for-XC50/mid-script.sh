#PBS -N 50-1.0-15-350-2e5 
#PBS -l nodes=2
#PBS -j oe
#PBS -q bulk-a
cd ${PBS_O_WORKDIR}
export OMP_NUM_THREADS=20
aprun -n 4 -N 2 -d ${OMP_NUM_THREADS} -cc depth ./sph.out midfile.bin FDPS-bsc-Rt=50km-Nt=2e5-q=1.0-theta=15-v=350ms-u=0-at-t=1000s > tmp-mid.txt

