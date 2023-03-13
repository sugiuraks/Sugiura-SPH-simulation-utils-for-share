###################################################
# timestamp: 2022/11/16 13:30 sugiura
# multiple-submit.shから呼び出されるスクリプトです.
# ただしmidfile.bin(前回の計算の最後に書き出される中間条件ファイル)
# から計算が始まるようになっており, 前回の計算の続きを担当するスクリプトです.
# 読み込む衝突条件ファイルがmidfile.binになっている以外は
# FDPS-basalt-sphere-collision-multiple.shと同じです.
###################################################
#PBS -N test
#PBS -l nodes=1
#PBS -q bulk-a
#PBS -j oe
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=20
aprun -n 2 -N 2 -d ${OMP_NUM_THREADS} -cc depth ./sph.out midfile.bin FDPS-bsc-Rt=50km-Nt=5e4-q=0.75-theta=${ARG1}-v=${ARG2}ms > tmp${PBS_JOBID%.*}.txt
