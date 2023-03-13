###################################################
# timestamp: 2022/11/16 13:30 sugiura
# multiple-submit.shから呼び出されるスクリプトです.
# 元々のmultiple-submit.shはこのスクリプトを呼び出す時に,
# 第一引数に衝突角度を, 第二引数に衝突速度を指定するようにしているので,
# それぞれ適切な箇所に引数を展開して適切な初期条件ファイルが使用されるようにしてください.
#
# #PBSで始まる行などのオプションなどの意味については天文台XC50のマニュアルなどを参照してください.
#
# aprunから始まる行が実際にジョブを投げる行です.
# オプションを除くと,
# 第一引数(./sph.out): 走らせる実行ファイルの名前
# 第二引数: 走らせる初期条件ファイル.binの名前
# 第三引数: 計算で書き出されるファイルなどが保存されるディレクトリ名
###################################################
#PBS -N test
#PBS -l nodes=1
#PBS -q bulk-a
#PBS -j oe
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=20
aprun -n 2 -N 2 -d ${OMP_NUM_THREADS} -cc depth ./sph.out FDPS-bsc-Rt=50km-Nt=5e4-q=0.75-theta=${ARG1}-v=${ARG2}ms.bin FDPS-bsc-Rt=50km-Nt=5e4-q=0.75-theta=${ARG1}-v=${ARG2}ms > tmp${PBS_JOBID%.*}.txt
