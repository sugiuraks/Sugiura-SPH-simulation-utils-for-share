#!/bin/bash
###################################################
# timestamp: 2022/11/16 13:30 sugiura
# グリッド上の衝突速度・衝突角度を持つ多数の衝突初期条件ファイルたちの数値計算を
# このスクリプトを使うだけで全て投げることができるスクリプトです.
# (スクリプトをうまく改造すれば衝突速度・角度以外のパラメータも自動で変えながら投げられるようになると思います)
# 衝突速度・衝突角度のグリッドの初期値・終了値・幅はスクリプトの頭のInit~・Fin~・Del~で指定します.
# またNumForSimで1つのパラメータの数値計算のジョブを連続して投げる回数-1を指定できます.
# (注: 天文台のXC50は一般的に12時間で計算が終了させられるので, それまでに終了して続きの計算を投げなければいけません)
# (その続きの計算を投げる回数-1をNumForSimで指定する, という意味です)
#
# 使い方:
# ./にあるFDPS-basalt-sphere-collision-multiple.shとmid-script-multiple.shを適切に書き換えた上で,
# $ ./multiple-submit.sh
###################################################

NumForSim=2 #number of jobs for each simulation - 1
j=0
flag=0

InitVImp=50
InitThetaImp=5

FinVImp=600
FinThetaImp=5

DelVImp=50
DelThetaImp=5

for ((VImp=$InitVImp ; VImp<=$FinVImp ; VImp=$[$VImp+$DelVImp]))
do
    for ((ThetaImp=$InitThetaImp ; ThetaImp<=$FinThetaImp ; ThetaImp=$[$ThetaImp+$DelThetaImp]))
    do
	if [ $flag -eq 0 ]
	then
	    qsub -v ARG1=${ThetaImp},ARG2=${VImp} ./FDPS-basalt-sphere-collision-multiple.sh > job_id.txt
	    sleep 0.1s
	    flag=1
	else
	    JobIDSTR=$(<./job_id.txt)
	    JobID=${JobIDSTR%.*}
	    echo $JobID
            qsub -W depend=afterok:$JobID -v ARG1=${ThetaImp},ARG2=${VImp} ./FDPS-basalt-sphere-collision-multiple.sh > job_id.txt
	    sleep 0.1s
	fi	
	for ((i=0 ; i<$NumForSim ; i++))
	do
	    JobIDSTR=$(<./job_id.txt)
            JobID=${JobIDSTR%.*}
	    echo $JobID
	    qsub -W depend=afterok:$JobID -v ARG1=${ThetaImp},ARG2=${VImp} ./mid-script-multiple.sh > job_id.txt
	    sleep 0.1s
	done
    done
done
