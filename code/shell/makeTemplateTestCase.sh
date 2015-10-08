# $1: target for CLEANING

make PARAMS=allParameters.template CLEAN=1 $1
make PARAMS=allParameters.template ALL | tee logFiles/makeLog_HP_$1.txt
make PARAMS=allParameters.template plotIndels | tee logFiles/makeLog_HP_Indels.txt
