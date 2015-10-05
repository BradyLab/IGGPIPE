# $1: target for CLEANING

make PARAMS=allParameters.test CLEAN=1 $1
make PARAMS=allParameters.test ALL | tee logFiles/makeLog_Test_$1.txt
make PARAMS=allParameters.test plotInDels | tee logFiles/makeLog_Test_Indels.txt
