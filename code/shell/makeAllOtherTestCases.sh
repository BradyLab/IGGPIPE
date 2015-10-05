# $1: target for CLEANING

make PARAMS=allParameters/allParameters.AtCL CLEAN=1 $1
make PARAMS=allParameters/allParameters.AtCL ALL | tee logFiles/makeLog_CL_$1.txt
make PARAMS=allParameters/allParameters.AtCL plotInDels | tee logFiles/makeLog_CL_Indels.txt

make PARAMS=allParameters/allParameters.HPT CLEAN=1 $1
make PARAMS=allParameters/allParameters.HPT ALL | tee logFiles/makeLog_HPT_$1.txt
make PARAMS=allParameters/allParameters.HPT plotInDels | tee logFiles/makeLog_HPT_Indels.txt

make PARAMS=allParameters/allParameters.HPT2.3 CLEAN=1 $1
make PARAMS=allParameters/allParameters.HPT2.3 ALL | tee logFiles/makeLog_HPT2.3_$1.txt
make PARAMS=allParameters/allParameters.HPT2.3 plotInDels | tee logFiles/makeLog_HPT2.3_Indels.txt

make PARAMS=allParameters/allParameters.Ta3 CLEAN=1 $1
make PARAMS=allParameters/allParameters.Ta3 ALL | tee logFiles/makeLog_Ta3_$1.txt
make PARAMS=allParameters/allParameters.Ta3 plotInDels | tee logFiles/makeLog_Ta3_Indels.txt
