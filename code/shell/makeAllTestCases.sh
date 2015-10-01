make PARAMS=allParameters.test CLEAN=1 ALL
make PARAMS=allParameters.test ALL | tee logFiles/makeLog_Test.txt
make PARAMS=allParameters.test plotInDels | tee logFiles/makeLog_Test_Indels.txt

make PARAMS=allParameters.template CLEAN=1 ALL
make PARAMS=allParameters.template ALL | tee logFiles/makeLog_HP.txt
make PARAMS=allParameters.template plotInDels | tee logFiles/makeLog_HP_Indels.txt

make PARAMS=allParameters/allParameters.AtCL CLEAN=1 ALL
make PARAMS=allParameters/allParameters.AtCL ALL | tee logFiles/makeLog_CL.txt
make PARAMS=allParameters/allParameters.AtCL plotInDels | tee logFiles/makeLog_CL_Indels.txt

make PARAMS=allParameters/allParameters.HPT CLEAN=1 ALL
make PARAMS=allParameters/allParameters.HPT ALL | tee logFiles/makeLog_HPT.txt
make PARAMS=allParameters/allParameters.HPT plotInDels | tee logFiles/makeLog_HPT_Indels.txt

make PARAMS=allParameters/allParameters.HPT2.3 CLEAN=1 ALL
make PARAMS=allParameters/allParameters.HPT2.3 ALL | tee logFiles/makeLog_HPT2.3.txt
make PARAMS=allParameters/allParameters.HPT2.3 plotInDels | tee logFiles/makeLog_HPT2.3_Indels.txt

make PARAMS=allParameters/allParameters.Ta3 CLEAN=1 ALL
make PARAMS=allParameters/allParameters.Ta3 ALL | tee logFiles/makeLog_Ta3.txt
make PARAMS=allParameters/allParameters.Ta3 plotInDels | tee logFiles/makeLog_Ta3_Indels.txt

make PARAMS=allParameters/allParameters.Hs.Soly CLEAN=1 ALL
make PARAMS=allParameters/allParameters.Hs.Soly ALL | tee logFiles/makeLog_Hs.Soly.txt
make PARAMS=allParameters/allParameters.Hs.Soly plotInDels | tee logFiles/makeLog_Hs.Soly_Indels.txt
