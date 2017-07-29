:DLEARN (STANDARD) RUN PROCEDURE
echo off
if exist input.dat del input.dat
if exist output.dat del output.dat
copy %1.dat input.dat
dlearn 
if exist %1.out del %1.out
if exist input.dat del input.dat
if exist output.dat rename output.dat %1.out
if exist rsinput.dat del rsinput.dat
if exist rsoutput.dat del rsoutput.dat
echo ********* end of dlearn (standard) *************************************
