@echo off
set MATLAB=C:\mat2010a
set MSVCDir=c:\program files (x86)\microsoft visual studio 9.0\VC

"C:\mat2010a\rtw\bin\win32\envcheck" INCLUDE "c:\program files (x86)\microsoft visual studio 9.0\VC\include"
if errorlevel 1 goto vcvars32
"C:\mat2010a\rtw\bin\win32\envcheck" PATH "c:\program files (x86)\microsoft visual studio 9.0\VC\bin"
if errorlevel 1 goto vcvars32
goto make

:vcvars32
set VSINSTALLDIR=c:\program files (x86)\microsoft visual studio 9.0
set VCINSTALLDIR=c:\program files (x86)\microsoft visual studio 9.0\VC
set FrameworkSDKDir=c:\program files (x86)\microsoft visual studio 9.0\SDK\v3.5
call "C:\mat2010a\toolbox\rtw\rtw\private\vcvars32_900.bat"

:make
cd .
nmake -f EssexAidProcessVFrameSwitchable_rtw.mk  GENERATE_REPORT=0 ADD_MDL_NAME_TO_GLOBALS=0
@if errorlevel 1 goto error_exit
exit /B 0

:error_exit
echo The make command returned an error of %errorlevel%
An_error_occurred_during_the_call_to_make
