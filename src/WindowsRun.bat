REM Check for a user defined environnment variable (PYTHON4SSNAKE) containing the location of python interpreter.
REM Else this script will start ssNake only if pyw.exe or pythonw.exe is found in current path search. 
REM Otherwise it will quit
@echo off
set mypath=%cd%
cd ..\..

IF DEFINED PYTHON4SSNAKE (
set PYTHONW=%PYTHON4SSNAKE%
goto run_ssnake
)

where /q pyw
IF ERRORLEVEL 1 (
    ECHO 'pyw' not in search path, try 'pythonw'
) ELSE (
    SET PYTHONW=pyw
    goto run_ssnake
)

where /q pythonw
IF ERRORLEVEL 1 (
    ECHO 'pythonw' also not in search path.
) ELSE (
    SET PYTHONW=pythonw
    goto run_ssnake
)
ECHO No python interpreter found: about to leave
PAUSE
EXIT

:run_ssnake
REM Let's try running ssNake !!!!
start "ssNake" %PYTHONW% "%mypath%"\ssNake.py

