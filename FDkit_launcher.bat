@ECHO off

:: 1/ Set path to your python install
set PYTHON_INSTALL=C:\Users\Charles\AppData\Local\Programs\Python\Python39\

:: 2/ Set path to ffmpeg install
set FFMPEG_INSTALL=%cd%\ffmpeg\bin\

:: 3/ Add the two new paths to your PATH environnement variable 
set PATH=%PYTHON_INSTALL%;%FFMPEG_INSTALL%;%PATH%

:: 4/ Run the program
python.exe FDkit.py

pause