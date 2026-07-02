@echo off
if "%1"=="h" goto start

:: Generate a temporary VBScript to execute run.bat in hidden mode safely with spaces in path
echo Set WshShell = CreateObject("WScript.Shell") > "%temp%\launch_promoter.vbs"
echo WshShell.Run "cmd.exe /c """ ^& "%~dpnx0" ^& """ h", 0, False >> "%temp%\launch_promoter.vbs"
wscript.exe "%temp%\launch_promoter.vbs"
del "%temp%\launch_promoter.vbs"
exit

:start
title Yeast Promoter Modifier Server Launcher

:: 1. Check Python installation
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo MsgBox "Python is not installed or not in PATH!", 16, "Error" > "%temp%\alert.vbs"
    wscript.exe "%temp%\alert.vbs"
    del "%temp%\alert.vbs"
    exit /b 1
)

:: 2. Verify required libraries and install if missing
python -c "import fastapi, uvicorn, requests, numpy, torch" >nul 2>&1
if %errorlevel% neq 0 (
    python -m pip install fastapi uvicorn requests numpy >nul 2>&1
    python -m pip install torch --index-url https://download.pytorch.org/whl/cpu >nul 2>&1
    if %errorlevel% neq 0 (
        echo MsgBox "Library installation failed! Check network and try again.", 16, "Error" > "%temp%\alert.vbs"
        wscript.exe "%temp%\alert.vbs"
        del "%temp%\alert.vbs"
        exit /b 1
    )
)

:: 3. Run Self-Diagnostic Tests before starting server
python "%~dp0test_prediction.py" > "%~dp0error_log.txt" 2>&1
if %errorlevel% neq 0 (
    echo MsgBox "Base prediction model diagnostic test failed! Check error_log.txt for details.", 16, "Error" > "%temp%\alert.vbs"
    wscript.exe "%temp%\alert.vbs"
    del "%temp%\alert.vbs"
    exit /b 1
)

python "%~dp0main_deep.py" --fast-test >> "%~dp0error_log.txt" 2>&1
if %errorlevel% neq 0 (
    echo MsgBox "Deep learning engine diagnostic test failed! Check error_log.txt for details.", 16, "Error" > "%temp%\alert.vbs"
    wscript.exe "%temp%\alert.vbs"
    del "%temp%\alert.vbs"
    exit /b 1
)

:: Clear error log on success
if exist "%~dp0error_log.txt" del "%~dp0error_log.txt" >nul 2>&1

:: 4. Clear port 8080 zombie process (Only LISTENING sockets to avoid side effects)
for /f "tokens=5" %%a in ('netstat -aon ^| findstr :8080 ^| findstr LISTENING') do (
    taskkill /F /PID %%a >nul 2>&1
)

:: 5. Launch Google Chrome Browser (Detect standard paths, fallback to default browser if missing)
if exist "C:\Program Files\Google\Chrome\Application\chrome.exe" (
    start "" "C:\Program Files\Google\Chrome\Application\chrome.exe" http://127.0.0.1:8080/
) else if exist "C:\Program Files (x86)\Google\Chrome\Application\chrome.exe" (
    start "" "C:\Program Files (x86)\Google\Chrome\Application\chrome.exe" http://127.0.0.1:8080/
) else if exist "%LocalAppData%\Google\Chrome\Application\chrome.exe" (
    start "" "%LocalAppData%\Google\Chrome\Application\chrome.exe" http://127.0.0.1:8080/
) else (
    start http://127.0.0.1:8080/
)

:: 6. Launch FastAPI server
python "%~dp0main.py"
