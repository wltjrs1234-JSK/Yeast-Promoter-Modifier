@echo off
if "%1"=="h" goto start

:: Self-backgrounding wrapper: re-run myself with 'h' argument in hidden mode (SW_HIDE = 0)
mshta vbscript:CreateObject("Wscript.Shell").Run("cmd.exe /c ""%~dpnx0"" h",0)(window.close)&&exit

:start
title Yeast Promoter Modifier Server Launcher

:: 1. Check Python installation
python --version >nul 2>&1
if %errorlevel% neq 0 (
    mshta vbscript:MsgBox("Python is not installed or not in PATH!", 16, "Error") (window.close)
    exit /b 1
)

:: 2. Verify required libraries and install if missing
python -c "import fastapi, uvicorn, requests, numpy, torch" >nul 2>&1
if %errorlevel% neq 0 (
    python -m pip install fastapi uvicorn requests numpy >nul 2>&1
    python -m pip install torch --index-url https://download.pytorch.org/whl/cpu >nul 2>&1
    if %errorlevel% neq 0 (
        mshta vbscript:MsgBox("Library installation failed! Please check your network and try again.", 16, "Error") (window.close)
        exit /b 1
    )
)

:: 3. Run Self-Diagnostic Tests before starting server
python test_prediction.py > error_log.txt 2>&1
if %errorlevel% neq 0 (
    mshta vbscript:MsgBox("Base prediction model diagnostic test failed! Check error_log.txt for details.", 16, "Error") (window.close)
    exit /b 1
)

python main_deep.py >> error_log.txt 2>&1
if %errorlevel% neq 0 (
    mshta vbscript:MsgBox("Deep learning engine diagnostic test failed! Check error_log.txt for details.", 16, "Error") (window.close)
    exit /b 1
)

:: Clear error log on success
del error_log.txt >nul 2>&1

:: 4. Clear port 8080 zombie process (Only LISTENING sockets to avoid side effects)
for /f "tokens=5" %%a in ('netstat -aon ^| findstr :8080 ^| findstr LISTENING') do (
    taskkill /F /PID %%a >nul 2>&1
)

:: 5. Launch Google Chrome Browser
start chrome http://127.0.0.1:8080/

:: 6. Launch FastAPI server
python main.py
