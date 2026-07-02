@echo off
title Yeast Promoter Modifier Server Launcher
color 0B

echo ===================================================
echo     Yeast Promoter Modifier System Launcher
echo ===================================================
echo.

:: 1. Check Python installation
echo [+] Checking Python environment...
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python is not installed or not in PATH!
    echo Please install Python and try again.
    pause
    exit /b 1
)

:: 2. Verify required libraries and install if missing
echo [+] Verifying library dependencies (fastapi, uvicorn, requests, numpy, torch)...
python -c "import fastapi, uvicorn, requests, numpy, torch" >nul 2>&1
if %errorlevel% neq 0 (
    echo [!] Missing required dependencies. Installing fastapi, uvicorn, requests, numpy, torch...
    python -m pip install fastapi uvicorn requests numpy
    python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
    if %errorlevel% neq 0 (
        echo [ERROR] Library installation failed! Please check your network.
        pause
        exit /b 1
    )
)
echo [+] Dependencies verified successfully.

:: 3. Run Self-Diagnostic Tests before starting server
echo [+] Running base prediction model diagnostic tests...
python test_prediction.py >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Base prediction model diagnostic test failed!
    python test_prediction.py
    pause
    exit /b 1
)

echo [+] Running deep learning & GA model diagnostic tests...
python main_deep.py >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Deep learning engine diagnostic test failed!
    python main_deep.py
    pause
    exit /b 1
)
echo [+] All self-diagnostic tests passed successfully!

:: 4. Clear port 8080 zombie process (Only LISTENING sockets to avoid side effects)
echo [+] Checking port 8080 status...
for /f "tokens=5" %%a in ('netstat -aon ^| findstr :8080 ^| findstr LISTENING') do (
    echo [!] Port 8080 is occupied. Terminating zombie process PID %%a...
    taskkill /F /PID %%a >nul 2>&1
)

:: 5. Launch Google Chrome Browser
echo [+] Launching Google Chrome browser at http://127.0.0.1:8080/ ...
start chrome http://127.0.0.1:8080/

:: 6. Launch FastAPI server
echo [+] Starting FastAPI server...
echo ---------------------------------------------------
python main.py
pause
