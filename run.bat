@echo off
setlocal

set "ROOT_DIR=%~dp0"
set "ANALYSIS_SCRIPT=%ROOT_DIR%tools\run_coronary_analysis.py"
set "PYTHON_CMD="
set "PYTHON_DISPLAY="

if not exist "%ANALYSIS_SCRIPT%" (
  echo [ERROR] Analysis entry not found: %ANALYSIS_SCRIPT%
  exit /b 1
)

if exist "%ROOT_DIR%.venv\Scripts\python.exe" (
  set "PYTHON_CMD="%ROOT_DIR%.venv\Scripts\python.exe""
  set "PYTHON_DISPLAY=%ROOT_DIR%.venv\Scripts\python.exe"
) else if exist "%ROOT_DIR%venv\Scripts\python.exe" (
  set "PYTHON_CMD="%ROOT_DIR%venv\Scripts\python.exe""
  set "PYTHON_DISPLAY=%ROOT_DIR%venv\Scripts\python.exe"
) else (
  where python >nul 2>nul
  if not errorlevel 1 (
    set "PYTHON_CMD=python"
    for /f "usebackq delims=" %%I in (`where python`) do (
      set "PYTHON_DISPLAY=%%I"
      goto python_ready
    )
  )

  where py >nul 2>nul
  if not errorlevel 1 (
    set "PYTHON_CMD=py -3"
    set "PYTHON_DISPLAY=py -3"
    goto python_ready
  )

  echo [ERROR] No usable Python interpreter found.
  echo [ERROR] Please create .venv\Scripts\python.exe or ensure python / py is available in PATH.
  exit /b 1
)

:python_ready
set "COMMAND=%PYTHON_CMD% "%ANALYSIS_SCRIPT%""

if /I "%~1"=="--dry-run" (
  echo %COMMAND%
  exit /b 0
)

echo [INFO] Running coronary analysis...
echo [INFO] Workspace: %ROOT_DIR%
echo [INFO] Python: %PYTHON_DISPLAY%
echo [INFO] Command: %COMMAND%
echo.

call %COMMAND%
set "EXIT_CODE=%ERRORLEVEL%"

if /I "%~1"=="--no-pause" (
  exit /b %EXIT_CODE%
)

echo.
if "%EXIT_CODE%"=="0" (
  echo [INFO] Analysis finished successfully.
) else (
  echo [ERROR] Analysis failed with exit code %EXIT_CODE%.
)
pause
exit /b %EXIT_CODE%
