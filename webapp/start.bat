@echo off
rem GOAI Explorer launcher: double-click to start, browser opens http://127.0.0.1:8765
rem Keep this window open while using the app; closing it stops the server.
cd /d "%~dp0.."
set PYTHONIOENCODING=utf-8

rem If server is already running, just open the page
curl -s -o nul http://127.0.0.1:8765/api/state && (
  echo Server already running, opening browser...
  start "" http://127.0.0.1:8765
  pause
  exit /b
)

echo Starting GOAI Explorer (keep this window open)...
start /b python -u webapp/server.py
timeout /t 3 /nobreak >nul
start "" http://127.0.0.1:8765
echo.
echo Server running: http://127.0.0.1:8765  (close this window to stop)
pause >nul
