$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root

if (-not (Test-Path ".api_pid")) {
    Write-Host "No saved DrugReflector UI process was found."
    exit 0
}

$serverPid = Get-Content ".api_pid"
if ($serverPid) {
    try {
        Stop-Process -Id ([int]$serverPid) -Force -ErrorAction Stop
        Write-Host "Stopped DrugReflector UI process $serverPid"
    } catch {
        Write-Host "Process $serverPid was not running."
    }
}

Remove-Item ".api_pid" -Force -ErrorAction SilentlyContinue
Remove-Item ".api_port" -Force -ErrorAction SilentlyContinue
