$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root

Write-Host "Starting DrugReflector UI from $root"

$needsInstall = $false
try {
    python -c "import fastapi, uvicorn, drugreflector" | Out-Null
} catch {
    $needsInstall = $true
}

if ($needsInstall) {
    Write-Host "Installing required dependencies..."
    python -m pip install -e ".[api]"
}

$frontendRoot = Join-Path $root "frontend"
$frontendDist = Join-Path $frontendRoot "dist\\index.html"
if ((Test-Path (Join-Path $frontendRoot "package.json")) -and -not (Test-Path $frontendDist)) {
    $npm = Get-Command npm -ErrorAction SilentlyContinue
    if ($npm) {
        Write-Host "Building React UI..."
        Push-Location $frontendRoot
        try {
            if (-not (Test-Path "node_modules")) {
                npm ci
            }
            npm run build
        } finally {
            Pop-Location
        }
    } else {
        Write-Host "npm was not found, so the script will fall back to the packaged static UI."
    }
}

if (Test-Path ".api_pid") {
    $existingPid = Get-Content ".api_pid"
    if ($existingPid) {
        try {
            Stop-Process -Id ([int]$existingPid) -Force -ErrorAction Stop
        } catch {
        }
    }
    Remove-Item ".api_pid" -Force -ErrorAction SilentlyContinue
}

$preferredPorts = @(8000, 8010, 8765)
$port = $null
foreach ($candidatePort in $preferredPorts) {
    $listener = Get-NetTCPConnection -LocalPort $candidatePort -State Listen -ErrorAction SilentlyContinue
    if (-not $listener) {
        $port = $candidatePort
        break
    }
}

if (-not $port) {
    throw "No free port was found in $($preferredPorts -join ', ')."
}

$server = Start-Process python -ArgumentList "-m","uvicorn","drugreflector.api:app","--host","127.0.0.1","--port",$port -PassThru
Set-Content -Path ".api_pid" -Value $server.Id
Set-Content -Path ".api_port" -Value $port

$ready = $false
for ($i = 0; $i -lt 30; $i++) {
    Start-Sleep -Seconds 1
    try {
        $health = Invoke-WebRequest -UseBasicParsing "http://127.0.0.1:$port/api/checkpoints"
        if ($health.StatusCode -eq 200) {
            $ready = $true
            break
        }
    } catch {
    }
}

if (-not $ready) {
    throw "DrugReflector UI did not become ready on http://127.0.0.1:$port/"
}

Start-Process "http://127.0.0.1:$port/"
Write-Host "DrugReflector UI is ready at http://127.0.0.1:$port/"
Write-Host "To stop it later, run .\stop_ui.ps1"
