param(
    [Parameter(Mandatory = $true)][string]$SourceDir,
    [Parameter(Mandatory = $true)][string]$IncludeDir,
    [Parameter(Mandatory = $true)][string]$OutputLibrary,
    [Parameter(Mandatory = $true)][string]$ObjectDir
)

$ErrorActionPreference = 'Stop'

if (-not (Get-Command cl.exe -ErrorAction SilentlyContinue)) {
    throw 'cl.exe is not on PATH. Start MATLAB from an x64 Native Tools Command Prompt for Visual Studio.'
}
if (-not (Get-Command lib.exe -ErrorAction SilentlyContinue)) {
    throw 'lib.exe is not on PATH. Start MATLAB from an x64 Native Tools Command Prompt for Visual Studio.'
}

New-Item -ItemType Directory -Force -Path $ObjectDir | Out-Null
New-Item -ItemType Directory -Force -Path (Split-Path $OutputLibrary) | Out-Null

$objects = @()
Get-ChildItem -Path $SourceDir -Filter '*.c' | ForEach-Object {
    $object = Join-Path $ObjectDir ($_.BaseName + '.obj')
    $objects += $object
    if ((-not (Test-Path $object)) -or ($_.LastWriteTimeUtc -gt (Get-Item $object).LastWriteTimeUtc)) {
        Write-Host "Compiling CSPICE: $($_.Name)"
        & cl.exe /nologo /c /O2 /MD /D NON_UNIX_STDIO "/I$IncludeDir" "/Fo$object" $_.FullName
        if ($LASTEXITCODE -ne 0) { throw "CSPICE compilation failed for $($_.FullName)" }
    }
}

$responseFile = Join-Path $ObjectDir 'cspice_objects.rsp'
$objects | ForEach-Object { '"' + $_ + '"' } | Set-Content -Encoding ASCII $responseFile
if (Test-Path $OutputLibrary) { Remove-Item -Force $OutputLibrary }
& lib.exe /nologo "/OUT:$OutputLibrary" "@$responseFile"
if ($LASTEXITCODE -ne 0) { throw 'Unable to create the CSPICE static library.' }
