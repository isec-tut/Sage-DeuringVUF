# Usage: .\scripts\sage.ps1 examples/eval.py [arguments...]
# Use $args directly so Python options such as --trials pass through unchanged.
$ErrorActionPreference = 'Stop'
if ($args.Count -eq 0) {
    Write-Host 'Usage: .\scripts\sage.ps1 <project-relative-file.py> [arguments...]'
    exit 1
}

$sageImage = 'sagemath/sagemath@sha256:19995db6194f4a4bab18ce9a88556fd15b9ed5e916b4504fefe618a7796ddbdb'
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$relativeFile = ([string]$args[0]).Replace('\', '/')
if ($relativeFile.StartsWith('./')) {
    $relativeFile = $relativeFile.Substring(2)
}
if ($relativeFile -cnotmatch '^[A-Za-z_][A-Za-z_0-9]*(/[A-Za-z_][A-Za-z_0-9]*)*\.py$') {
    throw 'Specify a Python file relative to the project root, e.g. examples/eval.py.'
}
if (-not (Test-Path -LiteralPath (Join-Path $projectRoot $relativeFile) -PathType Leaf)) {
    throw "Python file does not exist: $relativeFile"
}
$module = $relativeFile.Substring(0, $relativeFile.Length - 3).Replace('/', '.')
$pythonArgs = @()
if ($args.Count -gt 1) {
    $pythonArgs = @($args[1..($args.Count - 1)])
}

& docker run --rm --mount "type=bind,source=$projectRoot,target=/work,readonly" -w /work $sageImage sage -python -B -u -m $module @pythonArgs
exit $LASTEXITCODE
