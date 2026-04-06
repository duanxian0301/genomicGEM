param(
  [string[]]$Traits = @("AD", "PD", "LBD", "F1", "F2"),
  [int[]]$Chromosomes = @(1..22),
  [switch]$RunInBackground = $true
)

$python = "python"
$script = "D:\codex\GenomicSEM\scripts\04_postgwas\26_run_smr_single_trait.py"
$statusDir = "D:\文章\GS\postgwas\06_smr_f1f2_ndd\batch_status"
$statusFile = Join-Path $statusDir "brainmeta_batch_launch.tsv"

New-Item -ItemType Directory -Force -Path $statusDir | Out-Null

$rows = @()

foreach ($trait in $Traits) {
  foreach ($chr in $Chromosomes) {
    $args = @(
      $script,
      "--trait", $trait,
      "--panel", "brainmeta",
      "--chrom", "$chr"
    )

    if ($RunInBackground) {
      $proc = Start-Process -FilePath $python -ArgumentList $args -PassThru -WindowStyle Hidden
      $procId = $proc.Id
      $mode = "background"
    } else {
      & $python @args
      $procId = $PID
      $mode = "foreground"
    }

    $rows += [pscustomobject]@{
      trait = $trait
      chromosome = $chr
      mode = $mode
      pid = $procId
      launch_time = (Get-Date).ToString("s")
    }
  }
}

$rows | Export-Csv -Path $statusFile -NoTypeInformation -Delimiter "`t" -Encoding UTF8
Write-Output "saved_status`t$statusFile"
