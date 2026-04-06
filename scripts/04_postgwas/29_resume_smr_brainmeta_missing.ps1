param(
  [string[]]$Traits = @("AD", "PD", "LBD", "F1", "F2"),
  [int[]]$Chromosomes = @(1..22),
  [int]$MaxParallel = 4,
  [int]$PollSeconds = 15
)

$python = "python"
$script = "D:\codex\GenomicSEM\scripts\04_postgwas\26_run_smr_single_trait.py"
$outRoot = "D:\文章\GS\postgwas\06_smr_f1f2_ndd\single_trait"
$statusDir = "D:\文章\GS\postgwas\06_smr_f1f2_ndd\batch_status"
$statusFile = Join-Path $statusDir "brainmeta_resume_status.tsv"

New-Item -ItemType Directory -Force -Path $statusDir | Out-Null

$queue = New-Object System.Collections.Generic.List[object]
$launched = New-Object System.Collections.Generic.List[object]
$completed = New-Object System.Collections.Generic.List[object]

foreach ($trait in $Traits) {
  foreach ($chr in $Chromosomes) {
    $outDir = Join-Path $outRoot "$trait\brainmeta_v2_cortex"
    $outPrefix = Join-Path $outDir "${trait}_BrainMeta_chr${chr}"
    $smrFile = "${outPrefix}.smr"
    if (Test-Path -LiteralPath $smrFile) {
      continue
    }
    $queue.Add([pscustomobject]@{
      trait = $trait
      chromosome = $chr
      out_prefix = $outPrefix
    }) | Out-Null
  }
}

Write-Output ("missing_jobs`t{0}" -f $queue.Count)
if ($queue.Count -eq 0) {
  @() | Export-Csv -Path $statusFile -NoTypeInformation -Delimiter "`t" -Encoding UTF8
  Write-Output ("saved_status`t{0}" -f $statusFile)
  exit 0
}

$running = @{}

function Start-NextJob {
  param([object]$Job)
  $args = @(
    $script,
    "--trait", $Job.trait,
    "--panel", "brainmeta",
    "--chrom", "$($Job.chromosome)"
  )
  $proc = Start-Process -FilePath $python -ArgumentList $args -PassThru -WindowStyle Hidden
  $meta = [pscustomobject]@{
    trait = $Job.trait
    chromosome = $Job.chromosome
    pid = $proc.Id
    launch_time = (Get-Date).ToString("s")
    out_prefix = $Job.out_prefix
    status = "running"
  }
  $running[$proc.Id] = $meta
  $launched.Add($meta) | Out-Null
  Write-Output ("launched`t{0}`tchr{1}`tpid={2}" -f $Job.trait, $Job.chromosome, $proc.Id)
}

while ($queue.Count -gt 0 -or $running.Count -gt 0) {
  while ($queue.Count -gt 0 -and $running.Count -lt $MaxParallel) {
    $job = $queue[0]
    $queue.RemoveAt(0)
    Start-NextJob -Job $job
  }

  Start-Sleep -Seconds $PollSeconds

  foreach ($pid in @($running.Keys)) {
    $proc = Get-Process -Id $pid -ErrorAction SilentlyContinue
    if ($null -eq $proc) {
      $meta = $running[$pid]
      $smrFile = "{0}.smr" -f $meta.out_prefix
      $logFile = "{0}.run.log" -f $meta.out_prefix
      $meta.status = if (Test-Path -LiteralPath $smrFile) { "completed" } else { "failed_or_incomplete" }
      $meta.finish_time = (Get-Date).ToString("s")
      $meta.smr_exists = Test-Path -LiteralPath $smrFile
      $meta.log_exists = Test-Path -LiteralPath $logFile
      $completed.Add($meta) | Out-Null
      $running.Remove($pid)
      Write-Output ("finished`t{0}`tchr{1}`t{2}" -f $meta.trait, $meta.chromosome, $meta.status)
    }
  }
}

$completed | Sort-Object trait, chromosome | Export-Csv -Path $statusFile -NoTypeInformation -Delimiter "`t" -Encoding UTF8
Write-Output ("saved_status`t{0}" -f $statusFile)
