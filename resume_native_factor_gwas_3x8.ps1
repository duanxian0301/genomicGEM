$ErrorActionPreference = "Stop"

$jobs = @(
  @{ Tag = "r1"; StartChunk = 4; EndChunk = 32; Cores = 8 },
  @{ Tag = "r2"; StartChunk = 33; EndChunk = 62; Cores = 8 },
  @{ Tag = "r3"; StartChunk = 76; EndChunk = 122; Cores = 8 }
)

foreach ($job in $jobs) {
  $stdout = "/home/shenjing/gs_paths/gwas/step11_factor_gwas_native_results/run_$($job.Tag)_stdout.log"
  $stderr = "/home/shenjing/gs_paths/gwas/step11_factor_gwas_native_results/run_$($job.Tag)_stderr.log"
  $wslCommand = @"
cd /home/shenjing
export START_CHUNK=$($job.StartChunk)
export END_CHUNK=$($job.EndChunk)
export USE_CORES=$($job.Cores)
/usr/bin/Rscript /mnt/d/codex/GenomicSEM/step11_run_factor_gwas_native_wsl.R > $stdout 2> $stderr
"@

  $proc = Start-Process -FilePath "wsl.exe" -ArgumentList @("bash", "-lc", $wslCommand) -WindowStyle Hidden -PassThru
  [PSCustomObject]@{
    Tag = $job.Tag
    Pid = $proc.Id
    StartChunk = $job.StartChunk
    EndChunk = $job.EndChunk
    Cores = $job.Cores
  }
}
