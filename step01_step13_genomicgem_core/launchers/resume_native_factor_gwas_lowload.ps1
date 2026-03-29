$ErrorActionPreference = "Stop"

$wslCommand = @'
cd /home/shenjing
export START_CHUNK=1
export END_CHUNK=122
export USE_CORES=4
/usr/bin/Rscript /mnt/d/codex/GenomicSEM/step11_run_factor_gwas_native_wsl.R > /home/shenjing/gs_paths/gwas/step11_factor_gwas_native_results/run_resume_lowload_stdout.log 2> /home/shenjing/gs_paths/gwas/step11_factor_gwas_native_results/run_resume_lowload_stderr.log
'@

$proc = Start-Process -FilePath "wsl.exe" -ArgumentList @("bash","-lc",$wslCommand) -WindowStyle Hidden -PassThru
$proc.Id
