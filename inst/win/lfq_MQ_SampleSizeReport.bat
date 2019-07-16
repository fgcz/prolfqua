@echo off
:: start help with "lfq_MQ_SampleSizeReport --help"

SET args=%*
  Rscript %~dp0/../run_scripts/lfq_MQ_SampleSizeReport.R %args%
