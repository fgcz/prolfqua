# LFQService - for Label Free Quantification Services


The R package contains vignettes and functions for analyzing mass spec device based LFQ proteomics experiments at the [FGCZ](http://www.fgcz.ch/).


## 1. System Requirements  

a Windows/Linux/MacOSX x64 platform R 3.5.


## 2. Running R-scripts

Generate QC report from bat file.
First add `<LFQService_path>/win` to your path variable. Then you can generate a QC report from a maxquant QC by running.


```
lfq_MQ_SampleSizeReport.bat .\data\1296877_QC.zip
```
