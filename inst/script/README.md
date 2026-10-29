# Data in prolfqua

Most examples, tests and vignettes use simulated data (`sim_lfq_data_peptide_config()`, `sim_lfq_data_protein_config()`, `sim_lfq_data_2factor_config()`), created at run time. The package ships four small datasets in `data/`; this file says where they come from.

## data_ionstar

- Content: peptide intensities of 200 proteins, 20 samples (dilutions a to e of an E. coli spike-in into a human cell lysate, 4 replicates each), 26520 rows
- Source: IonStar benchmark dataset, PRIDE [PXD003881](https://www.ebi.ac.uk/pride/archive/projects/PXD003881) (Shen et al., J. Proteome Res. 2018)
- Processing: raw files quantified with MaxQuant; peptide intensities of a subset of 200 proteins kept
- Object: R6 object of class `IonstarData` with `$data`, `$config` (an `AnalysisConfiguration`) and the methods `filtered()`, `normalized()` and `subset_normalized()`
- Used by: examples in prolfqua, prolfquapp (`transform_lfqdata()`, `plot_abundance_vs_percent()`, QC vignette)

## data_checksummarizationrobust87, data_checksummarizerobust69

- Content: log2 peptide areas (`log2Area`, `peptide_Id`, `sampleName`) of a single protein, 20 rows each
- Source: one protein each, taken from an FGCZ label-free quantification dataset, as reference input for the robust summarisation (median polish, rlm)
- Used by: `tests/testthat/test-tidyMS_aggregation.R`

## x5463yzwer453bbb

- Content: sample annotation table of 36 samples with the columns `Relative Path`, `Name`, `Group`, `factor_A` (MI / MINOCA) and `factor_B` (sampling time T0 / T150 / T300)
- Source: sample annotation of an FGCZ proteomics order (raw files `20250429_*`); only file names and the two factors are kept
- Used by: examples of `annotation_add_contrasts()` and the article `ContrastsFromFactors.Rmd`
