# Fix Limma Sample/Annotation Alignment

## Plan

1. Fix `tidy_to_wide_config()` so sample columns in `data` are explicitly ordered to match rows in `annotation`.
2. Add a `data_wide()` regression test with intentionally interleaved sample order.
3. Add a limma-impute regression test showing reported `diff` matches observed treatment-minus-control group means after sample order shuffling.
4. Run focused `LFQData` and `ContrastsLimma` tests.

## Root Cause

`tidy_to_wide_config()` sorted the sample annotation table independently from the pivoted wide data columns. Limma used `wide$data` and `wide$annotation` together, so mismatched order assigned the wrong group labels to matrix columns.
