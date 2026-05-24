# Volcano Model Color Mapping

## Goal

Use a fast label-order hack so mixed imputed/non-imputed model results sort with the non-imputed model first.

The immediate observed case is:

- `WaldTest_moderated` should be black.
- `WaldTest_moderated_imputed` should not be black.

## Root Cause

The default volcano color scale in `prolfqua` uses an unnamed manual color vector such as `c("black", "green")`. ggplot assigns those colors by factor/discrete level order. Because `WaldTest_imputed_moderated` sorts before `WaldTest_moderated`, the imputed level can receive black and the non-imputed level receives green.

## Plan

1. Change the moderated contrast label suffixing so model names ending in `_imputed` become `_moderated_imputed` instead of `_imputed_moderated`.
2. Update tests and comments that assert or document the old `WaldTest_imputed_moderated` label.
3. Run the focused `prolfqua` imputation test file.

## Notes

This is intentionally a fast brittle hack. A named palette would be the cleaner long-term fix because it would not depend on string sorting.
