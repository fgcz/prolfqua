we did implement lm with imputation see @DONE_LM_Impute.md 
based on build_model_impute and the default Contrast R 6 class.

I would like to extend the limma in the same fashion, although here:

Two options:

Option 1:
a) use data with NA fit limma,
b) use data with LOD fit limma, 
determine failed models in a), and replace those with the fits from b).


Option 2:
