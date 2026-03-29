https://github.com/SmythLab/limpa
https://bioconductor.org/packages/release/bioc/html/limpa.html
also sorce code /Users/wolski/projects/prolfqua_fml/limpa

We also have both limpa papers in  /Users/wolski/projects/prolfqua_fml/addinfo

and at: ./.cursor/skills/adding-models-to-prolfqua/SKILL.md

- limpa does missing value imputation starting from precursor data, but only on protein level (?), or will it also do imputation, on precursor or peptide level?

- going from precursor to protein would possibly -> encourages to be one of the Aggregator.

- if precursor to precursor -> what is it then? Imputation?

Question relevant because we want both to support Peptide and protein level analysis.

- is imputation a separate step then protein abundance inference in limpa?
- is protein inference spearated from the modelling?

Can I use impuation and then fit a different then the limpa model to the data, e.g. limma or prolfqua lm?

If fitting limpa model can I disable moderation and apple DeqMS moderation?


Once inferred missigness is limpa any different then limma?

I am ready to introduce dependency on limpa similarily to the dependency on limma, The internal datasetructures are similar in limma and limpa so we might be reusing some code which we used for limma.