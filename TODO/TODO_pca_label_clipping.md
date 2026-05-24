# TODO: Fix PCA sample label clipping

## Problem

`plot_pca()` currently draws sample labels with `geom_text()` plus a fixed nudge. In small sample sets, labels at the
panel edge are clipped or overlap because the plot has no extra expansion and the text geometry does not repel labels
from each other or from the panel boundary.

This is visible in downstream `prolfquapp` Quarto reports for SAINT workflows where sample labels such as
`notag_control_*` and `PPE4_*` are partially cut off.

## Plan

- Change `plot_pca()` in `R/tidyMS_plotting.R` to draw PCA labels with `ggrepel::geom_text_repel()` when
  `add_txt = TRUE`.
- Keep the existing public arguments (`add_txt`, `nudge`) and preserve current color/shape mappings.
- Add plot expansion and disable coordinate clipping so labels have room to render.
- Add a focused test that `plot_pca(add_txt = TRUE)` uses a repel text layer and has clipping disabled.
- Run targeted plotting tests and format touched R files.
