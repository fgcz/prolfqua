# Volcano rendering: native `plot_ly` vs ggplot→`ggplotly`

**Date:** 2026-07-01
**Status:** Discussion + decision. **Decision: leave both implementations as-is** (see below).
No code change required; this note records why, so the question is not re-litigated.

## Question

The ecosystem has **two** volcano implementations. Which should we standardize on — our
hand-written native `plotly::plot_ly` volcano, or the ggplot→`ggplotly` one used by the R6 plotter?
Is either dead code? When/why did the DEA report change?

## The two (three) implementations — verified from source

| Implementation | Engine | Entry point | Used by |
|---|---|---|---|
| `.multigroup_volcano` (`R/utilities.R`) → `ggplotly` | **ggplot** | `ContrastsPlotter$volcano_plotly()` (R6) | the prolfquapp **DEA report** (`plotter$volcano_plotly()$FDR`) |
| `.volcano` (`R/utilities.R:377`) | **native `plot_ly`** | `volcano_plotly()` (free fn, `R/utilities.R:469`) | prolfquasaint SAINT report, prophosqua overview, `prolfqua/tests/test-plotting_functions.R`, `prolfquapp/tests/test-CMD_DEA_CD.R` |
| `.scatter` (`R/utilities.R:516`) | **native `plot_ly`** | `scatter_plotly()` (free fn, `R/utilities.R:578`) | prophosqua `_Overview_PhosphoAndIntegration_site.Rmd` (site-vs-protein Δ scatter) |

**Nothing here is dead code.** Verified caller chains:

- Free `volcano_plotly()` → maps `.volcano` per contrast (`purrr::map2(..., .volcano, ...)`,
  `utilities.R:497-511`). It is exercised by a real prolfqua test
  (`test-plotting_functions.R:328-345`) and called by the SAINT + phospho report templates. So
  both `volcano_plotly` and `.volcano` are live.
- Free `scatter_plotly()` → maps `.scatter` per contrast (`utilities.R:608-620`). One live caller:
  the prophosqua phospho overview report. Live.
- `ContrastsPlotter$volcano_plotly()` (R6) → `.volcano` (R6 private) →
  `prolfqua:::.multigroup_volcano(...)` (a **ggplot**) → `plotly::ggplotly(p, tooltip="subject_id")`
  (`ContrastsPlotter.R:396,415`). Live (the DEA report).

## History — the change that prompted this

The **R6** `ContrastsPlotter$volcano_plotly` has *always* been ggplot→`ggplotly` (it never used
native `plot_ly`; the R6 volcano file has never contained `plot_ly`). The change was in the
**report**, not the method:

- **`6ebeb9c` (2026-05-14, "Expand SE Quarto report and add SAINT model integration")** switched the
  DEA report's volcano from `xd <- prolfqua::volcano_plotly(...)` (native plot_ly free function) to
  `plotter$volcano_plotly()$FDR` (R6 ggplotly). Confirmed by the diff (removed the
  `prolfqua::volcano_plotly(` line, added `plotter$volcano_plotly()$FDR`).

So before `6ebeb9c` the DEA report used the native-plotly volcano; after it, the R6 ggplotly one.
The SAINT and phospho reports still call the native `volcano_plotly` directly.

## Is our native `plot_ly` volcano "better"? — No, it's roughly a wash

Both are plotly-based and **both support crosstalk linked-brushing** (R6:
`plotly::highlight_key(~subject_id, group=self$group)`; native:
`crosstalk::SharedData$new(..., group=group)`). Crosstalk is therefore *not* a differentiator.

| Dimension | native `volcano_plotly` (plot_ly) | R6 `.multigroup_volcano` → ggplotly |
|---|---|---|
| crosstalk | ✅ SharedData | ✅ highlight_key |
| multi-contrast layout | list-of-plots → `subplot()` | `facet_wrap` (one plot) |
| ablines / annotations | precise (`layout(shapes)`, `add_annotations`) | ggplot `geom_hline/vline` + strip labels, **translated by ggplotly (lossy)** |
| fidelity | direct — no translation loss | ggplotly translation can mangle legends/geoms |
| return shape | list-per-contrast (awkward contract) | single faceted plotly (clean) |
| consistency w/ codebase | odd-one-out (everything else is ggplot) | matches the ggplot-first R6 plotter family |

Slight edge to **native** on fidelity/control (no ggplotly translation loss, exact ablines); slight
edge to **ggplot** on maintainability + consistency. Neither wins clearly.

## Why the R6 plotter uses ggplot (not our plotly)

The R6 decorators (`LFQDataPlotter`, `ContrastsPlotter`) are a **ggplot-based family** — boxplots,
heatmaps, density, PCA (`pca_plotly` is literally `ggplotly(self$pca())`), MA, scores are all
ggplot. The volcano follows that house style. The native `volcano_plotly`/`scatter_plotly` are
**standalone report helpers** at a different layer. Each is internally consistent with its layer.

## Decision

**Leave both as-is.** Rationale:

- Neither implementation is clearly better (wash).
- They serve two consistent layers: R6 plotter = ggplot family; standalone `*_plotly` = native
  plotly report helpers.
- Making `ContrastsPlotter$volcano_plotly` delegate to the native `volcano_plotly` is real plumbing
  (return-shape reconciliation so the report's `plotter$volcano_plotly()$FDR` keeps working) + risk,
  and would make the R6 volcano the *inconsistent* one among its ggplot siblings — for a wash-level
  benefit.
- The only genuine cost of the status quo is that the volcano **drawing logic** (FDR flooring via
  `.floor_significance_values`, ablines, crosstalk) is coded twice (`.multigroup_volcano` and
  `.volcano`) — a small, contained duplication. That fix has already had to be applied in both once.

## If revisited later

Only worth acting on if the duplicated drawing logic actually bites again. Options, least-to-most
disruptive:

1. **Do nothing** (current decision).
2. **Consolidate the DEA report** onto the native `volcano_plotly` (revert `6ebeb9c`): one report
   change, no core churn, keeps the tested native impl; but leaves the R6 ggplot volcano for
   programmatic callers.
3. **Single shared ggplot volcano builder** that both layers call, dropping `.volcano`. Cleanest
   de-duplication, but: changes the native functions' return contract (list-per-contrast → single
   faceted), touches prolfqua + prophosqua + prolfquasaint reports, and adopts ggplotly's lossy
   translation for ablines/annotations. Needs visual verification of all three reports.

Do **not** do a blanket "drop native → ggplotly" without a concrete driver — it trades a working,
tested, precise plotly impl for a lossier translated one, across production reports.
