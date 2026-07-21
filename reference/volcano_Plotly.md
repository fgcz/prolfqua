# volcano plotly

volcano plotly

## Usage

``` r
volcano_plotly(
  .data,
  effect = "fc",
  significance = "BFDR",
  contrast = "condition",
  proteinID = "Prey",
  color = "modelName",
  palette = NULL,
  xintercept = c(-2, 2),
  yintercept = 0.1,
  minsignificance = NULL,
  title_size = 25,
  group = "BB"
)
```

## Arguments

- .data:

  data frame

- effect:

  effect size x-axis

- significance:

  significance

- contrast:

  column with contrast labels

- proteinID:

  column with protein ids

- color:

  column used for colouring points

- palette:

  named colour vector for the colour aesthetic

- xintercept:

  vertical abline at x

- yintercept:

  horizontal abline at y

- minsignificance:

  optional minimum significance value (floor for -log10 axis). If NULL,
  only exact zero and negative values are replaced with the smallest
  positive observed value in the same score column.

- title_size:

  font size of the subplot title annotation

- group:

  crosstalk group name for linked brushing

## Value

The requested plot, table, or transformed object.

## Examples

``` r
data <- data.frame(fc = c(-1,0,1,2,8), BFDR = c(0.01,1, 0.01, 0.005,0),
condition = rep("A",5), Prey = LETTERS[1:5], modelName = c("A","A","B","A","A"))

dataB <- data.frame(fc = c(-1,0,1,2,8), BFDR = c(0.01,1, 0.01, 0.005,0),
condition = rep("B",5), Prey = LETTERS[1:5],modelName = c("A","A","B","B","B"))
data <- dplyr::bind_rows(data, dataB)
bc <- volcano_plotly(data, xintercept = 1, yintercept= 0.01, palette = c(A = "black" , B = "red"))
bc |> plotly::subplot()

{"x":{"data":[{"x":[-1,0,2,8],"y":[2,-0,2.3010299956639813,2.3010299956639813],"mode":"markers","text":["A","B","D","E"],"showlegend":false,"type":"scatter","key":["A","B","D","E"],"set":"BB","name":"A","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y","_isNestedKey":false,"frame":null},{"x":[1],"y":[2],"mode":"markers","text":"C","showlegend":false,"type":"scatter","key":["C"],"set":"BB","name":"B","marker":{"color":"rgba(255,0,0,1)","line":{"color":"rgba(255,0,0,1)"}},"textfont":{"color":"rgba(255,0,0,1)"},"error_y":{"color":"rgba(255,0,0,1)"},"error_x":{"color":"rgba(255,0,0,1)"},"line":{"color":"rgba(255,0,0,1)"},"xaxis":"x","yaxis":"y","_isSimpleKey":true,"_isNestedKey":false,"frame":null},{"x":[-1,0],"y":[2,-0],"mode":"markers","text":["A","B"],"showlegend":false,"type":"scatter","key":["A","B"],"set":"BB","name":"A","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x2","yaxis":"y2","_isNestedKey":false,"frame":null},{"x":[1,2,8],"y":[2,2.3010299956639813,2.3010299956639813],"mode":"markers","text":["C","D","E"],"showlegend":false,"type":"scatter","key":["C","D","E"],"set":"BB","name":"B","marker":{"color":"rgba(255,0,0,1)","line":{"color":"rgba(255,0,0,1)"}},"textfont":{"color":"rgba(255,0,0,1)"},"error_y":{"color":"rgba(255,0,0,1)"},"error_x":{"color":"rgba(255,0,0,1)"},"line":{"color":"rgba(255,0,0,1)"},"xaxis":"x2","yaxis":"y2","_isNestedKey":false,"frame":null}],"layout":{"xaxis":{"domain":[0,0.47999999999999998],"automargin":true,"anchor":"y"},"xaxis2":{"domain":[0.52000000000000002,1],"automargin":true,"anchor":"y2"},"yaxis2":{"domain":[0,1],"automargin":true,"anchor":"x2"},"yaxis":{"domain":[0,1],"automargin":true,"anchor":"x"},"annotations":[{"text":"A","x":-1,"y":2.3010299956639813,"showarrow":false,"xanchor":"left","font":{"size":25},"xref":"x","yref":"y"},{"text":"B","x":-1,"y":2.3010299956639813,"showarrow":false,"xanchor":"left","font":{"size":25},"xref":"x2","yref":"y2"}],"shapes":[{"type":"line","y0":0,"y1":1,"yref":"paper","x0":1,"x1":1,"line":{"color":"green","dash":"dot"},"xref":"x"},{"type":"line","x0":0,"x1":0.47999999999999998,"xref":"paper","y0":2,"y1":2,"line":{"color":"red","dash":"dot"},"yref":"y"},{"type":"line","y0":0,"y1":1,"yref":"paper","x0":1,"x1":1,"line":{"color":"green","dash":"dot"},"xref":"x2"},{"type":"line","x0":0.52000000000000002,"x1":1,"xref":"paper","y0":2,"y1":2,"line":{"color":"red","dash":"dot"},"yref":"y2"}],"images":[],"margin":{"b":40,"l":60,"t":25,"r":10},"dragmode":"zoom","hovermode":"closest","showlegend":false},"attrs":{"2c83775087f5":{"x":{},"y":{},"mode":"markers","text":{},"showlegend":false,"color":{},"colors":["black","red"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"},"2c831af8bc40":{"x":{},"y":{},"mode":"markers","text":{},"showlegend":false,"color":{},"colors":["black","red"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0,"ctGroups":["BB"]},"subplot":true,"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
