# scatter plotly

scatter plotly

## Usage

``` r
scatter_plotly(
  .data,
  dx = "diff.protein",
  dy = "diff.site",
  contrast = "condition",
  proteinID = "protein_Id",
  color = "modelName",
  palette = NULL,
  title_size = 25,
  group = "BB"
)
```

## Arguments

- .data:

  data frame

- dx:

  column name for the x-axis difference

- dy:

  column name for the y-axis difference

- contrast:

  column with contrast labels

- proteinID:

  column with protein ids

- color:

  column used for colouring points

- palette:

  named colour vector for the colour aesthetic

- title_size:

  font size of the subplot title annotation

- group:

  crosstalk group name for linked brushing

## Value

The requested plot, table, or transformed object.

## Examples

``` r
data <- data.frame(diff.protein = c(-1,0,1,2,8), diff.site = c(0.01,1, 0.01, 0.005,0),
condition = rep("A",5), protein_Id = LETTERS[1:5], modelName = c("A","A","B","A","A"))

dataB <- data.frame(diff.protein = c(-1,0,1,2,8), diff.site = c(0.01,1, 0.01, 0.005,0),
condition = rep("B",5), protein_Id = LETTERS[1:5],modelName = c("A","A","B","B","B"))
data <- dplyr::bind_rows(data, dataB)
bc <- scatter_plotly(data, palette = c(A = "black" , B = "red"))
bc[[1]]

{"x":{"visdat":{"1b922ab3dfeb":["function () ","plotlyVisDat"]},"cur_data":"1b922ab3dfeb","attrs":{"1b922ab3dfeb":{"x":{},"y":{},"mode":"markers","text":{},"showlegend":false,"color":{},"colors":["black","red"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"annotations":[{"text":"A","x":-1,"y":0,"showarrow":false,"xanchor":"left","font":{"size":25}}],"xaxis":{"domain":[0,1],"automargin":true,"title":"diff.protein"},"yaxis":{"domain":[0,1],"automargin":true,"title":"diff.site"},"dragmode":"zoom","hovermode":"closest","showlegend":false},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-1,0,2,8],"y":[0.01,1,0.0050000000000000001,0],"mode":"markers","text":["A","B","D","E"],"showlegend":false,"type":"scatter","key":["A","B","D","E"],"set":"BB","name":"A","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y","_isNestedKey":false,"frame":null},{"x":[1],"y":[0.01],"mode":"markers","text":"C","showlegend":false,"type":"scatter","key":["C"],"set":"BB","name":"B","marker":{"color":"rgba(255,0,0,1)","line":{"color":"rgba(255,0,0,1)"}},"textfont":{"color":"rgba(255,0,0,1)"},"error_y":{"color":"rgba(255,0,0,1)"},"error_x":{"color":"rgba(255,0,0,1)"},"line":{"color":"rgba(255,0,0,1)"},"xaxis":"x","yaxis":"y","_isSimpleKey":true,"_isNestedKey":false,"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0,"ctGroups":["BB"]},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}bc[[2]]

{"x":{"visdat":{"1b925353e913":["function () ","plotlyVisDat"]},"cur_data":"1b925353e913","attrs":{"1b925353e913":{"x":{},"y":{},"mode":"markers","text":{},"showlegend":false,"color":{},"colors":["black","red"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"annotations":[{"text":"B","x":-1,"y":0,"showarrow":false,"xanchor":"left","font":{"size":25}}],"xaxis":{"domain":[0,1],"automargin":true,"title":"diff.protein"},"yaxis":{"domain":[0,1],"automargin":true,"title":"diff.site"},"dragmode":"zoom","hovermode":"closest","showlegend":false},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-1,0],"y":[0.01,1],"mode":"markers","text":["A","B"],"showlegend":false,"type":"scatter","key":["A","B"],"set":"BB","name":"A","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y","_isNestedKey":false,"frame":null},{"x":[1,2,8],"y":[0.01,0.0050000000000000001,0],"mode":"markers","text":["C","D","E"],"showlegend":false,"type":"scatter","key":["C","D","E"],"set":"BB","name":"B","marker":{"color":"rgba(255,0,0,1)","line":{"color":"rgba(255,0,0,1)"}},"textfont":{"color":"rgba(255,0,0,1)"},"error_y":{"color":"rgba(255,0,0,1)"},"error_x":{"color":"rgba(255,0,0,1)"},"line":{"color":"rgba(255,0,0,1)"},"xaxis":"x","yaxis":"y","_isNestedKey":false,"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0,"ctGroups":["BB"]},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}bc |> plotly::subplot()

{"x":{"data":[{"x":[-1,0,2,8],"y":[0.01,1,0.0050000000000000001,0],"mode":"markers","text":["A","B","D","E"],"showlegend":false,"type":"scatter","key":["A","B","D","E"],"set":"BB","name":"A","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y","_isNestedKey":false,"frame":null},{"x":[1],"y":[0.01],"mode":"markers","text":"C","showlegend":false,"type":"scatter","key":["C"],"set":"BB","name":"B","marker":{"color":"rgba(255,0,0,1)","line":{"color":"rgba(255,0,0,1)"}},"textfont":{"color":"rgba(255,0,0,1)"},"error_y":{"color":"rgba(255,0,0,1)"},"error_x":{"color":"rgba(255,0,0,1)"},"line":{"color":"rgba(255,0,0,1)"},"xaxis":"x","yaxis":"y","_isSimpleKey":true,"_isNestedKey":false,"frame":null},{"x":[-1,0],"y":[0.01,1],"mode":"markers","text":["A","B"],"showlegend":false,"type":"scatter","key":["A","B"],"set":"BB","name":"A","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x2","yaxis":"y2","_isNestedKey":false,"frame":null},{"x":[1,2,8],"y":[0.01,0.0050000000000000001,0],"mode":"markers","text":["C","D","E"],"showlegend":false,"type":"scatter","key":["C","D","E"],"set":"BB","name":"B","marker":{"color":"rgba(255,0,0,1)","line":{"color":"rgba(255,0,0,1)"}},"textfont":{"color":"rgba(255,0,0,1)"},"error_y":{"color":"rgba(255,0,0,1)"},"error_x":{"color":"rgba(255,0,0,1)"},"line":{"color":"rgba(255,0,0,1)"},"xaxis":"x2","yaxis":"y2","_isNestedKey":false,"frame":null}],"layout":{"xaxis":{"domain":[0,0.47999999999999998],"automargin":true,"anchor":"y"},"xaxis2":{"domain":[0.52000000000000002,1],"automargin":true,"anchor":"y2"},"yaxis2":{"domain":[0,1],"automargin":true,"anchor":"x2"},"yaxis":{"domain":[0,1],"automargin":true,"anchor":"x"},"annotations":[{"text":"A","x":-1,"y":0,"showarrow":false,"xanchor":"left","font":{"size":25},"xref":"x","yref":"y"},{"text":"B","x":-1,"y":0,"showarrow":false,"xanchor":"left","font":{"size":25},"xref":"x2","yref":"y2"}],"shapes":[],"images":[],"margin":{"b":40,"l":60,"t":25,"r":10},"dragmode":"zoom","hovermode":"closest","showlegend":false},"attrs":{"1b922ab3dfeb":{"x":{},"y":{},"mode":"markers","text":{},"showlegend":false,"color":{},"colors":["black","red"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"},"1b925353e913":{"x":{},"y":{},"mode":"markers","text":{},"showlegend":false,"color":{},"colors":["black","red"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0,"ctGroups":["BB"]},"subplot":true,"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
