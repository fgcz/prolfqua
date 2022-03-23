library(plotly)
library(tidyr)
library(crosstalk)
m <- tidyr::gather(mtcars, variable, value, -vs)
#View(m)
msd <- highlight_key(m, ~variable)
#gg <- ggplot(msd, aes(factor(vs), value))  + #geom_jitter(alpha = 0.3)
#  geom_boxplot()


bscols(
  widths = c(11, 11),
  filter_select("id", "Select a variable", msd, ~variable, multiple = FALSE),
  #ggplotly(gg, dynamicTicks = "y") %>% layout(margin = list(l = 30))
  #ggplotly(gg, dynamicTicks = "y")# %>% layout(margin = list(l = 30))
  plot_ly(msd, x = ~vs, y = ~value, type="box") #%>% add_markers(alpha = 0.3)
)


# These examples demonstrate ways to display binned/aggregated selections
library(plotly)

d <- highlight_key(mpg)
dots <- plot_ly(d, colors = "Set1", color = ~class, x = ~displ, y = ~cyl) %>%
  plotly::layout(
    xaxis = list(title = "Engine displacement"),
    yaxis = list(title = "Number of cylinders")
  )
boxs <- plot_ly(d, colors = "Set1", color = ~class, x = ~class, y = ~cty) %>%
  add_boxplot() %>%
  plotly::layout(
    xaxis = list(title = ""),
    yaxis = list(title = "Miles per gallon (city)")
  )
bars <- plot_ly(d, colors = "Set1", x = ~class, color = ~class)

subplot(dots, boxs, titleX = TRUE, titleY = TRUE) %>%
  subplot(bars, nrows = 2, titleX = TRUE, titleY = TRUE) %>%
  plotly::layout(
    title = "Dynamic 2-way ANOVA (click & drag on scatterplot)",
    barmode = "overlay",
    showlegend = FALSE
  ) %>%
  highlight("plotly_selected", opacityDim = 0.6)


# These examples demonstrate ways to display binned/aggregated selections
library(plotly)

d <- highlight_key(mtcars)
sp <- plot_ly(d, x = ~mpg, y = ~disp) %>%
  add_markers(color = I("black"))

# 'statistical trace types'
hist <- plot_ly(d, x = ~factor(cyl)) %>%
  add_histogram(color = I("black"))
box <- plot_ly(d, y = ~disp, color = I("black")) %>%
  add_boxplot(name = " ")
violin <- plot_ly(d, y = ~disp, color = I("black")) %>%
  add_trace(type = "violin", name = " ")

subplot(sp, box, violin, shareY = TRUE, titleX = TRUE, titleY = TRUE) %>%
  subplot(hist, widths = c(.75, .25), titleX = TRUE, titleY = TRUE) %>%
  plotly::layout(
    barmode = "overlay",
    title = "Click and drag scatterplot",
    showlegend = FALSE
  ) %>%
  highlight("plotly_selected")


# These examples demonstrate ways to display binned/aggregated selections
library(plotly)

tx <- highlight_key(txhousing, ~city)
p1 <- ggplot(tx, aes(date, median, group = city)) + geom_line() + xlab(NULL)
gg1 <- ggplotly(p1, tooltip = c("city", "date", "median"))
p2 <- plot_ly(tx, x = ~median, color = I("black")) %>%
  add_histogram(histnorm = "probability density")

subplot(gg1, p2, titleX = TRUE, titleY = TRUE) %>%
  plotly::layout(barmode = "overlay") %>%
  highlight(dynamic = TRUE, selected = attrs_selected(opacity = 0.3))
