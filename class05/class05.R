#' ---
#' title: 'Class 05: Data Visualization'
#' author: "Yuhan Zhang"
#' output: github_document
#' always_allow_html: true
#' ---

# Class 05: Data Visualization
# install.packages("ggplot2")
library(ggplot2)

head(cars)

# All ggplot have at lieast 3 layers,
# data + aes + geoms
p <- ggplot(data = cars) + 
  aes(x = speed, y = dist) +
  geom_point() +
  labs(title = "Stopping Distance of Old Cars", 
       x = "Speed (MPH)", 
       y = "Stopping Distance (ft)") + 
  geom_smooth(method = "lm", formula = 'y ~ x',  aes(x = speed, y = dist))

show(p)

# Side-note: ggplot is not the only graphics systems
# a very popular one is good old "base" R graphics
plot(cars)

# Plot some gene expression results.
# Dataset is online in tab separated format so we
# use the read.delim() function to import into R
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)

# Q. How many genes in this dataset
nrow(genes)

# Q. How many genes are "up"?
table(genes$State)

# Q. What % are up?
round(table(genes$State)/nrow(genes) * 100, 2)

p <- ggplot(genes, aes(x = Condition1, y = Condition2, col = State)) + 
  geom_point()

show(p)
p <- p + scale_color_manual(values = c("blue", "grey", "red")) + 
  labs(title = "Gene Expression Changes Upon Drug Treatment", 
       x = "Control (no drug", 
       y = "Drug Treatment")

show(p)

# Let's explore the gapminder dataset
# install.packages("gapminder")
library(gapminder)
head(gapminder)

# Let's make a new plot of year vs lifeExp
p <- ggplot(gapminder, aes(x = year, y = lifeExp, col = continent)) +
      geom_jitter(width = 0.3, alpha = 0.4) +
      geom_violin(aes(group = year), alpha = 0.2, 
                  draw_quantiles = 0.5) 

show(p)

# Install the plotly
# install.packages("plotly")
library(plotly)
ggplotly()




