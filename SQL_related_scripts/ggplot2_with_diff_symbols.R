library(ggplot2)

# Your data
data <- data.frame(
  x = rep(1:5, 5),
  y = rep(1:5, each = 5),
  symbol = sample(c("A", "B", "C", "D"), 25, replace = TRUE)
)

# Create a named vector for symbol mapping
symbol_map <- c(
  A = "▲",
  B = "■",
  C = "●",
  D = "◆"
)

# Plot
ggplot(data, aes(x, y)) +
  geom_tile(aes(fill = symbol), color = "white") +
  geom_text(aes(label = symbol_map[symbol]), size = 6) +
  scale_fill_discrete(name = "Symbol") +
  coord_equal() +
  theme_minimal()