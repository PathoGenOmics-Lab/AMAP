library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]

data <- read.table(input_file, header = FALSE, sep = "\t")

p <- ggplot(data, aes(x = V4, y = V2, fill = V4)) +
  geom_boxplot() +
  geom_point(aes(color = V4), position = position_jitter(width = 0.2)) +
  facet_wrap(~ V1, scales = "free_y") +
  theme_bw() +
  labs(x = "V4", y = "V2")

print(p)

ggsave("cases.png", plot = p, width = 10, height = 5)