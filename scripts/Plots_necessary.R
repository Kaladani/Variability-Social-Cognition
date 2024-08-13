library(ggplot2)
library(viridis)
library(ggpubr)

theme_set(theme_bw(base_size = 10,
                   base_family = "sans"))



rnormt <- function(n, range, mu, s) {
  
  # range is a vector of two values
  
  F.a <- pnorm(min(range), mean = mu, sd = s)
  F.b <- pnorm(max(range), mean = mu, sd = s)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qnorm(u, mean = mu, sd = s)
  
}

## Necessary condition


n <- 200

pre <- rep(0, n)

pre

X <- runif(n, 0, 1)

# Y <- pre + 1*X * rbinom(n, 1, 0.5)
# Change <- Y - pre

Y <- pre + runif(n, 0, X)

noise <- sapply(Y, function(y) rnormt(1, c(-y,1-y), 0, 0.05))
Y <- Y + noise

plot(X, Y)

Change <- Y - pre
plot(X, Change)

mydata_a <- data.frame(pre = pre, X = X, post = Y, delta = Change)
p_a <- mydata_a |> 
  ggplot(aes(x = X, y = delta)) + 
  geom_point(size=2, alpha=0.5) +
  # geom_smooth(method = "lm", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.5) +
  labs(x = "X", y = expression(Delta(Y))) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_a



## Sufficient condition

n <- 200

pre <- 0

pre

X <- runif(n, 0, 1)

Y <- runif(n, X, 1)
noise <- sapply(Y, function(y) rnormt(1, c(-y,1-y), 0, 0.05))
Y <- Y + noise

plot(X,Y)
  # rbinom(n, 1, 0.5)
Change <- Y - pre

plot(X, Change)

mydata_b <- data.frame(pre = pre, X = X, post = Y, delta = Change)
p_b <- mydata_b |> 
  ggplot(aes(x = X, y = delta)) + 
  geom_point(size=2, alpha=0.5) +
  # geom_smooth(method = "lm", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.5) +
  labs(x = "X", y = expression(Delta(Y))) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_b


## necessary & sufficient

n <- 200

pre <- 0

pre

X <- runif(n, 0, 1)

Y <- X
noise <- sapply(Y, function(y) rnormt(1, c(-y,1-y), 0, 0.05))
Y <- Y + noise

plot(X,Y)
# rbinom(n, 1, 0.5)
Change <- Y - pre

plot(X, Change)

mydata_c <- data.frame(pre = pre, X = X, post = Y, delta = Change)
p_c <- ggplot(data = mydata_c, aes(x = X, y = delta)) + 
  geom_point(size=2, alpha=0.5) +
  # geom_smooth(method = "lm", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.5) +
  labs(x = "X", y = expression(Delta(Y))) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_c


## neither

n <- 200

pre <- 0

pre

X <- runif(n, 0, 1)
Z <- runif(n, 0, 1)

Y <- Z
noise <- sapply(Y, function(y) rnormt(1, c(-y,1-y), 0, 0.05))
Y <- Y + noise

plot(X,Y)
# rbinom(n, 1, 0.5)
Change <- Y - pre

plot(X, Change)

mydata_d <- data.frame(pre = pre, X = X, post = Y, delta = Change)
p_d <- mydata_d |>
  ggplot(aes(x = X, y = delta)) +
  geom_point(size=2, alpha=0.5) +
  # geom_smooth(method = "lm", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.5) +
  labs(x = "X", y = expression(Delta(Y))) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_d


p_necessary <- ggarrange(p_a, p_c, p_b, p_d, ncol = 2, nrow=2, labels = c("(a)", "(c)", "(b)", "(d)"), 
                         vjust = 1.5, 
                         font.label = list(size = 8, color = "black", face = "bold", family = "sans"))

p_necessary

ggsave(file="plots/Plot_Necessary-Sufficient-Data.tiff", p_necessary, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/Plot_Necessary-Sufficient-Data.jpeg", p_necessary, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)
