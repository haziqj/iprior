gen_fbm <- function(n = 150, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  f <- function(x, truth = FALSE) {
    35 * dnorm(x, mean = 1, sd = 0.8) +
      65 * dnorm(x, mean = 4, sd = 1.5) +
      (x > 4.5) * (exp((1.25 * (x - 4.5))) - 1) +
      3 * dnorm(x, mean = 2.5, sd = 0.3)
  }
  x <- c(seq(0.2, 1.9, length = n * 5 / 8), seq(3.7, 4.6, length = n * 3 / 8))
  x <- sample(x, size = n)
  x <- x + rnorm(n, sd = 0.65)  # adding random fluctuation to the x
  x <- sort(x)
  y.err <- rt(n, df = 1)
  y <- f(x) + sign(y.err) * pmin(abs(y.err), rnorm(n, mean = 4.1))  # adding random
  data.frame(y = y, X = x)
}
