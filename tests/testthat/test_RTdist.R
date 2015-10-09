context("Checking Likelihood Function Sanity")

preds <- CFR_PCR(summarised=FALSE)
dist <- preds %>% group_by(class,obsOrder) %>% RTdist()

test_that("Total Likelihood sums to the number of groups/classes", {
  expect_equal(sum(dist$y), length(unique(dist$class)))
  }
)
