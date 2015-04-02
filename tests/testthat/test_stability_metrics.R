library(OmicsMarkeR)
context("Stability Metrics")

test_that("pairwise.stability requires appropriate inputs", {
    some.numbers <- seq(20)
    
    # matrix of Metabolites identified (e.g. 5 trials)
    features <- 
        replicate(5, paste("Metabolite", sample(some.numbers, 10), sep="_"))
    
    expect_error(pairwise.stability(features, "kuncheva"), 
                 'argument "nc" is missing, with no default')
})

test_that("pairwise.model.stability requires appropriate inputs", {
    some.numbers <- seq(20)
    
    plsda <- 
        replicate(5, paste("Metabolite", sample(some.numbers, 10), sep="_"))
    rf <-
        replicate(5, paste("Metabolite", sample(some.numbers, 10), sep="_"))
    
    features <- list(plsda=plsda, rf=rf)
    
    expect_error(pairwise.model.stability(features, "kuncheva"), 
                 'argument "nc" is missing, with no default')
})

test_that("RPT requires appropriate inputs", {
    expect_error(RPT("A", 0.80), "stability is not of type 'numeric'")
    expect_error(RPT(0.80, "B"), "performance is not of type 'numeric'")
    expect_error(RPT(0.80,0.80,"C"), "beta is not of type 'numeric'")
    expect_error(RPT(1.1,.80), "stability is not in range 0 to 1.")
    expect_error(RPT(.80,1.1), "performance is not in range 0 to 1.")
    expect_error(RPT(.80,.80,-1), "beta is non-positive.")
})