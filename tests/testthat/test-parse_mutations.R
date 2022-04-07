test_that("parse_mutations from various output of Art's function", {
    expect_equal(parse_mutation("~", "2832", "G"), "aa:orf1a:K856R")
    expect_equal(parse_mutation("~", 2832, "G"), "aa:orf1a:K856R")
    expect_equal(parse_mutation("-", "27782", "4"), "del:27783:4")
    expect_equal(parse_mutation("~", "22942", "C"), "T22942C")
    expect_equal(parse_mutation("+", "27966", "G"), "ins:27966:1")
    expect_equal(parse_mutation("~", "27785", "A"), "aa:orf7b:Y10*")
})
