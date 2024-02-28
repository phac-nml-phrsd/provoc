test_that("parse_mutations from various output of Art's function", {
    expect_equal(parse_mutation("~", "2832", "G"), "aa:orf1a:K856R")
    expect_equal(parse_mutation("~", 2832, "G"), "aa:orf1a:K856R")


    expect_equal(parse_mutation("~", "22942", "C"), "T22942C")
    expect_equal(parse_mutation("~", "5797", "G"), "T5797G")

    expect_equal(parse_mutation("~", "27785", "A"), "aa:orf7b:Y10*")
    expect_equal(parse_mutation("~", "27895", "G"), "aa:orf8:M1R")
    expect_equal(parse_mutation("~", "5755", "G"), "aa:orf1a:N1830K")

    expect_equal(parse_mutation("-", "27782", "4"), "del:27782:4")
    expect_equal(parse_mutation("-", "12339", "1"), "del:12339:1")
    expect_equal(parse_mutation("-", "22772", "1"), "del:22772:1")

    expect_equal(parse_mutation("+", "27966", "G"), "ins:27966:1")
    expect_equal(parse_mutation("+", "27804", "C"), "ins:27804:1")
    expect_equal(parse_mutation("+", "27953", "TAGTTTACAGTCATGTACTCAACACCAACCATATT"), "ins:27953:35")
})

test_that("parse_mutations returns the correct vector", {
  expect_equal(parse_mutations(c("~2832G", "~22942C", "~5797G")), c("aa:orf1a:K856R", "T22942C", "T5797G"))
  expect_equal(parse_mutations(c("~27785A","-27782.4", "+27966:G")), c("aa:orf7b:Y10*", "del:27782:4", "ins:27966:1"))
  expect_equal(parse_mutations(c("~27895G", "-12339.1", "+27804:C", "+27953:ABCD")), c("aa:orf8:M1R", "del:12339:1", "ins:27804:1", "ins:27953:4"))
})
