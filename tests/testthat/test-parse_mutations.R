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
