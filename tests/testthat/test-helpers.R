test_that("z_trans() adds z-scored columns", {
	D <- tibble::tibble(a = c(1, 2, 3), b = c(4, 5, 6))
	out <- z_trans(D, a)

	expect_equal(names(out), c("a", "b", "za"))
	expect_equal(out$za, as.numeric(scale(D$a)))
})

test_that("rescale helpers map values to the intended range", {
	expect_equal(rescale_unit(c(0, 5, 10)), c(0, .5, 1))
	expect_equal(rescale_zero_one(c(0, 5, 10)), c(0, .5, 1))
	expect_equal(rescale_centered(c(0, 10), scale = .5), c(2.5, 7.5))
})

test_that("reorder_levels() reorders a factor and validates positions", {
	x <- factor(c("a", "b", "c"))

	expect_equal(levels(reorder_levels(x, c(3, 1, 2))), c("c", "a", "b"))
	expect_error(reorder_levels(x, c(1, 2)))
	expect_error(reorder_levels(x, c(1, 2, 4)))
})

test_that("discard_redundant() drops constant columns", {
	D <- tibble::tibble(a = 1:3, b = c("x", "x", "x"), c = c(1, 0, 1))

	expect_equal(names(discard_redundant(D)), c("a", "c"))
	expect_equal(names(discard_redundant(D, except = "b")), c("a", "b", "c"))
})

test_that("discard_all_na() drops all-NA columns", {
	D <- tibble::tibble(a = c(1, NA), b = c(NA_real_, NA_real_), c = c(NA, 2))

	expect_equal(names(discard_all_na(D)), c("a", "c"))
})

test_that("update_by() mutates only the filtered rows", {
	D <- tibble::tibble(group = c(1, 1, 2, 2), value = c(4, 9, -4, -9))
	out <- update_by(D, group == 1, value = sqrt(value))

	expect_equal(out$value, c(2, 3, -4, -9))
	expect_equal(names(out), c("group", "value"))
})

test_that("left_union() appends the missing columns", {
	a <- tibble::tibble(x = 1:2)
	b <- tibble::tibble(x = 3:4, y = 5:6)

	expect_equal(names(left_union(a, b)), c("x", "y"))
	expect_equal(left_union(a, b)$y, 5:6)
	expect_error(left_union(a, tibble::tibble(x = 1)), "different number of rows")
})

test_that("go_first() and go_arrange() reorder columns", {
	D <- tibble::tibble(x = 1:3, y = 6:4, z = c(8, 9, 7))

	expect_equal(names(go_first(D, y)), c("y", "x", "z"))
	expect_equal(names(go_first(D, y:z)), c("y", "z", "x"))
	expect_equal(go_arrange(D, y)$y, c(4, 5, 6))
})

test_that("expand_grid_df() returns a tibble", {
	G <- expand_grid_df(a = 1:2, b = c("x", "y"))

	expect_s3_class(G, "tbl_df")
	expect_equal(names(G), c("a", "b"))
	expect_equal(nrow(G), 4)
})
