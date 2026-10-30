test_that("custom calls fail before formula evaluation, including masked and nested calls", {
  calls <- 0L
  custom <- function(x) { calls <<- calls+1L; x }
  data <- data.frame(x=seq_len(8L))
  formulas <- list(~custom(x), ~log(custom(x)), ~custom(x)-custom(x),
    ~I((function(x) custom(x))(x)), ~base::do.call(custom, list(x)),
    ~utils::head(x), ~custom(x):x)
  for (formula in formulas) {
    expect_error(rsdeseq2:::.rsdeseq2_model_frame(formula, data),
      "custom or unsupported formula function")
  }
  log <- custom
  expect_error(rsdeseq2:::.rsdeseq2_model_frame(~log(x), data),
    "custom or unsupported formula function 'log'")
  expect_identical(calls, 0L)
  for (kind in c("active", "delayed", "alias")) {
    scope <- new.env(parent=environment())
    if (kind=="active") makeActiveBinding("sqrt",function(value) {
      calls <<- calls+1L; base::sqrt
    },scope)
    if (kind=="delayed") delayedAssign("sqrt",{calls <<- calls+1L; base::sqrt},
      eval.env=environment(),assign.env=scope)
    if (kind=="alias") scope$sqrt <- base::sqrt
    formula <- stats::as.formula("~sqrt(x)",env=scope)
    expect_error(rsdeseq2:::.rsdeseq2_model_frame(formula,data),
      "custom or unsupported formula function 'sqrt'",info=kind)
  }
  expect_identical(calls,0L)
  # Explicit qualification resolves the established built-in even when masked.
  expect_identical(rsdeseq2:::.rsdeseq2_model_frame(~base::log(x),data),
    stats::model.frame(~base::log(x),data))
})

test_that("built-in formula metadata is unchanged by the call guard", {
  d <- data.frame(x=seq(.5,4,length.out=8), group=factor(rep(c("A","B"),4)))
  for (formula in list(~group*x,~group/x,~(group+x)^2,~0+group+offset(log(x)),
      ~I(x^2)+scale(x),~stats::poly(x,2),~relevel(group,"B"),~.,x~group)) {
    expect_identical(rsdeseq2:::.rsdeseq2_model_frame(formula,d), stats::model.frame(formula,d))
  }
})

test_that("public full and reduced designs apply the same custom-call policy", {
  skip_if_not_installed("DESeq2")
  d <- data.frame(condition=factor(rep(c("A","B"),each=3)),x=seq_len(6))
  counts <- matrix(rep(c(12L,19L,31L,26L,29L,42L),3),nrow=3,byrow=TRUE)
  object <- DESeq2::DESeqDataSetFromMatrix(counts,d,~condition)
  calls <- 0L
  custom <- function(x) { calls <<- calls+1L; x }
  changed <- object
  methods::slot(changed,"design") <- ~custom(x)
  expect_error(DESeq(changed),"custom or unsupported formula function 'custom'")
  expect_error(DESeq(object,test="LRT",reduced=~custom(x)),
    "custom or unsupported formula function 'custom'")
  expect_identical(calls,0L)
})
