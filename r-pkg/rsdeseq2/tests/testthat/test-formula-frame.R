test_that("formula transport retains every frozen independent primitive value", {
  for(file in c("formula_expression_r462.json","formula_expression_extended_r462.json")) {
  reference <- jsonlite::fromJSON(.formula_reference_path(test_path("data"),file),simplifyVector=FALSE)
  variables <- lapply(reference$variables,.formula_decode)
  expect_length(reference$cases,if(file=="formula_expression_r462.json")129L else 59L)
  for(case in reference$cases) {
    request <- list(operation="evaluate",expressions=list(.formula_ast_decode(case$expression)),
      variables=variables,collation=reference$metadata$collation)
    if(is.null(case$expected))expect_error(rsdeseq2:::.rsdeseq2_formula_request(request),info=case$name)
    else {
      actual <- rsdeseq2:::.rsdeseq2_formula_request(request)[[1L]]
      expected <- .formula_descriptor(.formula_decode(case$expected))
      expect_identical(.formula_canonical(.formula_descriptor(actual)),.formula_canonical(expected),info=case$name)
    }
  }
  }
})

test_that("shared expression frames and offsets match live original R values and metadata", {
  d <- data.frame(x=seq(.5,6,length.out=12),y=seq(-3,4,length.out=12),
    g=factor(rep(c("C","a","B"),4),levels=c("unused","B","C","a")),
    ch=rep(c("a","A","ä","b"),3),l=rep(c(TRUE,FALSE,NA),4))
  d$M <- I(cbind(first=d$x,second=d$y));rownames(d$M)<-paste0("matrix",seq_len(nrow(d)))
  d$xi <- seq_len(nrow(d))
  d$z <- d$x;d$z[c(3,5)] <- c(NA_real_,NaN)
  d$date <- as.Date("2020-01-01")+seq_len(nrow(d))
  coefficients <- attr(stats::poly(d$x,3),"coefs")
  formulas <- list(~1,~g*x,~g/x,~0+g+x,~ch,~l,~date+x,x~g,
    ~I(cbind(x,y)),~scale(M),~poly(x,3),~poly(x,3,raw=TRUE),
    ~stats::poly(x,3),~poly(x,3,coefs=coefficients),~scale(x,center=FALSE),
    ~factor(ch),~relevel(droplevels(g),"C"),~ordered(g),
    ~I(x^3),~ifelse(is.na(z),0,z),~z+g,
    ~g+offset(log(x))+offset(I(x-y)),~offset(M),
    ~poly(M,degree=2),~poly(x,y,degree=2),~poly(M,degree=2,raw=TRUE),
    ~ifelse(g=="a",x,y),~I(xi-sum(xi)),~I(x-mean(x,.2)),
    ~ifelse(is.na(poly(z,2,raw=TRUE)),0,poly(z,2,raw=TRUE)))
  old <- options("na.action");on.exit(options(old),add=TRUE)
  for (action in c("na.omit","na.exclude","na.pass","na.fail")) {
    options(na.action=action)
    for (formula in formulas) {
      expected <- tryCatch(stats::model.frame(formula,d),error=identity)
      if(inherits(expected,"error")) {
        expect_error(rsdeseq2:::.rsdeseq2_model_frame(formula,d),"missing|model variable")
      } else {
        actual <- rsdeseq2:::.rsdeseq2_model_frame(formula,d)
        expect_true(identical(actual,expected,num.eq=FALSE,single.NA=FALSE),info=paste(action,deparse(formula)))
        expect_true(identical(rsdeseq2:::.rsdeseq2_model_offsets(actual),stats::model.offset(expected),num.eq=FALSE,single.NA=FALSE),info=deparse(formula))
      }
    }
  }
})

test_that("formula evaluation executes no R numerical transform or frame function", {
  d <- data.frame(x=seq(.5,6,length.out=12),g=factor(rep(c("A","B"),6)),ch=rep(c("a","b"),6))
  formula <- ~g+ch+I(cbind(log(x),x^2))+offset(log(x))
  expectedFrame <- stats::model.frame(formula,d)
  expected <- stats::model.matrix(formula,expectedFrame)
  expectedOffsets <- stats::model.offset(expectedFrame)
  functions <- list(base=c("log","scale","factor","cbind"),stats=c("poly","model.frame","model.frame.default","model.offset"))
  on.exit(for(pkg in names(functions))for(name in functions[[pkg]])suppressMessages(untrace(name,where=asNamespace(pkg))),add=TRUE)
  for(pkg in names(functions))for(name in functions[[pkg]])suppressWarnings(suppressMessages(trace(name,where=asNamespace(pkg),print=FALSE,tracer=quote(stop("R numerical formula evaluation")))))
  frame <- rsdeseq2:::.rsdeseq2_model_frame(formula,d)
  expect_identical(frame,expectedFrame)
  expect_identical(rsdeseq2:::.rsdeseq2_model_matrix(formula,frame),expected)
  expect_identical(rsdeseq2:::.rsdeseq2_model_offsets(frame),expectedOffsets)
})
