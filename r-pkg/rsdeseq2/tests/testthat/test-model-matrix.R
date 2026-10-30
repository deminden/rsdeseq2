test_that("shared model-matrix assembly preserves independently frozen values and order", {
  r<-jsonlite::fromJSON(test_path("data","model_matrix_r462.json"),simplifyVector=FALSE)
  expect_length(r$cases,48L)
  for(case in r$cases) {
    input<-.model_reference_input(case)
    direct<-suppressWarnings(do.call(.Call,c(list("rsdeseq2_model_matrix"),input,list(PACKAGE="rsdeseq2"))))
    wrapped<-suppressWarnings(rsdeseq2:::.rsdeseq2_model_matrix(~1,.model_reference_frame(case)))
    expect_identical(as.character(colnames(direct)),as.character(unlist(case$expected_names)),info=case$name)
    expect_identical(attr(direct,"assign"),as.integer(unlist(case$assignments)),info=case$name)
    expect_identical(.row_reference_bits(as.vector(t(direct))),.row_reference_bits(.row_reference_values(case$expected)),info=case$name)
    expect_identical(unname(as.vector(wrapped)),unname(as.vector(direct)),info=case$name)
    expect_identical(colnames(wrapped),colnames(direct),info=case$name)
  }
  expect_length(r$average_contrasts,10L)
  for(case in r$average_contrasts) {
    run<-function().Call("rsdeseq2_average_contrast",as.integer(case$coefficients),
      as.integer(unlist(case$positive))+1L,as.integer(unlist(case$negative))+1L,PACKAGE="rsdeseq2")
    if(is.null(case$expected))expect_error(run(),"average contrast",info=case$name)
    else expect_identical(.row_reference_bits(run()),.row_reference_bits(.row_reference_values(case$expected)),info=case$name)
  }
})

test_that("formula adapters retain live matrix attributes and use no R matrix or contrast calculation", {
  d<-data.frame(a=factor(rep(c("A","B","C"),4)),b=factor(rep(c("X","Y"),6)),x=seq(-2,3,length.out=12))
  d[["cell type"]]<-d$a
  d$date<-as.Date("2020-01-01")+seq_len(nrow(d));d$y<-12:1;d$single<-factor(rep("only",nrow(d)))
  scenarios<-list(~a*b,~0+a*x,~x+a,~a/b,~factor(a),~I(cbind(x,x+1)),~`cell type`*x,~date+x,y~a,single~x)
  for(formula in scenarios) {
    frame<-stats::model.frame(formula,d)
    expected<-stats::model.matrix(formula,frame)
    functions<-c("model.matrix","model.matrix.default","contr.treatment","contr.sum","contr.helmert","contr.poly","contr.SAS")
    for(name in functions)suppressMessages(trace(name,where=asNamespace("stats"),print=FALSE,tracer=quote(stop("R numerical matrix calculation"))))
    actual<-tryCatch(rsdeseq2:::.rsdeseq2_model_matrix(formula,frame),finally={for(name in functions)suppressMessages(untrace(name,where=asNamespace("stats")))})
    expect_identical(actual,expected)
  }
  frame<-stats::model.frame(y~y*x,d)
  expected<-suppressWarnings(stats::model.matrix(y~y*x,frame))
  messages<-character()
  actual<-withCallingHandlers(rsdeseq2:::.rsdeseq2_model_matrix(y~y*x,frame),warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart("muffleWarning")})
  expect_identical(actual,expected)
  expect_length(messages,2L)
  expect_match(messages[[1]],"response appeared on the right-hand side")
})
