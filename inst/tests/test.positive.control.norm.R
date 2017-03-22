library(testthat);
suppressMessages(require(NanoStringNormCNV));

# We want to hide all futile logger calls while testing
# Testthat already has nice output messages describing errors
flog.threshold(FATAL);


data.raw <- read.markup.RCC(
	rcc.path = paste(getwd(),"/data",sep=''),
	rcc.pattern = "*.RCC"
	);
nano.raw <- data.raw$x;
test_value=positive.control.qc(nano.raw);


# Pull expected values
cs1x <- read.table(file="data/postive.control.norm.expect.txt",sep=' ',header=TRUE,stringsAsFactors=FALSE);
context("Positve Control Norm");



test_that("Normalized Data Frame is correct", {
	expect_equivalent(test_value, cs1x);
	});


