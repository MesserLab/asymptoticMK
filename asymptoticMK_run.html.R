#
#
#		asymptoticMK: Asymptotic McDonald-Kreitman Test web service
#
#		By Benjamin C. Haller and Philipp W. Messer
#		Copyright (C) 2017 Philipp Messer.
#
#		This web service should be available at http://benhaller.com/messerlab/asymptoticMK.html
#		The Github repository for asymptoticMK is at https://github.com/MesserLab/asymptoticMK
#
#

#  This file is part of asymptoticMK.
#
#  asymptoticMK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
#  asymptoticMK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with asymptoticMK.  If not, see <http://www.gnu.org/licenses/>.


#
#	Overrides from FastRWeb 1.1-2 update released by S. Urbanek on 2016-12-16 15:52:48
#	These would be unnecessary if I updated my FastRWeb install to the latest release
#

# From plots.R:

WebPlot <- function(width=640, height=480, type='png', ...) {
  file <- paste('tmp-',paste(sprintf('%x',as.integer(runif(4)*65536)),collapse=''),'.tmp',sep='')
  if (type %in% c('pdf', 'svg')) file <- paste(file, type, sep='.') ## some back-ends always append the extension, unfortuantely
  Cairo(width, height, type=type, file=file, ...)
  mime <- switch(type, png="image/png", pdf="application/pdf", jpg="image/jpeg", jpeg="image/jpeg", gif="image/gif", svg="image/svg+xml", "application/octet")
  structure(list(file=file,type=type,mime=mime,width=width,height=height), class="WebPlot")
}

as.character.WebPlot <- function(x, ...) {
  dev.off()
  tag <- "img"
  sz <- file.info(x$file)$size
  r <- readBin(x$file, raw(), sz)
  unlink(x$file)
  if (isTRUE(x$type == 'pdf')) tag <- "embed"
  paste0("<", tag, " src='",base64enc::dataURI(r, x$mime),"' width=",x$width," height=",x$height,">")
}

# From tools.R:

htmlEscape <- function(x) gsub(">","&gt;",gsub("<","&lt;",gsub("&","&amp;",x,fixed=TRUE),fixed=TRUE),fixed=TRUE)

otable <- function(..., tab='', tr='', cs='</td><td>', escape=TRUE) {
  a <- list(...)
  esc <- if (is.function(escape)) escape else if (isTRUE(escape)) htmlEscape else identity
  if (length(a)==1 && is.list(a[[1]])) a <- a[[1]]
  ml <- max(unlist(lapply(a,length)))
  m <- matrix(unlist(lapply(a,function(x) rep(as.character(x),length.out=ml))),ml)
  tout <- unlist(lapply(1:ml, function(x) paste("<tr",tr,"><td>",esc(paste(m[x,],collapse=cs)),"</td></tr>",sep='')))
  .e$out <- c(.e$out, paste("<table ",tab,">\n",sep=''), tout, '</table>')
}

ohead <- function(..., level=3, escape=TRUE) {
  esc <- if (is.function(escape)) escape else if (isTRUE(escape)) htmlEscape else identity 
  .e$out <- c(.e$out, paste("<h",level,">",esc(paste(...,sep='')),"</h",level,">",sep=''))
}

oprint <- function(..., sep='\n', escape=TRUE) {
  esc <- if (is.function(escape)) escape else if (isTRUE(escape)) htmlEscape else identity 
  .e$out <- c(.e$out, paste("<pre>",esc(paste(capture.output(print(...)),collapse=sep)),"</pre>",sep=''))
}

.opts <- function(..., disabled=FALSE, escape=TRUE) {
  esc <- if (is.function(escape)) escape else if (isTRUE(escape)) htmlEscape else identity
  l <- list(...)
  disabled <- if(isTRUE(disabled)) " disabled" else ""
  paste(
        if (length(l)) {
          n <- names(l)
          if (is.null(n) || any(n=="")) stop("Invalid unnamed argument")
          paste(unlist(lapply(seq.int(l), function(i) paste(" ",n[i],"=\"",gsub("\"","&quot;",esc(as.character(l[[i]])[1]),fixed=TRUE),"\"",sep=''))), collapse='')
        } else ""
        , disabled, sep='')
}

#
#	End code from FastRWeb 1.1-2  2016-12-16 15:52:48
#


# Get a CI using Monte Carlo simulation based upon a fitted model.  This is necessary because
# getting confidence intervals for non-linear models is a complicated business, apparently.
# Thanks to Andrej-Nikolai Spiess (http://www.dr-spiess.de) for this code.
# See: https://www.r-bloggers.com/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/
# Or, if that link goes stale: http://stats.stackexchange.com/a/251501/141766
predictNLS <- function(
object, 
newdata,
level = 0.95, 
nsim = 10000,
...
)
{
  require(MASS, quietly = TRUE)
   
  ## get right-hand side of formula
  RHS <- as.list(object$call$formula)[[3]]
  EXPR <- as.expression(RHS)
   
  ## all variables in model
  VARS <- all.vars(EXPR)
   
  ## coefficients
  COEF <- coef(object)
   
  ## extract predictor variable    
  predNAME <- setdiff(VARS, names(COEF))  
   
  ## take fitted values, if 'newdata' is missing
  if (missing(newdata)) {
    newdata <- eval(object$data)[predNAME]
    colnames(newdata) <- predNAME
  }
     
  ## check that 'newdata' has same name as predVAR
  if (names(newdata)[1] != predNAME) stop("newdata should have name '", predNAME, "'!")
   
  ## get parameter coefficients
  COEF <- coef(object)
     
  ## get variance-covariance matrix
  VCOV <- vcov(object)
   
  ## augment variance-covariance matrix for 'mvrnorm' 
  ## by adding a column/row for 'error in x'
  NCOL <- ncol(VCOV)
  ADD1 <- c(rep(0, NCOL))
  ADD1 <- matrix(ADD1, ncol = 1)
  colnames(ADD1) <- predNAME
  VCOV <- cbind(VCOV, ADD1)
  ADD2 <- c(rep(0, NCOL + 1))
  ADD2 <- matrix(ADD2, nrow = 1)
  rownames(ADD2) <- predNAME
  VCOV <- rbind(VCOV, ADD2) 
         
  ## iterate over all entries in 'newdata' as in usual 'predict.' functions
  NR <- nrow(newdata)
  respVEC <- numeric(NR)
  seVEC <- numeric(NR)
  varPLACE <- ncol(VCOV)   
   
  ## define counter function
  counter <- function (i) 
  {
    if (i%%10 == 0) 
      cat(i)
    else cat(".")
    if (i%%50 == 0) 
      cat("\n")
    flush.console()
  }
   
  outMAT <- NULL 
   
  for (i in 1:NR) {
    counter(i)
     
    ## get predictor values and optional errors
    predVAL <- newdata[i, 1]
    if (ncol(newdata) == 2) predERROR <- newdata[i, 2] else predERROR <- 0
    names(predVAL) <- predNAME  
    names(predERROR) <- predNAME  
     
    ## create mean vector for 'mvrnorm'
    MU <- c(COEF, predVAL)
     
    ## create variance-covariance matrix for 'mvrnorm'
    ## by putting error^2 in lower-right position of VCOV
    newVCOV <- VCOV
    newVCOV[varPLACE, varPLACE] <- predERROR^2
     
    ## create MC simulation matrix
    simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
     
    ## evaluate expression on rows of simMAT
    EVAL <- try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
    if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
     
    ## collect statistics
    PRED <- data.frame(predVAL)
    colnames(PRED) <- predNAME   
    FITTED <- predict(object, newdata = data.frame(PRED))
    MEAN.sim <- mean(EVAL, na.rm = TRUE)
    SD.sim <- sd(EVAL, na.rm = TRUE)
    MEDIAN.sim <- median(EVAL, na.rm = TRUE)
    MAD.sim <- mad(EVAL, na.rm = TRUE)
    QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2))
    RES <- c(FITTED, MEAN.sim, SD.sim, MEDIAN.sim, MAD.sim, QUANT[1], QUANT[2])
    outMAT <- rbind(outMAT, RES)
  }
   
  colnames(outMAT) <- c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
  rownames(outMAT) <- NULL
   
  cat("\n")
   
  return(outMAT)  
}

# core code for two-step nls2() model fit at a given level of precision (res)
fitMKmodel <- function(alpha_trimmed, f_trimmed, res)
{
	require(nls2)
	
	mod <- tryCatch({
			st <- expand.grid(const_a=seq(-1,1,length.out=res + 1), const_b=seq(-1,1,length.out=res), const_c=seq(1,10,length.out=res + 1))
			nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start=st, algorithm="brute-force", control=nls.control(maxiter=NROW(st)))
		},
		error=function(cond) {})
	
	if (length(mod) == 0)
		return(NULL)
	
	mod2 <- tryCatch({
			nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start = mod, control=nls.control(maxiter=200))
		},
		error=function(cond) {})
	
	if (length(mod2) == 0)
		return(NULL)
	
	return(mod2)
}

run <- function(...)
{
	require(nls2)
	require(MASS)
	
	#	The raw request is in the object "request"
	#
	if (request$uri != "/cgi-bin/R/asymptoticMK_run.html")
		return(as.WebResult("Request initiated from invalid URI."))
	if (request$method != "POST")
		return(as.WebResult("Request must be HTTP POST."))
	
	#	A parsed request is obtained from parse.multipart()
	#
	req <- parse.multipart()
	
	if (is.null(req) || is.na(req))
		return(as.WebResult("Request was not parseable with parse.multipart()."))
	
	#	Look for a "reply" field, used with curl to indicate that the user does not want a full HTML page response
	#
	reply_type <- "html"
	
	if (!is.null(req[["reply"]]))
	{
		reply <- req$reply
		
		if (!is.na(reply) && !is.null(reply))
		{
			if (reply == "table")
				reply_type <- "table"
			else
				return(as.WebResult("Unrecognized type for 'reply'."))
		}
		else
		{
			return(as.WebResult("Bad value for 'reply'."))
		}
	}
	
	#	Check the presence and type of all required response fields
	#
	if (is.null(req[["d0"]]))
		return(as.WebResult("Missing d0; please supply all parameters."))
	if (is.null(req[["d"]]))
		return(as.WebResult("Missing d; please supply all parameters."))
	if (is.null(req[["xlow"]]))
		return(as.WebResult("Missing xlow; please supply all parameters."))
	if (is.null(req[["xhigh"]]))
		return(as.WebResult("Missing xhigh; please supply all parameters."))
	
	d0 <- as.numeric(req$d0)
	d <- as.numeric(req$d)
	xlow <- as.numeric(req$xlow)
	xhigh <- as.numeric(req$xhigh)
	
	if (is.na(d0) || is.null(d0))
		return(as.WebResult("Malformed d0 (must be numeric)."))
	if (is.na(d) || is.null(d))
		return(as.WebResult("Malformed d (must be numeric)."))
	if (is.na(xlow) || is.null(xlow))
		return(as.WebResult("Malformed xlow (must be numeric)."))
	if (is.na(xhigh) || is.null(xhigh))
		return(as.WebResult("Malformed xhigh (must be numeric)."))
	
	if (is.null(req[["datafile"]]) || is.na(req[["datafile"]]))
		return(as.WebResult("Missing datafile; please supply this as a file upload."))
	
	datafile <- req$datafile
	
	if (!is.list(datafile))
		return(as.WebResult("Malformed datafile meta-information in response; not a list."))
	
	if (is.null(datafile[["filename"]]) || is.null(datafile[["tempfile"]]) || is.null(datafile[["content_type"]]))
		return(as.WebResult("Malformed datafile meta-information in response."))
	
	filename <- datafile$filename
	tempfile <- datafile$tempfile
	content_type <- datafile$content_type
	
	if (is.null(filename) || is.na(filename) || !is.character(filename))
		return(as.WebResult("Malformed datafile meta-information in response; filename is not a string."))
	if (is.null(tempfile) || is.na(tempfile) || !is.character(tempfile))
		return(as.WebResult("Malformed datafile meta-information in response; tempfile is not a string."))
	
	if (nchar(filename) < 1)
		return(as.WebResult("Malformed datafile meta-information in response; zero-length filename."))
	if (nchar(tempfile) < 1)
		return(as.WebResult("Malformed datafile meta-information in response; zero-length tempfile."))
	
	#	Bounds-check response variables
	#
	if (d0 <= 0)
		return(as.WebResult("d0 must greater than zero."))
	if (d <= 0)
		return(as.WebResult("d must greater than zero."))
	if ((xlow < 0.0) || (xlow > 1.0))
		return(as.WebResult("xlow must be in the interval [0,1]."))
	if ((xhigh < 0.0) || (xhigh > 1.0))
		return(as.WebResult("xhigh must be in the interval [0,1]."))
	if (xlow >= xhigh)
		return(as.WebResult("xlow must be less than xhigh."))
	
	if (!file_test("-f", tempfile))
		return(as.WebResult(paste("Uploaded file", tempfile, "does not exist.")))
	if (file.size(tempfile) == 0)
		return(as.WebResult("Uploaded file is zero-length."))
	
	#	Read in the file and check its format
	#
	df <- read.csv(tempfile, sep="	")
	
	if (NCOL(df) != 3)
		return(as.WebResult("Uploaded file does not contain exactly three tab-separated columns."))
	if (NROW(df) <= 0)
		return(as.WebResult("Uploaded file contains no data rows."))
	
	cols <- names(df)
	
	if (!is.na(as.numeric(cols[1])) || !is.na(as.numeric(cols[2])) || !is.na(as.numeric(cols[3])))
		return(as.WebResult("Uploaded file has a numeric column name; probably the required header row is missing."))
	
	f <- df[[1]]
	p <- df[[2]]
	p0 <- df[[3]]
	
	if ((length(f) != length(p)) || (length(f) != length(p0)))
		return(as.WebResult("Uploaded file has columns of unequal length."))
	if (!is.numeric(f))
		return(as.WebResult("The first column of the uploaded file, frequency, is not numeric."))
	if (!is.numeric(p))
		return(as.WebResult("The second column of the uploaded file, p, is not numeric."))
	if (!is.numeric(p0))
		return(as.WebResult("The third column of the uploaded file, p0, is not numeric."))
	if (any(is.na(f)))
		return(as.WebResult("The first column of the uploaded file, frequency, contains NA values (not allowed)."))
	if (any(is.na(p)))
		return(as.WebResult("The second column of the uploaded file, p, contains NA values (not allowed)."))
	if (any(is.na(p0)))
		return(as.WebResult("The third column of the uploaded file, p0, contains NA values (not allowed)."))
	if (any(f < 0.0) || any(f > 1.0))
		return(as.WebResult("The first column of the uploaded file, frequency, contains values out of the required range [0,1]."))
	if (any(p < 0))		# note that zero is allowed, although not recommended
		return(as.WebResult("The second column of the uploaded file, p, contains values < 0 (not allowed)."))
	if (all(p == 0))		# not all can be zero, however
		return(as.WebResult("The second column of the uploaded file, p, contains all values == 0 (not allowed)."))
	if (any(p0 <= 0))
		return(as.WebResult("The third column of the uploaded file, p0, contains values <= 0 (not allowed)."))
	
	if (NROW(df) < 3)
		return(as.WebResult("At least three data rows are required, to constrain the fit."))
	
	#	Compute alpha values and trim
	#
	alpha <- 1 - (d0/d) * (p/p0)
	cutoff_f1 <- xlow
	cutoff_f2 <- xhigh
	
	trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))
	
	if (sum(trim) < 3)
		return(as.WebResult("At least three data points are required after trimming the frequency range, to constrain the fit."))
	
	f_trimmed <- f[trim]
	alpha_trimmed <- alpha[trim]
	
	#	Compute the original McDonald-Kreitman alpha; we decided to use the trimmed data for this.
	#
	alpha_nonasymp <- 1 - (d0/d) * (sum(p[trim])/sum(p0[trim]))			# from trimmed data
	#alpha_nonasymp <- 1 - (d0/d) * (sum(p)/sum(p0))						# from untrimmed data
	
	#	Fit models
	#
	mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 10)
	
	if (length(mod1) == 0)
	{
		# try a deeper scan for a decent fit
		mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 20)
	}
	
	tryCatch({
			mod2 <- lm(alpha_trimmed ~ f_trimmed)
		},
		error=function(cond) {})
	
	linear_better <- FALSE
	
	if ((length(mod1) == 0) || (AIC(mod2) < AIC(mod1)))
		linear_better <- TRUE
	
	if (!linear_better)
	{
		# if we're leaning toward the exponential model, check for ridiculously wide confidence intervals; sometimes
		# we should reject the exponential model for that reason, because it is basically just a linear model with
		# a "cheat" of a swing up or down to fit one additional data point perfectly, which is lame :->
		ci_pred <- predictNLS(mod1, newdata=data.frame(f_trimmed=1.0))
		alpha_1_low <- ci_pred[6]
		alpha_1_high <- ci_pred[7]
		
		if ((alpha_1_low < -100) || (alpha_1_high > 100))
			linear_better <- TRUE
	}
	
	# Prepare for output and plotting
	full_seq <- seq(from=min(f), to=max(f), by=0.001)
	trimmed_seq <- seq(from=min(f_trimmed), to=max(f_trimmed), by=0.001)
	
	if (linear_better)
	{
		alpha_1_est <- predict(mod2, newdata=data.frame(f_trimmed=1.0))
		ci_pred <- predict(mod2, newdata=data.frame(f_trimmed=1.0), interval="confidence")	# we want confidence, not prediction
		alpha_1_low <- ci_pred[2]
		alpha_1_high <- ci_pred[3]
		const_a <- coef(mod2)["(Intercept)"]
		const_b <- coef(mod2)["f_trimmed"]
		const_c <- NA
		
		full_predicts <- predict(mod2, newdata=data.frame(f_trimmed=full_seq))
		trimmed_predicts <- predict(mod2, newdata=data.frame(f_trimmed=trimmed_seq))
		fit_color <- "red"
	}
	else
	{
		alpha_1_est <- predict(mod1, newdata=data.frame(f_trimmed=1.0))
		const_a <- coef(mod1)["const_a"]
		const_b <- coef(mod1)["const_b"]
		const_c <- coef(mod1)["const_c"]
		
		full_predicts <- predict(mod1, newdata=data.frame(f_trimmed=full_seq))
		trimmed_predicts <- predict(mod1, newdata=data.frame(f_trimmed=trimmed_seq))
		fit_color <- "red"
	}
	
	#	Output alternate reply data, if requested
	#
	if (reply_type == "table")
	{
		out("# Asymptotic McDonald&ndash;Kreitman Evaluator: Results")
		out("#")
		out("#    from: http://benhaller.com/messerlab/asymptoticMK.html")
		out("# ")
		out("# Analysis dataset:")
		out("# ")
		out(paste0("#    d0 = ", format(d0, scientific=FALSE)))
		out(paste0("#    d = ", format(d, scientific=FALSE)))
		out(paste0("#    input file = ", filename))
		out(paste0("#    x interval = [", format(xlow, digits=3, nsmall=3, scientific=FALSE), ", ", format(xhigh, digits=3, nsmall=3, scientific=FALSE), "]"))
		out("# ")
		
		if ((length(mod1) == 0) || linear_better)
			out("# Fitted model: linear, alpha(x) = a + bx")
		else
			out("# Fitted model: exponential, alpha(x) = a + b * exp(−cx)")
		
		out("")
		
		out(paste0("a\t", format(const_a, digits=5, nsmall=5, scientific=FALSE)));
		out(paste0("b\t", format(const_b, digits=5, nsmall=5, scientific=FALSE)));
		out(paste0("c\t", format(const_c, digits=5, nsmall=5, scientific=FALSE)));
		out(paste0("alpha_asymptotic\t", format(alpha_1_est, digits=5, nsmall=5, scientific=FALSE)));
		out(paste0("95% CI(lower)\t", format(alpha_1_low, digits=5, nsmall=5, scientific=FALSE)));
		out(paste0("95% CI(upper)\t", format(alpha_1_high, digits=5, nsmall=5, scientific=FALSE)));
		out(paste0("alpha_original\t", format(alpha_nonasymp, digits=5, nsmall=5, scientific=FALSE)));
		
		return(done())
	}
	
	#	BEGIN OUTPUT; from here on we should not use WebResult directly, since we are assembling a page
	#
	out("<!DOCTYPE html PUBLIC \"-//w3c//dtd html 4.0 transitional//en\">
<html>
<head>
<title>Asymptotic McDonald-Kreitman Test</title>
<style type=\"text/css\">
	body, td { 
		font:normal normal 100%/1.0 optima, times new roman, verdana, serif;
		line-height: 130%;
	}
	.math {
		font-family: times new roman, serif;
		white-space: nowrap;
	}
	sub { vertical-align:baseline; position:relative; top:0.3em; line-height:0; }
</style>
</head>
<body topmargin=50 leftmargin=50 rightmargin=50 bottommargin=50>
")
	
	out("
<h2 style=\"margin-bottom:5px;\">Asymptotic McDonald&ndash;Kreitman Test: Results</h2>

")
	
	#	Output analysis details
	#
	out(paste0("<p style=\"margin-top: 30px;\"><b>Analysis dataset:</b></p>

<blockquote><table cellspacing=0 cellpadding=0 border=0>
<tr><td height=25 valign=\"middle\"><span class=\"math\"><i>d</i><sub><small>0</small></sub></span></td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">", format(d0, scientific=FALSE), "</span></td></tr>
<tr><td height=25 valign=\"middle\"><span class=\"math\"><i>d</i></span></td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">", format(d, scientific=FALSE), "</span></td></tr>
<tr><td height=25 valign=\"middle\">Input file</td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><pre>", filename, "</pre></td></tr>
<tr><td height=25 valign=\"middle\"><span class=\"math\"><i>x</i></span> interval</td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">[", format(xlow, digits=3, nsmall=3, scientific=FALSE), ", ", format(xhigh, digits=3, nsmall=3, scientific=FALSE), "]</span></td></tr>
</table></blockquote>

<p style=\"margin-top: 30px;\"><b>Plots:</b></p>

"))
	
	#	PLOT 1: Frequency spectrum: p and p0 versus x
	#
	plot1 <- WebPlot(width=4 * 75, height=4 * 75, bg="white")						# PNG
	#plot1 <- WebPlot(width=4 * 75, height=4 * 75, bg="white", type="pdf", dpi=72)	# PDF
	par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
	plot(x=c(0,1), y=range(c(p,p0)), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab="polymorphism counts")
	points(x=f, y=p0, col="black", pch=19, cex=0.7)
	points(x=f, y=p, col="red", pch=19, cex=0.7)
	legend(x="topright", legend=c("p0 : neutral region", "p : test region"), col=c("black", "red"), pch=19, cex=0.9, pt.cex=0.7)
	
	out("<div>")
	out(plot1)		# Output of plots 1 and 2 are joined, see below
	
	#	PLOT 2: Frequency spectrum: p and p0 versus x
	#
	normalized_p0 <- p0 / sum(p0)
	normalized_p <- p / sum(p)
	
	plot2 <- WebPlot(width=4 * 75, height=4 * 75, bg="white")						# PNG
	#plot2 <- WebPlot(width=4 * 75, height=4 * 75, bg="white", type="pdf", dpi=72)	# PDF
	par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
	plot(x=c(0,1), y=range(c(normalized_p,normalized_p0)), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab="normalized SFS")
	points(x=f, y=normalized_p0, col="black", pch=19, cex=0.7)
	points(x=f, y=normalized_p, col="red", pch=19, cex=0.7)
	legend(x="topright", legend=c("p0 : neutral region", "p : test region"), col=c("black", "red"), pch=19, cex=0.9, pt.cex=0.7)
	
	#out("<div style=\"margin-top: 30px;\">")
	out("&nbsp;&nbsp;&nbsp;")
	out(plot2)
	out("</div>
<p>Fig. 1 [left].&nbsp;&nbsp;Polymorphism levels in the test region (<span class=\"math\"><i>p</i></span>, red points) and in the neutral reference region (<span class=\"math\"><i>p</i><sub><small>0</small></sub></span>, black points), as a function of derived allele frequency <span class=\"math\"><i>x</i></span>.</p>")
	out("
<p>Fig. 2 [right].&nbsp;&nbsp;Normalized site frequency spectra (SFS) in the test region (red points) and the neutral reference region (black points).  This plot shows the data from Fig. 1, normalized such that <span class=\"math\">&Sigma;<sub><small><i>x</i></small></sub>&nbsp;<i>p</i><sub><small>0</small></sub>(<i>x</i>)&nbsp;=&nbsp;&Sigma;<sub><small><i>x</i></small></sub>&nbsp;<i>p</i>(<i>x</i>)&nbsp;=&nbsp;1</span> for purposes of comparison.</p>")
	
	#	PLOT 3: alpha(x) ~ x
	#
	plot3 <- WebPlot(width=4 * 75, height=4 * 75, bg="white")						# PNG
	#plot3 <- WebPlot(width=4 * 75, height=4 * 75, bg="white", type="pdf", dpi=72)	# PDF
	par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
	plot(x=c(0,1), y=range(alpha), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab=expression(paste("MK ", alpha, "(", x, ")")))
	points(x=f, y=alpha, col="black", pch=19, cex=0.7)
	
	out("<div style=\"margin-top: 30px;\">")
	out(plot3)		# Output of plots 3 and 4 are joined, see below
	
	#	PLOT 4: alpha(x) ~ x plus fit to that data
	#
	yr <- na.omit(c(alpha, alpha_1_est, max(0.0, alpha_1_low), min(1.0, alpha_1_high), alpha_nonasymp))
	
	plot4 <- WebPlot(width=4 * 75, height=4 * 75, bg="white")						# PNG
	#plot4 <- WebPlot(width=4 * 75, height=4 * 75, bg="white", type="pdf", dpi=72)	# PDF
	par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
	plot(x=c(0,1), y=range(yr), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab=expression(paste("MK ", alpha, "(", x, ")")))
	polygon(x=c(-1, -1, 2.0, 2.0), y=c(alpha_1_low, alpha_1_high, alpha_1_high, alpha_1_low), col="#DDDDDD", border=NA)
	
	if (!is.na(alpha_nonasymp))
		abline(h=alpha_nonasymp, lty=3, col="#777777")
	
	points(x=f, y=alpha, col=ifelse(trim, "black", "#999999"), pch=19, cex=0.7)
	abline(v=cutoff_f1, col="#5555FF")
	abline(v=cutoff_f2, col="#5555FF")
	
	lines(x=full_seq, full_predicts, col="#333333", lwd=1)
	lines(x=trimmed_seq, trimmed_predicts, col=fit_color, lwd=2)
	abline(h=alpha_1_est, lty=2, col=fit_color)
	
	box()
	
	#out("<div style=\"margin-top: 30px;\">")
	out("&nbsp;&nbsp;&nbsp;")
	out(plot4)
	out("</div>
<p>Fig. 3 [left].&nbsp;&nbsp;McDonald&ndash;Kreitman <span class=\"math\"><i>&alpha;</i>(<i>x</i>) = 1 &minus; (<i>d</i><sub><small>0</small></sub> / <i>d</i>) (<i>p</i>(<i>x</i>) / <i>p</i><sub><small>0</small></sub>(<i>x</i>))</span> versus <span class=\"math\"><i>x</i></span>.</p>")
	out("
<p>Fig. 4 [right].&nbsp;&nbsp;The asymptotic McDonald&ndash;Kreitman test results.  This plot shows the data from Fig. 3, with fitting information superimposed.  The blue vertical lines indicate the cutoff interval for the polymorphism data; points outside of the cutoff interval are plotted in gray, indicating that they were not used in the fit.  ")
	if (linear_better)
		out("The red curve shows the best fit to the data within the cutoff interval for a function <span class=\"math\"><i>&alpha;</i><sub><small>fit</small></sub>(<i>x</i>) = <i>a</i> + <i>bx</i></span>.  The dashed red horizontal line ")
	else
		out("The red curve shows the best fit to the data within the cutoff interval for a function <span class=\"math\"><i>&alpha;</i><sub><small>fit</small></sub>(<i>x</i>) = <i>a</i> + <i>b</i> exp(&minus;<i>cx</i>)</span>.  The dashed red horizontal line ")
	out("shows the estimate of <span class=\"math\"><i>&alpha;</i><sub><small>asymptotic</small></sub></span> from the fitted function, and the gray band indicates the 95% confidence interval around that estimate.  Finally, the dotted gray horizontal line shows <span class=\"math\"><i>&alpha;</i><sub><small>original</small></sub></span>, the estimate from the original non-asymptotic McDonald&ndash;Kreitman test (also using only the data within the cutoff interval), for comparison; note that use of this value is not recommended.</p>")
	
	#	Output fit details
	#
	out("<p style=\"margin-top: 30px;\"><b>Fitted <span class=\"math\"><i>&alpha;</i>(<i>x</i>):</b></p>

");
	
	if (length(mod1) == 0)
		out("<p>The exponential fit failed to converge (usually because the data are not exponential in shape); the linear model is therefore reported here.</p>")
	else if (linear_better)
		out("<p>The linear model <span class=\"math\"><i>&alpha;</i>(<i>x</i>) = <i>a</i> + <i>bx</i></span> was better (by AIC) than the exponential model, and is therefore reported here.</p>")
	else
		out("<p>The exponential model <span class=\"math\"><i>&alpha;</i>(<i>x</i>) = <i>a</i> + <i>b</i> exp(&minus;<i>cx</i>)</span> was better (by AIC) than the linear model, and is therefore reported here.</p>")
	
	out(paste0("<blockquote><table cellspacing=0 cellpadding=0 border=0>
<tr><td height=25 valign=\"middle\"><span class=\"math\"><i>a</i></span></td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">", format(const_a, digits=5, nsmall=5, scientific=FALSE), "</span></td></tr>
<tr><td height=25 valign=\"middle\"><span class=\"math\"><i>b</i></span></td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">", format(const_b, digits=5, nsmall=5, scientific=FALSE), "</span></td></tr>
<tr><td height=25 valign=\"middle\"><span class=\"math\"><i>c</i></span></td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">", format(const_c, digits=5, nsmall=5, scientific=FALSE), "</span></td></tr>
</table></blockquote>
"))
	
	out("

<p style=\"margin-top: 30px;\"><b>Estimates of <span class=\"math\"><i>&alpha;</i>:</b></p>

");
	out("<p>The result of the asymptotic McDonald&ndash;Kreitman test is given by <span class=\"math\"><i>&alpha;</i><sub><small>asymptotic</small></sub></span>; this value is obtained by extrapolating the above fitted function to <span class=\"math\"><i>x</i>&nbsp;=&nbsp;1</span>.  The 95% confidence interval around the estimated value of <span class=\"math\"><i>&alpha;</i><sub><small>asymptotic</small></sub></span> is also shown.  The value of <span class=\"math\"><i>&alpha;</i><sub><small>original</small></sub></span>, from the original non-asymptotic McDonald&ndash;Kreitman test, is also given here for comparison but its use is not recommended.  Both <span class=\"math\"><i>&alpha;</i></span> estimates are derived from the polymorphism frequency data within the supplied cutoff interval for <span class=\"math\"><i>x</i></span>.</p>")
	
	out(paste0("<blockquote><table cellspacing=0 cellpadding=0 border=0>
<tr><td height=25 valign=\"middle\"><span class=\"math\"><i>&alpha;</i><sub><small>asymptotic</small></sub></span></td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">", format(alpha_1_est, digits=5, nsmall=5, scientific=FALSE), "</span></td><td>&nbsp;&nbsp;&nbsp;&nbsp;</td></tr>
<tr><td height=25 valign=\"middle\"><span class=\"math\">95% CI(lower)</span></td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">", format(alpha_1_low, digits=5, nsmall=5, scientific=FALSE), "</span></td></tr>
<tr><td height=25 valign=\"middle\"><span class=\"math\">95% CI(upper)</span></td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">", format(alpha_1_high, digits=5, nsmall=5, scientific=FALSE), "</span></td></tr>
<tr style=\"color: #666666\"><td height=45 valign=\"middle\"><span class=\"math\"><i>&alpha;</i><sub><small>original</small></sub></span></td><td>&nbsp;&nbsp;=&nbsp;&nbsp;</td><td><span class=\"math\">", format(alpha_nonasymp, digits=5, nsmall=5, scientific=FALSE), "</span></td><td>&nbsp;&nbsp;&nbsp;&nbsp;</td></tr>
</table></blockquote>
"))
	
	#	Output final page content
	#
	
	out("<p style=\"margin-top: 30px;\"><b>Citation:</b></p>

<p>If you use this service, please cite our paper:</p>
<blockquote><i>[not yet published, please check back for a citation...]</i></blockquote>

<hr width=\"75%\" style=\"margin-top: 30px; margin-bottom: 30px;\">

<p>Please let us know of any issues with this service at <i>philipp {dot} messer &lt;at&gt; gmail [dot] com</i>.  Thanks!</p>

</body>
</html>
")
	done()
}
















