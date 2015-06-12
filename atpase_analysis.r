require(RColorBrewer)

colfun = colorRampPalette(c('purple','cornflowerblue','seagreen', 'yellow','orange','tomato2'), space='Lab')

getSlopev2 <- function(data) { 
	t <- data[,1]
	plot(data[,1], data[,2], pch=16, col='#3399cc')
	print('Pick region to fit')
	fitregion <- locator(2)$x
  print(fitregion)
	trange <- c(findClosest(t, fitregion[1]), findClosest(t, fitregion[2]))
	ty_begin <- which(t==trange[1])
	ty_end <- which(t==trange[2])
	fit_y <- lm(data[,2][ty_begin:ty_end] ~ t[ty_begin:ty_end])
	lines(t[ty_begin:ty_end], predict(fit_y), col='red', lty=2)
	sprintf('Slope : %1.4f µM/sec', (coef(fit_y)[2]/6.22 * 1000))
	
	}

getSlopeAll = function(data) {
	t = data[,1]
	w = data[,2:length(data)]
	datarange = c(min(w), max(w))
	cols = colfun(length(data)-1)
	
	plot(range(t), datarange, type='n', las=1, xlab='Time', ylab='NADH (uM)')
	for (i in 2:length(data))
	{
		points(t, data[,i], pch=16, col=cols[i-1])
	}	
	
	print('Pick region to fit')
	fitregion = locator(2)$x
	print(fitregion)
	trange = c(findClosest(t, fitregion[1]), findClosest(t, fitregion[2]))
	ty_begin = which(t==trange[1])
	ty_end   = which(t==trange[2])
	wt = t[ty_begin:ty_end]
	fits = c()
	for (i in 2:length(data))
	{
		wy = data[,i][ty_begin:ty_end]
		wfit = lm(wy ~ wt)
		lines(wt, predict(wfit), col='black', lwd=2)
		fits = append(fits, (coef(wfit)[[2]]))	
	}
	return(fits) 
}	

findClosest <- function(numlist, num) {
	
	distnum <- abs(numlist - num)
	whichmin <- which(distnum == min(distnum))
	return(numlist[whichmin])
	
}

getSlopeRaw <- function(data) { 
	t <- data[,1]
	plot(data[,1], data[,2], pch=16, col='#3399cc')
	#print('Pick region to fit')
	fitregion <- locator(2)$x
	trange <- c(findClosest(t, fitregion[1]), findClosest(t, fitregion[2]))
	ty_begin <- which(t==trange[1])
	ty_end <- which(t==trange[2])
	fit_y <- lm(data[,2][ty_begin:ty_end] ~ t[ty_begin:ty_end])
	lines(t[ty_begin:ty_end], predict(fit_y), col='red', lty=2)
	#sprintf('%1.4f µM/sec', coef(fit_y)[2])
	return(-coef(fit_y)[[2]])
	}
