# DWNN: Distance-Weighted Nearest Neighbors

dwnn <- function(query, dataset, sigma) {
	if (!is.matrix(dataset)) {
		return ("dataset must be a data frame.")
	}
	classId = ncol(dataset)
	# Conjunto de entrada
	X = matrix(dataset[,1:(classId-1)], ncol=classId-1)
	# Conjunto de saída
	Y = dataset[,classId]
	E = apply(X, 1, function(row) { sqrt(sum((row - query)^2)) } )
	act = exp(-E^2 / (2*sigma^2))

	return ( act%*%Y / sum(act) )
}

test.identity <- function(sigma) {

	# função identidade
	x = cbind(-5:5,-5:5)
	plot(x)

	points = NULL
	for (query in seq(-5,5,length=100)) {
		y = dwnn(query, x, sigma)
		points = rbind(points, cbind(query, y))
	}
	points(points, col=2)

}

test.identity2 <- function(sigma) {

	# função identidade
	x = cbind(-5:5,-5:5)
	plot(x, xlim=c(-10,10))

	points = NULL
	for (query in seq(-10,10,length=100)) {
		y = dwnn(query, x, sigma)
		points = rbind(points, cbind(query, y))
	}
	points(points, col=2)

}

sin.test <- function(sigma=0.00025) {

	time = seq(0,9,length=1000)
	series = sin(2*pi*time)

	# AMI -> tseriesChaos::mutual -> d=4
	# False-Nearest Neighbors -> tseriesChaos::false.nearest -> m=2

	## Dimensões para o Teorema de Takens
	m = 2
	d = 4

	dataset = tseriesChaos::embedd(series, m=m, d=d)
	train.set = 0.7
	last = floor(nrow(dataset)*train.set)

	train = dataset[1:last,]
	test = dataset[(last+1):nrow(dataset),]
	buffer = train
	# x1 x5
	# x2 x6
	# x3 x7
	# x4 x8
	# x5 x9
	# x6 x10
	# x7 x11
	# x8 x12
	#=======
	# x9 x13'
	# x10 x14'

	for (i in (last+1):nrow(dataset)) {
		query = buffer[i - d, 2:ncol(buffer)]
		y = dwnn(query, train, sigma)
		buffer = rbind(buffer, cbind(query, y))
	}

	plot(dataset[,ncol(dataset)])
	points(buffer[,ncol(dataset)], col=2) # aproximação

	ret = list()
	ret$dataset = dataset
	ret$train = train
	ret$buffer = buffer

	return (ret)
}

logistic.test <- function(sigma) {

	series = c()
	x = 0.5
	for (i in 1:1000) {
		series = c(series, x)
		x = 3.8 * x * (1 - x)
	}

	## Dimensões para o Teorema de Takens
	m = 2
	d = 1

	dataset = tseriesChaos::embedd(series, m=m, d=d)
	train.set = 0.7
	last = floor(nrow(dataset)*train.set)

	train = dataset[1:last,]
	test = dataset[(last+1):nrow(dataset),]
	buffer = train
	# x1 x5
	# x2 x6
	# x3 x7
	# x4 x8
	# x5 x9
	# x6 x10
	# x7 x11
	# x8 x12
	#=======
	# x9 x13'
	# x10 x14'

	for (i in (last+1):nrow(dataset)) {
		query = buffer[i - d, 2:ncol(buffer)]
		y = dwnn(query, train, sigma)
		buffer = rbind(buffer, cbind(query, y))
	}

	plot(dataset[(last+1):(last+1+50),ncol(dataset)], t="l")
	lines(buffer[(last+1):(last+1+50),ncol(dataset)], col=2) # aproximação

	ret = list()
	ret$dataset = dataset
	ret$train = train
	ret$buffer = buffer

	return (ret)
}

mypredict <- function(series, m, d, sigma, train.set=0.7, plot=F, decomposition=F, pca=F) {

	if (decomposition) {
		cat("Finding the deterministic component...\n")
		res = EMD::emd(series, boundary="wave")
		deterministic = rowSums(res$imf[,4:res$nimf]) + res$residue

		cat("Reconstructing the space phase....\n")
		dataset = tseriesChaos::embedd(deterministic , m=m, d=d)
	
	} else {
		cat("Reconstructing the space phase....\n")
		dataset = tseriesChaos::embedd(series , m=m, d=d)
	}
	
	if (pca) {
		dataset.pca = prcomp(dataset[,1:ncol(dataset)-1], center = TRUE, scale = TRUE)
		print((summary(dataset.pca)))
	}

	#plot(dataset)
	last = floor(nrow(dataset)*train.set)

	cat("Separating train and test sets....\n")
	train = dataset[1:last,]
	test = dataset[(last+1):nrow(dataset),]
	buffer = train

	cat("Predicting....\n")
	for (i in (last+1):nrow(dataset)) {
		query = as.numeric(buffer[i - d, 2:ncol(buffer)])
		y = as.numeric(dwnn(query, train, sigma))

		buffer = rbind(buffer, c(query, y))
	}

	cat("Plot....\n")

	if (plot) {
		plot(dataset[(last+1):(last+1+50),ncol(dataset)], t="l")
		lines(buffer[(last+1):(last+1+50),ncol(dataset)], col=2) # aproximação
	}

	cat("Error...\n")

	error = dtw::dtw(dataset[(last+1):(last+1+100),ncol(dataset)],buffer[(last+1):(last+1+100),ncol(dataset)])$normalizedDistance
	cat(error,"\n")

	cat("End\n")
	ret = list()
	ret$last = last
	ret$dataset = dataset
	ret$train = train
	ret$buffer = buffer

	return (ret)
}

mypredict.tune <- function(series, m, d, train.set,
			 sigma.start, sigma.end, len) {

	M = NULL
	for (sigma in seq(sigma.start, sigma.end, len=len)) {
		#cat("sigma = ", sigma, "\n")
		res = mypredict(series, m, d, sigma, train.set)
		classId = ncol(res$dataset)
		expected = 
		  res$dataset[(res$last+1):nrow(res$dataset),classId]
		obtained = 
		  res$buffer[(res$last+1):nrow(res$dataset),classId]
	  	#error = sqrt(sum((expected - obtained)^2)) / (nrow(res$dataset) - res$last)
		error = dtw::dtw(expected,obtained)$normalizedDistance
		M = rbind(M, cbind(sigma, error))
	}

	return (M)
}

lag.tune <- function(lag.initial,lag.end,series,sigma,m,train.set=0.7){
	
	M = NULL

	for (d in lag.initial:lag.end){
		res = mypredict(series, m, d, sigma, train.set,plot=T)
		classId = ncol(res$dataset)
		expected = 
		  res$dataset[(res$last+1):nrow(res$dataset),classId]
		obtained = 
		  res$buffer[(res$last+1):nrow(res$dataset),classId]
		error = dtw::dtw(expected,obtained)$normalizedDistance
		M = rbind(M, cbind(d, error))
	}

	return (M)
}














