library(flexmix)
library(circlus)

mix <- rbind(rPKBD_ACG(300, 0.95, c(1,0,0)), rPKBD_ACG(300, 0.9, c(-1,0,0)))
m1 <- flexmix(mix ~ 1, k = 2, model = circlus::PKBDNN_clust())
m1

#df = readRDS("~/Documents/GitHub/PKBD---code/ExamplesAndTesting/df_final.RDS")
#load("~/Documents/GitHub/PKBD---code/data/Abstracts.rda")
data("Abstracts")
#df = data.frame(df)
head(Abstracts)
names(Abstracts)[c(279,280,281,282)]


GTE <-  t(sapply(Abstracts[,279], function(x) x))
GTE <- GTE/sqrt(rowSums(GTE^2))
OAI <-  t(sapply(Abstracts[,280], function(x) x))
OAI512 <-  t(sapply(Abstracts[,281], function(x) x))
OAI256 <-  t(sapply(Abstracts[,282], function(x) x))
pages = matrix(Abstracts$pages, ncol = 1)

PKBD_abstract_6 <- flexmix(OAI256 ~ 1, k = 6, model = circlus::PKBD_clust())
PKBD_abstract_6
PKBD_abstract_6@logLik
sapply(PKBD_abstract_6@components, function(x) x[[1]]@parameters$rho)

PKBDNN_abstract_6 <- flexmix(OAI ~ 1, k = 6, model = circlus::PKBDNN_clust(LR= 0.01, free_iter = 5),
                             control = list(verbose = 1))
PKBDNN_abstract_6
PKBDNN_abstract_6@logLik

PKBDNN_abstract_8 <- flexmix(OAI ~ 1, k = 8, model = circlus::PKBDNN_clust(LR = 0.05),
                              control = list(verbose = 1))

PKBDNN_abstract_6b <- flexmix(OAI ~ 1 + pages, k = 6, model = circlus::PKBDNN_clust(LR = 0.01),
                              control = list(verbose = 1, nrep = 10))
#saveRDS(PKBDNN_abstract_6b, file = "ExamplesAndTesting/PKBD6to5withpagesand107kloglik.RDS")



PKBDNN_abstract_6b_adam <- flexmix(OAI ~ 1 + pages, k = 6, 
                                   model = circlus::PKBDNN_clust_adam(LR = 0.1, EPOCHS = 100, free_iter =  5),
                                   control = list(verbose = 1))

reduced_OAI <- readRDS("~/Documents/GitHub/PKBD---code/ExamplesAndTesting/df_reduced.RDS")
PKBDNN_abstract_6 <- flexmix(GTE ~ 1, k = 6, model = circlus::PKBDNN_clust(LR= 0.1, free_iter = 5),
                             control = list(verbose = 1))
PKBD_abstract_6 <- flexmix(GTE ~ 1, k = 6, model = circlus::PKBD_clust())
PKBD_abstract_6
PKBD_abstract_6@logLik

PKBDNN_abstract_6b_adam <- flexmix(GTE ~ 1 + pages, k = 6, 
                                   model = circlus::PKBDNN_clust(LR = 0.01, EPOCHS = 100, free_iter =  3),
                                   control = list(verbose = 1))

PKBDNN_abstract_6b <- flexmix(GTE ~ 1 + pages, k = 6, model = circlus::PKBDNN_clust(LR = 0.01, free_iter = 3),
                              control = list(verbose = 1, nrep = 10))


#####################################################################

set.seed(1)
SC_abstract_8 <- flexmix(OAI256 ~ 1, k = 8, model = circlus::SCauchy_clust())
SC_abstract_8
SC_abstract_8@logLik
sapply(SC_abstract_8@components, function(x) x[[1]]@parameters$rho)


set.seed(1)
torch::torch_manual_seed(1)
SCNN_abstract_8 <- flexmix(OAI256 ~ 1, k = 8, model = circlus::SCauchyNN_clust(LR= 0.02, free_iter = 5),
                             control = list(verbose = 1))
#1 Genetic Influences on Psychiatric and Behavioral Disorders
#2 Advancements in Model Validation and Benchmarking Techniques
#3 Travel Behavior, Environmental Impact, and Social Factors in Transportation and Tourism
#4 Environmental and Biological Effects of Agrochemicals
#5 Advancing Market Segmentation Techniques and Applications
#6 Advancing Biopharmaceutical Production Through Data-Driven Approaches
#7 Advancing Finite Mixture Models and Their Applications
#8 Advancing Clustering Techniques in Data Analysis


library(tm)
library(wordcloud)
library(reshape)
library(tm)

comparison.cloud <- function (term.matrix, scale = c(4, 0.5), max.words = 300, random.order = FALSE, 
                              rot.per = 0.1, colors = brewer.pal(max(3, ncol(term.matrix)), 
                                                                 "Dark2"), use.r.layout = FALSE, title.size = 3, title.colors = NULL, 
                              match.colors = FALSE, title.bg.colors = "grey90", ...) {
  ndoc <- ncol(term.matrix)
  thetaBins <- seq(from = 0, to = 2 * pi, length = ndoc + 1)
  for (i in 1:ndoc) {
    term.matrix[, i] <- term.matrix[, i]/sum(term.matrix[, 
                                                         i])
  }
  mean.rates <- rowMeans(term.matrix)
  for (i in 1:ndoc) {
    term.matrix[, i] <- term.matrix[, i] - mean.rates
  }
  group <- apply(term.matrix, 1, function(x) which.max(x))
  words <- rownames(term.matrix)
  freq <- apply(term.matrix, 1, function(x) max(x))
  tails <- "g|j|p|q|y"
  last <- 1
  nc <- length(colors)
  overlap <- function(x1, y1, sw1, sh1) {
    if (!use.r.layout) 
      return(wordcloud:::is_overlap(x1, y1, sw1, sh1, boxes))
    s <- 0
    if (length(boxes) == 0) 
      return(FALSE)
    for (i in c(last, 1:length(boxes))) {
      bnds <- boxes[[i]]
      x2 <- bnds[1]
      y2 <- bnds[2]
      sw2 <- bnds[3]
      sh2 <- bnds[4]
      if (x1 < x2) 
        overlap <- x1 + sw1 > x2 - s
      else overlap <- x2 + sw2 > x1 - s
      if (y1 < y2) 
        overlap <- overlap && (y1 + sh1 > y2 - s)
      else overlap <- overlap && (y2 + sh2 > y1 - s)
      if (overlap) {
        last <<- i
        return(TRUE)
      }
    }
    FALSE
  }
  ord <- rank(-freq, ties.method = "random")
  words <- words[ord <= max.words]
  freq <- freq[ord <= max.words]
  group <- group[ord <= max.words]
  if (random.order) {
    ord <- sample.int(length(words))
  } else {
    ord <- order(freq, decreasing = TRUE)
  }
  words <- words[ord]
  freq <- freq[ord]
  group <- group[ord]
  thetaStep <- 0.05
  rStep <- 0.05
  plot.new()
  op <- par("mar")
  par(mar = c(0, 0, 0, 0))
  plot.window(c(0, 1), c(0, 1), asp = 1)
  normedFreq <- freq/max(freq)
  size <- (scale[1] - scale[2]) * normedFreq + scale[2]
  boxes <- list()
  docnames <- colnames(term.matrix)
  if (!is.null(title.colors)) {
    title.colors <- rep(title.colors, length.out = ndoc)
  }
  title.bg.colors <- rep(title.bg.colors, length.out = ndoc)
  adj <- -0.03
  rr <- c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE)
  for (i in 1:length(words)) {
    rotWord <- runif(1) < rot.per
    r <- 0
    theta <- runif(1, 0, 2 * pi)
    x1 <- 0.5
    y1 <- 0.5
    wid <- strwidth(words[i], cex = size[i])
    ht <- strheight(words[i], cex = size[i])
    if (grepl(tails, words[i])) 
      ht <- ht + ht * 0.2
    if (rotWord) {
      tmp <- ht
      ht <- wid
      wid <- tmp
    }
    isOverlaped <- TRUE
    while (isOverlaped) {
      inCorrectRegion <- theta > thetaBins[group[i]] && 
        theta < thetaBins[group[i] + 1]
      if (inCorrectRegion && !overlap(x1 - 0.5 * wid, y1 - 
                                      0.5 * ht, wid, ht) && x1 - 0.5 * wid > 0 && y1 - 
          0.5 * ht > 0 && x1 + 0.5 * wid < 1 && y1 + 0.5 * 
          ht < 1) {
        text(x1, y1, words[i], cex = size[i], offset = 0, 
             srt = rotWord * 90, col = colors[group[i]])
        boxes[[length(boxes) + 1]] <- c(x1 - 0.5 * wid, 
                                        y1 - 0.5 * ht, wid, ht)
        isOverlaped <- FALSE
      }
      else {
        if (r > 5) {
          warning(paste(words[i], "could not be fit on page. It will not be plotted."))
          isOverlaped <- FALSE
        }
        theta <- theta + thetaStep
        if (theta > 2 * pi) 
          theta <- theta - 2 * pi
        r <- r + rStep * thetaStep/(2 * pi)
        x1 <- 0.5 + r * cos(theta)
        y1 <- 0.5 + r * sin(theta)
      }
    }
  }
  par(mar = op)
  invisible()
}

dataset_labels <- SC_abstract_8@cluster
dataset_s <- sapply(1:8 ,function(label) list( Abstracts[dataset_labels %in% label,1] ) )
names(dataset_s) <- c("Genetic Influences on Psychiatric and Behavioral Disorders",
                      "Advancements in Model Validation and Benchmarking Techniques",
                      "Travel Behavior and Environmental Impact in Transportation and Tourism",
                      "Environmental and Biological Effects of Agrochemicals",
                      "Market Segmentation Techniques and Applications",
                      "Biopharmaceutical Production Through Data-Driven Approaches",
                      "Finite Mixture Models and Their Applications",
                      "Clustering Techniques in Data Analysis")

dataset_corpus <- lapply(dataset_s, function(x) VCorpus(VectorSource( toString(x) )))
dataset_corpus_all <- dataset_corpus[[1]]
for (i in 2:8) { dataset_corpus_all <- c(dataset_corpus_all,dataset_corpus[[i]]) }
dataset_corpus_all <- tm_map(dataset_corpus_all, removePunctuation)
dataset_corpus_all <- tm_map(dataset_corpus_all, removeNumbers)
dataset_corpus_all <- tm_map(dataset_corpus_all, function(x) removeWords(x,stopwords("english")))
words_to_remove <- c("said","from","what","told","over","more","other","have","last","with","this","that","such","when","been","says","will","also","where","why","would","today")
dataset_corpus_all <- tm_map(dataset_corpus_all, removeWords, words_to_remove)
document_tm <- TermDocumentMatrix(dataset_corpus_all)
document_tm_mat <- as.matrix(document_tm)
colnames(document_tm_mat) <- names(dataset_s)
index <- as.logical(sapply(rownames(document_tm_mat), function(x) (nchar(x)>3) ))
document_tm_clean_mat_s <- document_tm_mat[index,]
comparison.cloud(document_tm_clean_mat_s,max.words=400,random.order=FALSE,c(4,0.5), title.size=0.6, use.r.layout = TRUE)
legend("bottomright", legend=names(dataset_s), cex=0.5, text.col = brewer.pal(8, "Dark2"))



#####



authors <- t(aggregate(Abstracts[,c(7:278)], by = list(SC_abstract_8@cluster), sum)[,-1])
authornames = gsub("\\.(?!\\.)", " ",  rownames(authors), perl = TRUE)
authornames = sapply(strsplit(authornames, " "), function(x) paste0(substring(head(x, 1),1,1) , ". ", tail(x, 1)))
rownames(authors) <- authornames
authors <- authors[-which(authornames == "F. Leisch"),]  
colnames(authors) <- NULL
wordcloud::comparison.cloud(authors ,max.words=Inf,random.order=FALSE,c(2.3,0.23), use.r.layout = TRUE)
legend("bottomright", legend=names(dataset_s), cex=0.5, text.col = brewer.pal(8, "Dark2"))


library(tidyverse)

Abstracts %>%
  mutate(num_of_coauthors = rowSums(Abstracts[,c(7:278)]) - 1, Clusters = factor(names(dataset_s)[SC_abstract_8@cluster], levels = names(dataset_s)) ) %>%
  ggplot() + geom_boxplot(aes(group = Clusters, y = num_of_coauthors, fill = Clusters)) + ylab("Number of co-authors")
  
num_of_coauthors = rowSums(Abstracts[,c(7:278)]) - 1

set.seed(1)
torch::torch_manual_seed(1)
SCNN_abstract_8b <- flexmix(OAI256 ~ 1 + num_of_coauthors, k = 8, model = circlus::SCauchyNN_clust_adam(EPOCHS = 200 ,LR = 0.02, free_iter = 10),
                              control = list(verbose = 1))
#saveRDS(SCNN_abstract_8b, "SCNN_abstract_8b.RDS")
table(SCNN_abstract_8b@cluster, SC_abstract_8@cluster)


Abstracts %>%
  mutate(num_of_coauthors = rowSums(Abstracts[,c(7:278)]) - 1, Clusters = factor(SCNN_abstract_8b@cluster) ) %>%
  ggplot() + geom_boxplot(aes(group = Clusters, y = num_of_coauthors, fill = Clusters)) + ylab("Number of coauthors")




dataset_labels <- SCNN_abstract_8b@cluster
dataset_s <- sapply(1:6 ,function(label) list( Abstracts[dataset_labels %in% label,1] ) )
names(dataset_s) <- 1:6
dataset_corpus <- lapply(dataset_s, function(x) VCorpus(VectorSource( toString(x) )))
dataset_corpus_all <- dataset_corpus[[1]]
for (i in 2:6) { dataset_corpus_all <- c(dataset_corpus_all,dataset_corpus[[i]]) }
dataset_corpus_all <- tm_map(dataset_corpus_all, removePunctuation)
dataset_corpus_all <- tm_map(dataset_corpus_all, removeNumbers)
dataset_corpus_all <- tm_map(dataset_corpus_all, function(x) removeWords(x,stopwords("english")))
words_to_remove <- c("said","from","what","told","over","more","other","have","last","with","this","that","such","when","been","says","will","also","where","why","would","today")
dataset_corpus_all <- tm_map(dataset_corpus_all, removeWords, words_to_remove)
document_tm <- TermDocumentMatrix(dataset_corpus_all)
document_tm_mat <- as.matrix(document_tm)
colnames(document_tm_mat) <- names(dataset_s)
index <- as.logical(sapply(rownames(document_tm_mat), function(x) (nchar(x)>3) ))
document_tm_clean_mat_s <- document_tm_mat[index,]
comparison.cloud(document_tm_clean_mat_s,max.words=200,random.order=FALSE,c(2,0.2), title.size=0.6, use.r.layout = TRUE)

document_tm_clean_mat_s[apply(document_tm_clean_mat_s>30,1, any)==1, ]



p = predict(SCNN_abstract_8b, newdata = data.frame(num_of_coauthors = num_of_coauthors[sample(129,129)]))

W = round(SCNN_abstract_8b@posterior$scaled)




howsitlooking <- function(model){
  p = predict(model)
  W = round(model@posterior$scaled)
  ll <- function(i){
    x = p[[i]]
    w = W[,i]
    w = w/sum(w)
    mu = x[,1]$mu
    rho = x[,1]$rho
    d = dim(mu)[2]
    sum(w* ((d-1)*log(1 + rho^2 - 2*rho*diag(tcrossprod(mu, OAI256)) ) -(d-1)*log(1-rho^2))  ) 
  }
  true <- sapply(seq_len(model@k), ll)
  
  simmed <- matrix(0, nrow = 10000, ncol = model@k)
  for(i in 1:10000){
    p <- predict(model, newdata = data.frame(num_of_coauthors = num_of_coauthors[sample(129,129)]))
    simmed[i,] <- sapply(seq_len(model@k), ll)
  }
  op <- par(mfrow = c(2,3))
  for(i in seq_len(model@k)){
    d <- density(simmed[,i])
    plot(d)
    pv = sum(simmed[,i] < true[i])/10000
    text(x = max(d$x) , y = max(d$y)*0.9, pv, pos = 2)
    abline(v = true[i], col = 2)
  }
  par(op)
}

SCNN_abstract_4b <- flexmix(OAI256 ~ 1 + num_of_coauthors, k = 4, model = circlus::SCauchyNN_clust_adam(EPOCHS = 200 ,LR = 0.02, free_iter = 10),
                            control = list(verbose = 1))
table(SCNN_abstract_4b@cluster, SC_abstract_8@cluster)

