library("maSigPro")
library("igraph")
library("RCytoscape")
library("pvclust")
library("multiClust")
dir <- setwd("~/directory")

alldesign <- read.table("design", header = TRUE)
alldata <- read.table("allcounts", header = TRUE)
expdesign <- alldesign[27:52,]
dataset <- alldata[,27:52]
dataset <- dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]  #Removing rows with 0 counts in all parallels
design <- make.design.matrix(expdesign, degree = 2)
fit <- p.vector(dataset, design, Q = 0.05, MT.adjust = "BH", min.obs = 3, counts=TRUE, theta=10)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
dataset <- subset(dataset, !(row.names(dataset) %in% colnames(tstep$influ.info))) # Removing outliers

#-----Pulling outliers-----
int <- intersect(sigs$summary[,2],sigs$summary[,3]) # Intersection
left <- setdiff(sigs$summary[,2],sigs$summary[,3])  # Only in 1st
right <- setdiff(sigs$summary[,3],sigs$summary[,2]) # Only in 2nd

dataset_l <- subset(dataset, row.names(dataset) %in% left)
fit_l <- p.vector(dataset_l, design, Q = 0.05, MT.adjust = "BH", min.obs = 3, counts=TRUE, theta=10)
tstep_l <- T.fit(fit_l, step.method = "backward", alfa = 0.05)

dataset <- subset(dataset, row.names(dataset) %in% right)
fit_r <- p.vector(dataset_r, design, Q = 0.05, MT.adjust = "BH", min.obs = 3, counts=TRUE, theta=10)
tstep_r <- T.fit(fit_r, step.method = "backward", alfa = 0.05)
#-----Plotting-----
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
suma2Venn(sigs$summary[, c(2:3)])
plot <- see.genes(data, main = "Exs vs Axi", show.fit = T,
                  dis =design$dis, cluster.method="hclust", cluster.data = 1, k = 9)
PlotGroups(dataset,  edesign = expdesign)

matrix <- as.matrix(read.table("../../19.04.16_clusters/MSP/base100corr", header=TRUE, row.names = 1, as.is=TRUE))
candidate <- readLines("../../19.10.15/axex_S_2")
step1 <- subset(matrix, colnames(matrix) %in% candidate)
step2 <- t(step1)
step3 <- subset(step2, rownames(step2) %in% candidate)

mincond <- c(rep("AxSS", 3), rep("AxMS", 3), rep("AxLS", 3))
mincounts <- subset(ncounts, row.names(ncounts) %in% candidate)
minnorm <- sweep(mincounts, 1, unlist(mincounts[,53]),  "/") # normalizing by baseMean value
mincounts <- t(mincounts)
mincounts <- subset(mincounts, colnames(ncounts) %in% mincond)
mincor <- cor(mincounts)
write.table(mincor,     file = "../../oxi_tf.csv", col.names = TRUE, row.names = TRUE, quote = FALSE)

mm <- matrix(1, dim(step3), dim(step3))
distance <- mm - abs(step3)
graph <- graph.adjacency(step3, weighted=TRUE, mode="upper")
sgraph <- simplify(graph, remove.multiple = F, remove.loops = T) 
plot(sgraph, edge.arrow.size=.4,vertex.label=NA, vertex.size=7)
layout <- layout.reingold.tilford(sgraph)
plot(sgraph, layout = layout)
write.graph(graph = sgraph, file = "../../graph.gml", format = "gml")
#---Manual---
#--As in maSigPro--
clusterdata <- sigs$sig.genes$ExsertavsAxillaris[[1]]
dcorrel <- matrix(rep(1, nrow(clusterdata)^2), 
                  nrow(clusterdata), nrow(clusterdata)) - cor(t(clusterdata), 
                  use = "pairwise.complete.obs")
clust <- hclust(as.dist(dcorrel), method = "ward.D2")
plot(clust, labels = NULL, xlab = NULL)
#--With pvclust - sample outliers--
result <- pvclust(ncounts, nboot=1000, method.dist="correlation", method.hclust = "ward.D2")
plot(result)
pvrect(result, alpha=0.95)
seplot(result)
#--Approximating cluster numbers--
cluster_num <- number_clusters(data.exp=ncounts, Fixed=NULL,
                               gap_statistic=TRUE)
#-----Writing data to files-----
write.table(colnames(tstep$influ.info), file = "../../outliers", col.names = FALSE, row.names = TRUE, quote = FALSE)
write.table(sigs$summary, file = "../../genelist", col.names = FALSE, row.names = TRUE, quote = FALSE)
write.table(vstc, file = "../../vstc.csv", col.names = TRUE, row.names = TRUE, quote = FALSE, sep = ",")

#---Edited see.genes - absolute correlation value---
"sea.genes" <-
  function (data, edesign = data$edesign, time.col = 1, repl.col = 2, 
            group.cols = c(3:ncol(edesign)), names.groups = colnames(edesign)[3:ncol(edesign)], 
            cluster.data = 1, groups.vector = data$groups.vector, k = 9, m = 1.45, 
            cluster.method = "hclust", distance = "cor", agglo.method = "ward.D", 
            show.fit = FALSE, dis = NULL, step.method = "backward", min.obs = 3, 
            alfa = 0.05, nvar.correction = FALSE, show.lines = TRUE, iter.max = 500, 
            summary.mode = "median", color.mode = "rainbow", cexlab = 1, legend = TRUE, 
            newX11 = TRUE,  ylim = NULL, main = NULL, ...) 
  {
    time = edesign[, time.col]
    repvect = edesign[, repl.col]
    groups = edesign[, group.cols]
    narrays <- length(time)
    if (!is.null(dim(data))) {
      dat <- as.data.frame(data)
      clusterdata <- data
    }
    else {
      clusterdata <- data[[cluster.data]]
      dat <- as.data.frame(data$sig.profiles)
    }
    clusterdata <- clusterdata
    if (nrow(dat) > 1) {
      dat <- as.data.frame(dat[, (ncol(dat) - length(time) + 
                                    1):ncol(dat)])
      count.na <- function(x) length(x[is.na(x)])
      NAs <- apply(as.matrix(dat), 1, count.na)
      count.noNa <- function(x) (length(x) - length(x[is.na(x)]))
      dat <- dat[which(apply(as.matrix(dat), 1, count.noNa) >= 
                         2), ]
    }
    else {
      NAs <- 1
    }
    kdata <- NULL
    out <- TRUE
    
    if (nrow(dat) > 1) {
      if (cluster.data != 1 || cluster.data != "sig.profiles") {
        if (any(is.na(clusterdata))) 
          clusterdata[is.na(clusterdata)] <- 0
      }
      else if (is.na(all(dist(clusterdata) > 0)) || (cluster.method == 
                                                     "kmeans" & any(is.na(clusterdata))) || (distance == 
                                                                                             "cor" & any(sd(t(clusterdata), na.rm = TRUE) == 0))) {
        if (!is.null(kdata)) {
          clusterdata <- kdata
        }
        else {
          clusterdata <- NULL
        }
      }
      clusterdata <- clusterdata
      if (!is.null(clusterdata)) {
        k <- min(k, nrow(dat), na.rm = TRUE)
        if (cluster.method == "hclust") {
          if (distance == "cor") {
            dcorrel <- matrix(rep(1, nrow(clusterdata)^2), 
                              nrow(clusterdata), nrow(clusterdata)) - abs(cor(t(clusterdata)), #Edit here
                                                                          use = "pairwise.complete.obs")
            clust <- hclust(as.dist(dcorrel), method = agglo.method)
            c.algo.used = paste(cluster.method, "cor", 
                                agglo.method, sep = "_")
          }
          else {
            clust <- hclust(dist(clusterdata, method = distance), 
                            method = agglo.method)
            c.algo.used = paste(cluster.method, distance, 
                                agglo.method, sep = "_")
          }
          cut <- cutree(clust, k = k)
        }
        else if (cluster.method == "kmeans") {
          cut <- kmeans(clusterdata, k, iter.max)$cluster
          c.algo.used = paste("kmeans", k, iter.max, sep = "_")
        }
        else if (cluster.method == "mfuzz") { 
          n<-dim(clusterdata)[2]
          clusterdata[is.na(clusterdata)]<-0
          temp <- tempfile()
          write.table(clusterdata, temp, quote = FALSE, sep = "\t", row.names =TRUE, col.names = TRUE)
          signif <- readExpressionSet(temp)
          cl <- mfuzz(signif, c = k, m = m)
          clus<-acore(signif,cl=cl,min.acore=(1/k))
          for(i in 1:k){
            clus[[i]]<-transform(clus[[i]],cluster= i )
          }
          cut0<-clus[[1]][,c(1,3)]
          for(i in 2:k){
            cut0<-rbind(cut0, clus[[i]][,c(1,3)])
          }
          cut<-transform(clusterdata, name="")
          cut<-transform(cut, cluster=0)
          cut<-cut[,c(n+1,n+2)]
          cut[,1]<-rownames(cut)
          for(i in 1:dim(clusterdata)[1]){
            cut[i,2]<-cut0[cut[i,1],2]
          }
          cut<-cut[,2]
          c.algo.used = paste("mfuzz", k, m, sep = "_")
        }
        else stop("Invalid cluster algorithm")
        if (newX11) 
          X11()
        groups <- as.matrix(groups)
        colnames(groups) <- names.groups
        if (k <= 4) 
          par(mfrow = c(2, 2))
        else if (k <= 6) 
          par(mfrow = c(3, 2))
        else if (k > 6) 
          par(mfrow = c(3, 3))
        for (i in 1:(k)) {
          PlotProfiles(data = dat[cut == i, ], repvect = repvect, 
                       main = i, ylim = ylim, color.mode = color.mode, 
                       cond = rownames(edesign), ...)
        }
        if (newX11) 
          X11()
        if (k <= 4) {
          par(mfrow = c(2, 2))
          cexlab = 0.6
        }
        else if (k <= 6) {
          par(mfrow = c(3, 2))
          cexlab = 0.6
        }
        else if (k > 6) {
          par(mfrow = c(3, 3))
          cexlab = 0.35
        }
        for (j in 1:(k)) {
          PlotGroups(data = dat[cut == j, ], show.fit = show.fit, 
                     dis = dis, step.method = step.method, min.obs = min.obs, 
                     alfa = alfa, nvar.correction = nvar.correction, show.lines = show.lines, time = time, 
                     groups = groups, repvect = repvect, summary.mode = summary.mode, 
                     xlab = "time", main = paste("Cluster", j, sep = " "), 
                     ylim = ylim, cexlab = cexlab, legend = legend, 
                     groups.vector = groups.vector, ...)
        }
      }
      else {
        print("warning: impossible to compute hierarchical clustering")
        c.algo.used <- NULL
        cut <- 1
      }
    }
    else if (nrow(dat) == 1) {
      if (newX11) 
        X11()
      PlotProfiles(data = dat, repvect = repvect, main = NULL, 
                   ylim = ylim, color.mode = color.mode, cond = rownames(edesign), 
                   ...)
      if (newX11) 
        X11()
      PlotGroups(data = dat, show.fit = show.fit, dis = dis, 
                 step.method = step.method, min.obs = min.obs, alfa = alfa, nvar.correction = nvar.correction,
                 show.lines = show.lines, time = time, groups = groups, 
                 repvect = repvect, summary.mode = summary.mode, xlab = "time", 
                 main = main, ylim = ylim, cexlab = cexlab, legend = legend, 
                 groups.vector = groups.vector, ...)
      c.algo.used <- NULL
      cut <- 1
    }
    else {
      print("warning: NULL data. No visualization possible")
      c.algo.used <- NULL
      cut <- NULL
    }
    OUTPUT <- list(cut, c.algo.used, groups)
    names(OUTPUT) <- c("cut", "cluster.algorithm.used", "groups")
    OUTPUT
  }
