# ---2/05, updated 11/05
#included in total code 16/06, added condition again and regularization of var/cov matrices

#use K-means to initialize the parameters: proportions and variances

############--------- initialization with K-means ------------------------


# ----1. read the continuous data ---------------------

#-----2. launch k means ------------------------------


kmeansME=function(data.cont.keep){

continue = TRUE
cl <- kmeans(data.cont.keep, 2, nstart = 5)

if(is.null(nrow(data.cont.keep[cl$cluster==1,])) | is.null(nrow(data.cont.keep[cl$cluster==2,]))){continue=FALSE}else{if(nrow(data.cont.keep[cl$cluster==1,]) <= 1 | nrow(data.cont.keep[cl$cluster==2,]) <= 1){continue = FALSE}}

if(continue == TRUE){
clust1 = data.cont.keep[cl$cluster==1,]
clust2 = data.cont.keep[cl$cluster==2,]

# ----3. store outputs -------------------------------

prop.kmeans = t(cl$size/sum(cl$size))
means.kmeans = cl$centers
var.kmeans = array(0, dim=c(ncol(data.cont.keep),ncol(data.cont.keep),2 ))
var.kmeans[,,1] = var(clust1)
var.kmeans[,,2] = var(clust2)

} #end continue = TRUE

return(invisible(list(
prop.kmeans = prop.kmeans,
means.kmeans = means.kmeans,
var.kmeans = var.kmeans
)))

}

