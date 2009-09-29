# Copyright (C) 2009 
# Kim-Anh Lê Cao, ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# and Queensland Facility for Advanced Bioinformatics, The University of Queensland, Australia
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.




kmeans.init=function(data.cont){

continue = TRUE
cl <- kmeans(data.cont, 2, nstart = 5)

if(is.null(nrow(data.cont[cl$cluster==1,])) | is.null(nrow(data.cont[cl$cluster==2,]))){continue=FALSE}else{if(nrow(data.cont[cl$cluster==1,]) <= 1 | nrow(data.cont[cl$cluster==2,]) <= 1){continue = FALSE}}

if(continue == TRUE){
clust1 = data.cont[cl$cluster==1,]
clust2 = data.cont[cl$cluster==2,]

# ----3. store outputs -------------------------------

prop.kmeans = t(cl$size/sum(cl$size))
means.kmeans = cl$centers
var.kmeans = array(0, dim=c(ncol(data.cont),ncol(data.cont),2 ))
var.kmeans[,,1] = var(clust1)
var.kmeans[,,2] = var(clust2)

} #end continue = TRUE

return(invisible(list(
prop.kmeans = prop.kmeans,
means.kmeans = means.kmeans,
var.kmeans = var.kmeans,
continue = continue
)))

}

