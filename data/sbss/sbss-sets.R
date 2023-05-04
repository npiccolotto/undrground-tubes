results[[1]]$s

all.s <- do.call(cbind, purrr::map(results, function(x){return(x$s)}))

k <- ncol(all.s)

dist_euclid <- function(x,y) {
  return(sqrt( sum((x-y))^2 ))
}

sum(c(0,0,1) - c(0,2,0))

dist_euclid(c(0,0,1), c(0,2,0))

dist(matrix(c(0,0,1,0,2,0), nrow=2, byrow=T))[1]

D.s <- outer(1:k, 1:k, Vectorize(function(i,j) {
  return(dist( rbind(all.s[,i], all.s[,j]) )[1])
}))

coords <-
  sp::spTransform(sp::SpatialPoints(results[[1]]$coords, CRS_FLAT), WGS84)
D.points <- 1/as.matrix(sp::spDists(coords))^2
diag(D.points) <- 0

feat.moransi <- apply(all.s, 2, function(x){return(moransi(x, D.points))})
breaks.moransi <- c( BAMMtools::getJenksBreaks(feat.moransi, 6))
sets.moransi <- cut(feat.moransi, breaks=breaks.moransi, labels=1:5, include.lowest = T)

feat.kurtosis <- apply(all.s, 2, function(x){return(moments::kurtosis(x))})
breaks.kurtosis <- c( BAMMtools::getJenksBreaks(feat.kurtosis, 6))
sets.kurtosis <- cut(feat.kurtosis, breaks=breaks.kurtosis, labels=1:5, include.lowest = T)

feat.skewness <- apply(all.s, 2, function(x){return(moments::skewness(x))})
breaks.skewness <- c( BAMMtools::getJenksBreaks(feat.skewness, 6))
sets.skewness <- cut(feat.skewness, breaks=breaks.skewness, labels=1:5, include.lowest = T)



colnames(all.s) <- paste("Map_",1:k, sep='')

element.names <- paste("Map_",1:k, sep='')


as_json <- list(
  E = element.names,
  EA = D.s,
  SR = list(
    moransi = list(
      bin = sets.moransi,
      breaks = breaks.moransi,
      values = feat.moransi
    ),
    kurtosis = list(
      bin = sets.kurtosis,
      breaks = breaks.kurtosis,
      values = feat.kurtosis
    ),
    skewness = list(
      bin = sets.skewness,
      breaks = breaks.skewness,
      values = feat.skewness
    )
  )
)


setwd('/Users/npiccolotto/Desktop')
jsonlite::write_json(as_json, 'ensemble_set.json')

bbox <- sp::bbox(P.spatial)
gmap <- ggmap::get_map(location=bbox, source='osm', zoom = 8, color = 'bw', force = TRUE)
coords <-
  sp::spTransform(sp::SpatialPoints(results[[1]]$coords, CRS_FLAT), WGS84)@coords

bbox_flat <- sp::bbox(P)
bbox_ar <- bbox_flat[,2]-bbox_flat[,1]
bbox_ar <- bbox_ar[1]/bbox_ar[2]
i<-1
for (i in 1:k) {
  x <-
    sp::SpatialPointsDataFrame(coords = coords, data = as.data.frame(all.s[,i]))
  
  d <-
    as.data.frame(cbind(
      x@coords,
      abs(x@data),
      x@data < 0
    ))
  
  colnames(d) <-
    c('longitude',
      'latitude',
      'c',
      'is_neg')
  
  jpeg(paste(element.names[i], 'jpg', sep = '.'),
       width = 1500*bbox_ar,
       height = 1500)
  
  plot(ggmap::ggmap(gmap) +
         ggplot2::geom_point(
           data = d,
           stroke = 0,
           shape=23,
           alpha=0.7,
           ggplot2::aes(
             x = longitude,
             y = latitude,
             size = c,
             fill=is_neg
           )
         ) +
         ggplot2::scale_fill_manual(values=c('red','blue'))+
         ggplot2::scale_size_continuous(range = c(0,48), limits = c(0,abs(max(all.s)))) +
         ggplot2::theme(
           legend.position = 'none',
           axis.ticks = ggplot2::element_blank(),
           axis.text = ggplot2::element_blank(),
           axis.title = ggplot2::element_blank()
         ) 
  )
  
  dev.off()
}

