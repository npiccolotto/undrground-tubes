setwd('/Users/npiccolotto/Projects/cvast/bssvis/sbss-app-gh/server/app/data')
library(ggplot2)
library(ggmap)

util_get_bbox_polygon <- function(bbox, margin = 1.05, proj4string=WGS84) {
  dataset.extent <-
    (bbox[, 2] - bbox[, 1])
  dataset.center <- dataset.extent / 2 + bbox[, 1]
  sidelength <- dataset.extent * margin
  dataset.bbox <- matrix(
    c(
      dataset.center - c(sidelength[1] / 2, sidelength[2] / 2),
      dataset.center + c(sidelength[1] / 2,-sidelength[2] / 2),
      dataset.center + c(sidelength[1] / 2, sidelength[2] / 2),
      dataset.center + c(-sidelength[1] / 2, sidelength[2] / 2)
    ),
    ncol = 2,
    byrow = T
  )
  colnames(dataset.bbox) <- c('longitude','latitude')
  dataset.extent <-
    c(max(dataset.bbox[, 1]) - min(dataset.bbox[, 1]),
      max(dataset.bbox[, 2]) - min(dataset.bbox[, 2]))
  dataset.bbox <-
    sp::SpatialPolygons(list(sp::Polygons(list(
      sp::Polygon(dataset.bbox, hole = F)
    ), ID = 'bbox')), proj4string = proj4string)
  return(dataset.bbox)
}

moss <- read.csv('kola.csv')
moss

WGS84 <- sp::CRS('+proj=longlat +datum=WGS84 +no_defs')
UTM <- sp::CRS('+proj=utm +zone=35 +datum=WGS84 +units=m +no_defs')
WEB_MERC <-
  sp::CRS(
    '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs'
  )

coords <- moss[,1:2]
coords.spatial <- sp::SpatialPoints(coords, proj4string = WGS84)
coords.flat <- sp::spTransform(coords.spatial, UTM)
ilrBase <- compositions::ilrBase(moss[,3:ncol(moss)])
moss.ilr <- log(as.matrix(moss[,3:ncol(moss)]), base=10) %*% ilrBase
moss.ilr.spatial <- sp::SpatialPointsDataFrame(coords=coords,data=as.data.frame(moss.ilr))

kernels <- c(25000, 50000, 75000, 100000)
results <- list()

i <- 1

for (kernel in kernels) {
  results[[i]]<-SpatialBSS::sbss(x=moss.ilr,coords=coords.flat@coords,kernel_type='ball',kernel_parameters = c(kernel))
  i <- i+1
}

components <- results[[1]]$s
for (j in seq(2,i-1)) {
  components <- cbind(components, results[[j]]$s)
}
colnames(components) <- 1:120

loadings <- results[[1]]$w %*% t(ilrBase)
for (j in seq(2,i-1)){
  loadings <- rbind(loadings, results[[j]]$w %*% t(ilrBase))
}

elements <- colnames(moss[,3:ncol(moss)])
colnames(loadings) <- elements
rownames(loadings) <- 1:120

sp::spDists(as.matrix(coords))
coords

D.points <- 1/as.matrix(sp::spDists(as.matrix(coords)))^2
diag(D.points) <- 0

moransi <- function(x) {
  return(ape::Moran.I(x, D.points)$observed)
}

x <- setNames(-3:3, c('A', 'B', 'C', 'D', 'E', 'F', 'G'))
"F" %in% names(sort(abs(x), decreasing=T)[1:3])

make_feats <- function(M, feats = c()) {
  n.feats <- length(names(feats))
  # if there's at least 1 feature, apply it per col, categorize per-col values
  if (n.feats > 1) {
    M.feats <- list()
    for (f in names(feats)) {
      M.feat <- apply(M, 2, feats[[f]])
      M.break <- BAMMtools::getJenksBreaks(M.feat, 4)
      M.label <- cut(M.feat, breaks=M.break, include.lowest=T)
      M.label.int <- cut(M.feat, breaks=M.break, include.lowest=T, labels=1:3)
      levels(M.label.int) <- c(
        paste0(f, ': ', 'low'),
        paste0(f, ': ', 'mid'),
        paste0(f, ': ', 'high')
      )
      M.feats[[f]] <- list(
        'raw'=unlist(M.feat),
        'breaks'=M.break,
        'labels-breaks'=M.label,
        'labels'=M.label.int
      )
    }
    return(M.feats)
  } else {
    # otherwise do no feature derivation per col and handle columns directly
    # TODO do it differently
    # since analysts look at rows, not cols of the W matrix
    # let's classify the top/bottom N as high/low loadings
    # e.g., top 3 loadings ni, s, ag -> high ni, high s, high ag. bot 3 loadings Mg, Zn,  U -> low these. rest "unspecific".
    M.feats <- list()
    
    M.feat <- apply(M, 1, function(row) {
      top3 <- names(sort(row,decreasing =F)[1:3])
      bot3 <- names(sort(row,decreasing =T)[1:3])
      M.label.int <- purrr::map(colnames(M), function(cname) {
        if (cname %in% top3) {
          return(paste0(cname, ': high'))
        }
        if (cname %in% bot3) {
          return(paste0(cname, ': low'))
        }
        return(paste0(cname, ': mid'))
      })
      return(M.label.int)
    })
    
    for (col in colnames(M)) {
        M.feats[[col]] <- list(
          labels=purrr::map(M.feat, function(complabels) {
            if (paste0(col, ': high') %in% complabels) {
              return(paste0(col, ': high'))
            }
            if (paste0(col, ': mid') %in% complabels) {
              return(paste0(col, ': mid'))
            }
            if (paste0(col, ': low') %in% complabels) {
              return(paste0(col, ': low'))
            }
          })
        )
    }
    
    return(M.feats)
  }
}

loads.feats <- make_feats(loadings)
comps.feats <- make_feats(components, list('kurtosis'=moments::kurtosis, 'skewness'=moments::skewness, 'moransi'=moransi))

kernel_labels <- c(rep(25, 30), rep(50, 30), rep(75, 30), rep(100,30))
comps.feats[['kernel']] <- list(
  'raw'=kernel_labels,
  'breaks'=c(),
  'labels-breaks'=c(),
  'labels'=factor(paste('kernel:', kernel_labels))
)

# Image plotting

get_breaks <- function(x) {
  return(quantile(unlist(x), c(0, 0.05, 0.25, 0.75, 0.95, 1)))
}

fin_nor_rus <- spMaps::getSpMaps(countries=NULL, states=c('FIN','NOR','RUS', 'SWE'))
expanded_bbox <- coords.spatial@bbox * matrix(c(0.95, 1.05, 0.995, 1.005), ncol=2, byrow=T)
fin_nor_rus_map <- raster::intersect(fin_nor_rus, coords.spatial@bbox)

bbox.flat <- sp::spTransform(coords.spatial, WEB_MERC)@bbox
bbox_ar <- bbox.flat[,2]-bbox.flat[,1]
bbox_ar <- bbox_ar[1]/bbox_ar[2]

comps <- sp::SpatialPointsDataFrame(coords=coords, data=as.data.frame(components), proj4string = WGS84)

path_rel <- 'data/sbss/img'
img_path <-  '/Users/npiccolotto/Projects/cvast/bssvis/ensemble-set-rendering/'
setwd(paste0(img_path, path_rel))
for (j in 1:120) {
  comp.interp <- interp::interp(comps, y=NULL, z=paste(j),output='grid', nx=40, ny=40)
  comp.df <- as.data.frame(cbind(comp.interp@coords, comp.interp@data))
  colnames(comp.df) <- c('x', 'y', 'z')
  
  jpeg(paste(paste(j), 'jpg', sep = '.'),
       width = 500*bbox_ar,
       height = 500)
  
  plot(
    ggplot2::ggplot(comp.df) + 
         ggplot2::stat_contour_filled(breaks=get_breaks(comp.interp@data), mapping=ggplot2::aes(x=x,y=y,z=z)) +
         ggplot2::scale_fill_brewer(type='div', palette='RdBu',direction = -1) +
         ggplot2::geom_point(as.data.frame(coords), mapping=ggplot2::aes(x=longitude,y=latitude), size=.2) +
         ggplot2::geom_path(ggplot2::fortify(fin_nor_rus_map), mapping=ggplot2::aes(x = long, y = lat, group = group), colour='black') +
         ggplot2::coord_map() +
         ggplot2::theme_void() +
         ggplot2::theme(
           legend.position='none'
         )
  )
  
  dev.off()
}

comps.df <- as.data.frame(cbind(coords, abs(comps$`10`)))
colnames(comps.df) <- c('x','y','z')

# TODO experiment with other plots
ggplot2::ggplot(comp.df) + 
  #ggplot2::stat_contour_filled(breaks=get_breaks(comp.interp@data), mapping=ggplot2::aes(x=x,y=y,z=z)) +
  #ggplot2::scale_fill_brewer(type='div', palette='RdBu',direction = -1) +
  ggplot2::geom_point(comps.df, mapping=ggplot2::aes(x=x,y=y, size=z)) +
  ggplot2::scale_size_continuous(range=c(0.01,3)) +
  ggplot2::geom_path(ggplot2::fortify(fin_nor_rus_map), mapping=ggplot2::aes(x = long, y = lat, group = group), colour='black') +
  ggplot2::coord_map() +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position='none'
  )



# Components distance matrix

D.s <- outer(1:120, 1:120, Vectorize(function(i,j) {
  return(1-abs(cor(components[,i],components[,j])))
}))

# Write all the things

# Edit these as necessary, rest should be automatic
#S <- c('kernel: 25', 'kernel: 100')
#SC <- c('#005824', '#ef6548')

S <- c('Ni: low', 'Ni: high')
SC <- c('#2166ac', '#b2182b')

cl.feats <- c(comps.feats,loads.feats)
# make binary matrix of set memberships
SM <- matrix(data=0, ncol=length(S), nrow=120)
for (j in 1:length(S)) {
  set <- S[j]
  category <- strsplit(set, ': ')[[1]][1]
  SM[,j] <- as.integer(cl.feats[[category]]$labels == set)
}
SA <- as.matrix(dist(SM, method='binary'))
SR <- list()
for (j in 1:120) {
  in_sets_idx <- which(SM[j,] == 1)
  SR[[j]] <- unlist(purrr::map(in_sets_idx, function(idx) {return(S[idx])}))
}
E <- paste('SBSS',1:120)
EA <- as.matrix(D.s)

as_json <- list(
  'glyph_ids'=paste(paste(path_rel, 1:120, sep='/'),'jpg',sep='.'),
  'E'=E,
  'S'=S,
  'SR'=SR,
  'SC'=SC,
  'EA'=EA,
  'SA'=SA,
)

jsonlite::write_json(as_json, paste(img_path, 'moss.json', sep='/'), simplifyVector = F)

