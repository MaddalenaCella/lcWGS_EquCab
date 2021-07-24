#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 03-07-2021
# Last Modified: 03-07-2021
# Description: script that plots PCA graphs


library(dplyr)
library(ggplot2)
library(tools)
library(RColorBrewer)
library(ggforce)
require(gghighlight)
require(ggpubr)
require(scales)
library(grid)

#devtools::install_github("thomasp85/ggforce", ref = '4008a2e')
facet_zoom2 <- function(x, y, xy, zoom.data, xlim = NULL, ylim = NULL, 
                        split = FALSE, horizontal = TRUE, zoom.size = 2, 
                        show.area = TRUE, shrink = TRUE) {
  x <- if (missing(x)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(x)
  y <- if (missing(y)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(y)
  zoom.data <- if (missing(zoom.data)) NULL else lazyeval::lazy(zoom.data)
  if (is.null(x) && is.null(y) && is.null(xlim) && is.null(ylim)) {
    stop("Either x- or y-zoom must be given", call. = FALSE)
  }
  if (!is.null(xlim)) x <- NULL
  if (!is.null(ylim)) y <- NULL
  ggproto(NULL, FacetZoom2,
          shrink = shrink,
          params = list(
            x = x, y = y, xlim = xlim, ylim = ylim, split = split, zoom.data = zoom.data,
            zoom.size = zoom.size, show.area = show.area,
            horizontal = horizontal
          )
  )
}

FacetZoom2 <- ggproto(
  "FacetZoom2",
  ggforce::FacetZoom,
  
  compute_layout = function(data, params) {
    layout <- rbind( # has both x & y dimension
      data.frame(name = 'orig', SCALE_X = 1L, SCALE_Y = 1L),
      data.frame(name = 'x', SCALE_X = 2L, SCALE_Y = 1L),
      data.frame(name = 'y', SCALE_X = 1L, SCALE_Y = 2L),
      data.frame(name = 'full', SCALE_X = 2L, SCALE_Y = 2L),
      data.frame(name = 'orig_true', SCALE_X = 1L, SCALE_Y = 1L),
      data.frame(name = 'zoom_true', SCALE_X = 1L, SCALE_Y = 1L)
    )
    if (is.null(params$y) && is.null(params$ylim)) { # no y dimension
      layout <- layout[c(1,2, 5:6),]
    } else if (is.null(params$x) && is.null(params$xlim)) { # no x dimension
      layout <- layout[c(1,3, 5:6),]
    }
    layout$PANEL <- seq_len(nrow(layout))
    layout
  },
  
  draw_panels = function(panels, layout, x_scales, y_scales, ranges, coord,
                         data, theme, params) {
    
    if (is.null(params$x) && is.null(params$xlim)) {
      params$horizontal <- TRUE
    } else if (is.null(params$y) && is.null(params$ylim)) {
      params$horizontal <- FALSE
    }
    if (is.null(theme[['zoom']])) {
      theme$zoom <- theme$strip.background
    }
    if (is.null(theme$zoom.x)) {
      theme$zoom.x <- theme$zoom
    }
    if (is.null(theme$zoom.y)) {
      theme$zoom.y <- theme$zoom
    }
    axes <- render_axes(ranges, ranges, coord, theme, FALSE)
    panelGrobs <- ggforce:::create_panels(panels, axes$x, axes$y)
    panelGrobs <- panelGrobs[seq_len(length(panelGrobs) - 2)]
    if ('full' %in% layout$name && !params$split) {
      panelGrobs <- panelGrobs[c(1, 4)]
    }
    
    # changed coordinates in indicator / lines to zoom from 
    # the opposite horizontal direction
    if ('y' %in% layout$name) {
      if (!inherits(theme$zoom.y, 'element_blank')) {
        zoom_prop <- scales::rescale(
          y_scales[[2]]$dimension(ggforce:::expansion(y_scales[[2]])),
          from = y_scales[[1]]$dimension(ggforce:::expansion(y_scales[[1]])))
        indicator <- polygonGrob(
          x = c(0, 0, 1, 1), # was x = c(1, 1, 0, 0), 
          y = c(zoom_prop, 1, 0), 
          gp = gpar(col = NA, fill = alpha(theme$zoom.y$fill, 0.5)))
        lines <- segmentsGrob(
          x0 = c(1, 1), x1 = c(0, 0), # was x0 = c(0, 0), x1 = c(1, 1)
          y0 = c(0, 1), y1 = zoom_prop,
          gp = gpar(col = theme$zoom.y$colour,
                    lty = theme$zoom.y$linetype,
                    lwd = theme$zoom.y$size,
                    lineend = 'round'))
        indicator_h <- grobTree(indicator, lines)
      } else {
        indicator_h <- zeroGrob()
      }
    }
    
    if ('x' %in% layout$name) {
      if (!inherits(theme$zoom.x, 'element_blank')) {
        zoom_prop <- scales::rescale(x_scales[[2]]$dimension(ggforce:::expansion(x_scales[[2]])),
                                     from = x_scales[[1]]$dimension(ggforce:::expansion(x_scales[[1]])))
        indicator <- polygonGrob(c(zoom_prop, 1, 0), c(1, 1, 0, 0), 
                                 gp = gpar(col = NA, fill = alpha(theme$zoom.x$fill, 0.5)))
        lines <- segmentsGrob(x0 = c(0, 1), y0 = c(0, 0), x1 = zoom_prop, y1 = c(1, 1), 
                              gp = gpar(col = theme$zoom.x$colour,
                                        lty = theme$zoom.x$linetype,
                                        lwd = theme$zoom.x$size,
                                        lineend = 'round'))
        indicator_v <- grobTree(indicator, lines)
      } else {
        indicator_v <- zeroGrob()
      }
    }
    
    if ('full' %in% layout$name && params$split) {
      space.x <- theme$panel.spacing.x
      if (is.null(space.x)) space.x <- theme$panel.spacing
      space.x <- unit(5 * as.numeric(convertUnit(space.x, 'cm')), 'cm')
      space.y <- theme$panel.spacing.y
      if (is.null(space.y)) space.y <- theme$panel.spacing
      space.y <- unit(5 * as.numeric(convertUnit(space.y, 'cm')), 'cm')
      
      # change horizontal order of panels from [zoom, original] to [original, zoom]
      # final <- gtable::gtable_add_cols(panelGrobs[[3]], space.x)
      # final <- cbind(final, panelGrobs[[1]], size = 'first')
      # final_tmp <- gtable::gtable_add_cols(panelGrobs[[4]], space.x)
      # final_tmp <- cbind(final_tmp, panelGrobs[[2]], size = 'first')
      final <- gtable::gtable_add_cols(panelGrobs[[1]], space.x)
      final <- cbind(final, panelGrobs[[3]], size = 'first')
      final_tmp <- gtable::gtable_add_cols(panelGrobs[[2]], space.x)
      final_tmp <- cbind(final_tmp, panelGrobs[[4]], size = 'first')
      
      final <- gtable::gtable_add_rows(final, space.y)
      final <- rbind(final, final_tmp, size = 'first')
      final <- gtable::gtable_add_grob(final, list(indicator_h, indicator_h),
                                       c(2, 6), 3, c(2, 6), 5,
                                       z = -Inf, name = "zoom-indicator")
      final <- gtable::gtable_add_grob(final, list(indicator_v, indicator_v), 
                                       3, c(2, 6), 5, 
                                       z = -Inf, name = "zoom-indicator")
      heights <- unit.c(
        unit(max_height(list(axes$x[[1]]$top, axes$x[[3]]$top)), 'cm'),
        unit(1, 'null'),
        unit(max_height(list(axes$x[[1]]$bottom, axes$x[[3]]$bottom)), 'cm'),
        space.y,
        unit(max_height(list(axes$x[[2]]$top, axes$x[[4]]$top)), 'cm'),
        unit(params$zoom.size, 'null'),
        unit(max_height(list(axes$x[[2]]$bottom, axes$x[[4]]$bottom)), 'cm')
      )
      
      # swop panel width specifications according to the new horizontal order
      widths <- unit.c(
        # unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
        # unit(params$zoom.size, 'null'),
        # unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm'),
        # space.x,
        # unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
        # unit(1, 'null'),
        # unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')        
        unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
        unit(1, 'null'),
        unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm'),
        space.x,
        unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
        unit(params$zoom.size, 'null'),
        unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm')
        
      )
      final$heights <- heights
      final$widths <- widths
    } else {
      if (params$horizontal) {
        space <- theme$panel.spacing.x
        if (is.null(space)) space <- theme$panel.spacing
        space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
        heights <- unit.c(
          unit(max_height(list(axes$x[[1]]$top, axes$x[[2]]$top)), 'cm'),
          unit(1, 'null'),
          unit(max_height(list(axes$x[[1]]$bottom, axes$x[[2]]$bottom)), 'cm')
        )
        
        # change horizontal order of panels from [zoom, original] to [original, zoom]
        # first <- gtable::gtable_add_cols(panelGrobs[[2]], space)
        # first <- cbind(final, panelGrobs[[1]], size = 'first')
        final <- gtable::gtable_add_cols(panelGrobs[[1]], space) 
        final <- cbind(final, panelGrobs[[2]], size = "first") 
        
        final$heights <- heights
        
        # swop panel width specifications according to the new horizontal order
        # unit(c(params$zoom.size, 1), 'null')
        final$widths[panel_cols(final)$l] <- unit(c(1, params$zoom.size), 'null') 
        
        final <- gtable::gtable_add_grob(final, indicator_h, 2, 3, 2, 5, 
                                         z = -Inf, name = "zoom-indicator")
      } else {
        space <- theme$panel.spacing.y
        if (is.null(space)) space <- theme$panel.spacing
        space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
        widths <- unit.c(
          unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
          unit(1, 'null'),
          unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')
        )
        final <- gtable::gtable_add_rows(panelGrobs[[1]], space)
        final <- rbind(final, panelGrobs[[2]], size = 'first')
        final$widths <- widths
        final$heights[panel_rows(final)$t] <- unit(c(1, params$zoom.size), 'null')
        final <- gtable::gtable_add_grob(final, indicator_v, 3, 2, 5, 
                                         z = -Inf, name = "zoom-indicator")
      }
    }
    final
  }
)

breeds <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F)
pcs <-read.table("results/ancestry/PCA_analysis.eigenvec")
eigenvals<- read.table("results/ancestry/PCA_analysis.eigenval") ##for calculating PCs variance
plot( pcs[,3], pcs[,4], xlab = "PC1", ylab = "PC2" )

##variance explained by eac PC
PVE <- eigenvals[1] / sum(eigenvals$V1)
PC1<-PVE[1,] ##0.1432354
PC2<-PVE[2,] ##0.1112072

#vector<-c("Benson",1,"unknown")
#add column called "indv" to dataframe pcs
vec <- seq(1,nrow(pcs))
pcs_new <- pcs %>%
  mutate(indv = vec) %>% relocate(indv, .before = V1)

#pcs file with breed names
pcs_with_breed <- merge(pcs_new, breeds, by.x = 1, by.y = 1, all.x = T) %>% relocate(CLUSTER, .after = indv)
pcs_with_breed[167,2]<-"Unknown" #< add Benson's breed
#cannot bind by ID for Benson--> need to fix this
colnames(pcs_with_breed)[5:9] = c( "PC1","PC2","PC3","PC4","PC5" )
pairs( pcs_with_breed[, 5:9] )
pcs_with_breed[,2] = as.factor(pcs_with_breed[,2])
pcs_with_breed$colour = as.integer(pcs_with_breed[,2])
pairs(pcs_with_breed[,3:7], col=pcs_with_breed$colour, pch=20, lower.panel=NULL)
legend(0.1, 0.5, col=1:max(pcs_with_breed$colour), legend = levels(pcs_with_breed[,2]), pch=20, xpd=NA )

#plot of just PC1 and PC2
plot( pcs_with_breed[,5], pcs_with_breed[,6], col=pcs_with_breed$colour, xlab = "PC1 (14.32%)", ylab = "PC2 (11.12%)" )
text( pcs_with_breed[,5], pcs_with_breed[,6], labels = pcs_with_breed[,2], cex = 0.7, pos= 3) 
##try ggplot2
# counts
popCount <- length(unique(pcs_with_breed$CLUSTER))

# colours
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")
# shapes
shapes <- c(15:18, 3, 4, 8 ,10)

pdf(file='results/ancestry/PCA_with_Przewalski.pdf', width=12, height=7)
ggplot(data=pcs_with_breed, aes(x=PC1, y=PC2, 
                         color=CLUSTER, fill=CLUSTER)) + 
  geom_point(size = 2,  alpha = 0.85, stroke = 0.7, aes(shape=CLUSTER)) +
  annotate(
    geom = "curve", x = 0.0, y = 0.1, xend = -0.037920100, yend = 0.09533110,
    curvature = .3, arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(geom = "text", x = 0.0, y = 0.1, label = "Benson", hjust = "left")+
  xlab("PC1 (14.32%)") + ylab("PC2 (11.12%)") +
  scale_shape_manual(values = rep(shapes, len = popCount)) + 
  scale_colour_manual(values = rep(cbPalette, len = popCount))+
  scale_x_continuous(breaks=trans_breaks(identity, identity, n = 5))+
  theme_pubr()+
  labs(shape = "Breeds", color = "Breeds", fill = "Breeds") +
  theme(legend.position="bottom",
        legend.key.size = unit(3, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
dev.off()
  
highlight_df <- subset(pcs_with_breed, CLUSTER == "Unknown")
##now i need to filter the PCA eigenvector file to remove Prez.
pcs_no_Prze <- subset(pcs_with_breed, CLUSTER != "Przewalski" & CLUSTER != "Przewalski-hybrid")

pdf(file='results/ancestry/PCA_no_Przewalski.pdf', width=12, height=7)
ggplot(data=pcs_no_Prze, aes(x=PC1, y=PC2, 
                                                 color=CLUSTER, fill=CLUSTER)) + 
  geom_point(size = 2,  alpha = 0.85, stroke = 0.7, aes(shape=CLUSTER)) +
  annotate(
    geom = "curve", x = -0.02, y = 0.1, xend = -0.037920100, yend = 0.09533110,
    curvature = .3, arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(geom = "text", x = -0.02, y = 0.1, label = "Benson", hjust = "left")+ 
  xlab("PC1 (14.32%)") + ylab("PC2 (11.12%)") +
  scale_shape_manual(values = rep(shapes, len = popCount)) + 
  scale_colour_manual(values = rep(cbPalette, len = popCount))+
  scale_x_continuous(breaks=trans_breaks(identity, identity, n = 5))+
  theme_pubr()+
  labs(shape = "Breeds", color = "Breeds", fill = "Breeds") +
  theme(legend.position="bottom",
        legend.key.size = unit(3, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) +
  facet_zoom2(xlim=c(-0.06,0.00), ylim=c(0.04,0.12), zoom.size=2.5, horizontal = T) ##zoom in particular range in y axis
dev.off()

## same thing for the PCA without Prz
pcs_noPrz <-read.table("results/ancestry/PCA_no_Prz.eigenvec")
eigenvals_noPrz<- read.table("results/ancestry/PCA_no_Prz.eigenval") ##for calculating PCs variance
plot( pcs_noPrz[,3], pcs_noPrz[,4], xlab = "PC1", ylab = "PC2" )

##variance explained by eac PC
PVE_noPrz <- eigenvals_noPrz[1] / sum(eigenvals_noPrz$V1)
PC1_noPrz<-PVE_noPrz[1,] ##0.1307375
PC2_noPrz<-PVE_noPrz[2,] ##0.0983190

vec_noPrz <- seq(1,nrow(pcs_noPrz))
pcs_new_noPrz <- pcs_noPrz %>%
  mutate(indv = vec_noPrz) %>% relocate(indv, .before = V1)


breeds_no_prz<- breeds[-c(80,81,82,83,84,110,111),]
vec_noPrz_br <- seq(1,nrow(breeds_no_prz))
breeds_new_noPrz <- breeds_no_prz %>%
  mutate(ID = vec_noPrz_br)

#pcs file with breed names
pcs_with_breed_noPrz <- merge(pcs_new_noPrz, breeds_new_noPrz, by.x = 1, by.y = 1, all.x = T) %>% relocate(CLUSTER, .after = indv)
pcs_with_breed_noPrz[160,2]<-"Unknown" #< add Benson's breed

##try ggplot2
# counts
popCount_noPrz <- length(unique(pcs_with_breed_noPrz$CLUSTER))

# colours
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")
# shapes
shapes <- c(15:18, 3, 4, 8 ,10)

pdf(file='results/ancestry/PCA_noPrz.pdf', width=12, height=7)
ggplot(data=pcs_with_breed_noPrz, aes(x=PC1, y=PC2, 
                                color=CLUSTER, fill=CLUSTER)) + 
  geom_point(size = 2,  alpha = 0.85, stroke = 0.7, aes(shape=CLUSTER)) +
  annotate(
    geom = "curve", x = 0.05, y = 0.09, xend = 0.09884140, yend = 0.061974900,
    curvature = .3, arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(geom = "text", x = 0.04, y = 0.1, label = "Benson", hjust = "left")+
  xlab("PC1 (13.07%)") + ylab("PC2 (9.83%)") +
  scale_shape_manual(values = rep(shapes, len = popCount_noPrz)) + 
  scale_colour_manual(values = rep(cbPalette, len = popCount_noPrz))+
  scale_x_continuous(breaks=trans_breaks(identity, identity, n = 5))+
  theme_pubr()+
  labs(shape = "Breeds", color = "Breeds", fill = "Breeds") +
  theme(legend.position="bottom",
        legend.key.size = unit(3, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
dev.off()


pdf(file='results/ancestry/PCA_noPrz_zoom.pdf', width=12, height=7)
ggplot(data=pcs_with_breed_noPrz, aes(x=PC1, y=PC2, 
                                      color=CLUSTER, fill=CLUSTER)) + 
  geom_point(size = 2,  alpha = 0.85, stroke = 0.7, aes(shape=CLUSTER)) +
  annotate(
    geom = "curve", x = 0.07, y = 0.08, xend = 0.09884140, yend = 0.061974900,
    curvature = .3, arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(geom = "text", x = 0.07, y = 0.082, label = "Benson", hjust = "left")+
  xlab("PC1 (13.07%)") + ylab("PC2 (9.83%)") +
  scale_shape_manual(values = rep(shapes, len = popCount_noPrz)) + 
  scale_colour_manual(values = rep(cbPalette, len = popCount_noPrz))+
  scale_x_continuous(breaks=trans_breaks(identity, identity, n = 5))+
  theme_pubr()+
  labs(shape = "Breeds", color = "Breeds", fill = "Breeds") +
  theme(legend.position="bottom",
        legend.key.size = unit(3, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )+
  facet_zoom2(xlim=c(0.05,0.12), ylim=c(0.0,0.12), zoom.size=2.5, horizontal = T) ##zoom in particular range in y axis

dev.off()
