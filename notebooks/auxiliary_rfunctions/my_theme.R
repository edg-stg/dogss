theme_Ed <- function(base_size=14) { #
  (theme_foundation(base_size=base_size) #
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, vjust = 2),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.margin = margin(0, unit="cm"),
            legend.title = element_text(face="italic"),
            legend.text = element_text(size=rel(1)),
            plot.margin=unit(c(5,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

Ed_palette <- c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#984ea3","#ffff33","#fb9a99", "#484848")

# Ed_palette <- c("#003366", "#99cc00", "#007777", "#cc0022", "#ee9900", "#aaccee")

scale_fill_Ed <- function(...){
  # library(scales)
  discrete_scale("fill","Ed",scales::manual_pal(values = Ed_palette), ...)

}

scale_colour_Ed <- function(...){
  # library(scales)
  discrete_scale("colour","Ed",scales::manual_pal(values = Ed_palette), ...)

}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right"), which.legend=1) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[which.legend]] +
                    theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x +
                 theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,ncol = 1,
                                            heights = unit.c(unit(0.98, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)
}
