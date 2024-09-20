## =================================================================== ##
## Set layout for plotting
## =================================================================== ##
plot_layout <- function(nrow, ncol, heights = NULL, widths = NULL) {
  
  # Function to test if sth is a interger
  is.wholenumber <-
    function(x, tol = .Machine$double.eps ^ 0.5)  abs(x - round(x)) < tol
  
  # Check if nrow or ncol are "null unit"
  if (is.null(heights)) {
    heights = unit(rep_len(1, nrow), "null")
  } 
  if (is.null(widths)) {
    widths = unit(rep_len(1, ncol), "null")
  }
  
  # Set layput for later plotting
  grid.newpage()
  layout <- grid.layout(nrow = nrow, ncol = ncol, 
                        heights = heights,
                        widths = widths)
  vp_overall <- viewport(name = "vp_overall", width = unit(0.95, "npc"), height = unit(0.95, "npc"), layout = layout)
  pushViewport(vp_overall)
}


###===================================================================###
## function to draw intra chromosomal interactions
###===================================================================###
plot_intra_chr_interaction <- function(interaction, cyto, peaks, return = FALSE) {
  #interaction: the filtered cluster output
  #cyto: stnadard cytoband file
  #peaks peak bedgragh file with scores at the 4th column
  
  # Set color scheme for each "stain" in cytoband data
  colnames(cyto) <- c("chrom", "start", "end", "name", "stain")
  cyto_color <- data.frame(stain = c("gpos100", "gpos75", "gpos66", "gpos50", "gpos33", "gpos25", "gneg", "acen", "gvar", "stalk"),
                           col = c("#4D4D4D", "#969696", "#AAAAAA", "#C3C3C3", "#C2C2C2", "#E6E6E6", "white", "red", "#F3E6C4", "#5A81AA"))
  cyto <- left_join(cyto, cyto_color, by = "stain", keep =F)
  
  # Order chromosome from chr1 to chr12
  #idx <- grepl("[0-9]+", cyto$chrom)
  #chrom_no <- paste("Chr", unique(c(sort(as.numeric(gsub("Chr", "", cyto$chrom[idx]))), gsub("Chr", "", sort(cyto$chrom[!idx])))), sep = "")
  chrom_no <- unique(cyto$chrom)%>%sort() #This is only compatible for seqnames like "Chr01".

  
  # Whether return color columns for RCiscos plot
  # interaction <- interaction[, 1:7] 
  interaction$PlotColor <- apply(interaction, 1, function(x){x[7] <- as.numeric(x[7]); if(x[7] == 4){a = "#FCB519FF"}
  else if(x[7] == 5){a = "#ED6925FF"}
  else if(x[7] == 6 | x[7] == 7){a = "#BB3754FF"}
  else if(x[7] >=8 & x[7] <= 12){a = "#781C6DFF"}
  else if(x[7] >= 13 & x[7] <= 20){a = "#330A5FFF"}
  else(a = "#000004FF")})
  
  ###define the loop colors by loop span
  
  #interaction$PlotColor <- apply(interaction, 1, function(x)
  #{
  #span <- as.numeric(x[6])-as.numeric(x[2])
 # 
 # if(span < 10^4){a = alpha("#b5d336",0.5)}
 # else if(span < 10^5){a = alpha("#36b779",0.5)}
 # else if(span < 10^6){a = alpha("#31698f",0.5)}
 # else(a = alpha("#441854",0.5))
 # 
 # }
 # )
  
  
  if (return == TRUE) {
    #col_div <- split(chrom_no, ceiling(seq_along(chrom_no)/4))
    #         names(col_div) <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    #                             "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")[1:length(col_div)]
    #        names(col_div) <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
    #                            "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")[1:length(col_div)]
    #names(col_div) <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    #                     "#66A61E", "#E6AB02", "#A6761D", "#666666")[1:length(col_div)]
    ##freq <- unique(as.numeric(interaction[[7]]))[order(unique(as.numeric(interaction[[7]])))]
    ##color = colorRampPalette(c("#E92F4A", "#F499A6"))(length(freq)) 
    #names(col_1) <- color
    ##freq_color <- data.frame(PlotColor = color, freq = freq)
    
    #         interaction$PlotColor <- apply(interaction, 1, function(x){if((x[1] %in% col_div[[1]]) | (x[4] %in% col_div[[1]])) {a = names(col_div)[1]}
    #                                                                    else if((x[1] %in% col_div[[2]]) | (x[4] %in% col_div[[2]])) {a = names(col_div)[2]}
    #                                                                    else if((x[1] %in% col_div[[3]]) | (x[4] %in% col_div[[3]])) {a = names(col_div)[3]}
    #                                                                    else if((x[1] %in% col_div[[4]]) | (x[4] %in% col_div[[4]])) {a = names(col_div)[4]}
    #                                                                    else if((x[1] %in% col_div[[5]]) | (x[4] %in% col_div[[5]])) {a = names(col_div)[5]}
    #                                                                    else {a = names(col_div)[6]}
    #                                                                   }
    #                                        )
    
    
    #interaction$PlotColor <- apply(interaction, 1, function(x){a = col_1[which(names(col_1) == x[7])]})
    ## interaction <- merge(interaction, freq_color, by.x = "V7", by.y = "freq", all.x = T)[, c("V1", "V2", "V3", "V4", "V5", "V6", "PlotColor")]
    #interaction$PlotColor <- apply(interaction, 1, function(x){freq_color[which(freq_color$freq == as.character(x[7])), ][["color"]]})
    #interaction <- interaction[, c(1:6, 8)]
    return(interaction)
  } else {
    #  interaction <- interaction[, c(1:6, 8)]
    # Order cytoband data
    cyto$chrom_no_f <- factor(cyto$chrom, levels = chrom_no)
    cyto <- cyto[with(cyto, order(chrom_no_f, start)), ]
    # Seperate centromere data from cytoband data
    centro <- cyto[cyto$stain == "acen", ]
    centro <- centro[, c("chrom", "start", "end")]
    cyto <- cyto[cyto$stain != "acen", ]
    cyto <- cyto[, c("chrom", "start", "end", "col")]
    
    # Convert variables to correct type
    # interaction <- interaction[which(as.character(interaction[[1]]) == as.character(interaction[[4]])), ]
    cyto[, c(1, 4)] <- sapply(cyto[, c(1, 4)], as.character)
    cyto[, 2:3] <- sapply(cyto[, 2:3], as.numeric)
    centro[, 1] <- as.character(centro[, 1])
    centro[, 2:3] <- sapply(centro[, 2:3], as.numeric)
    
    peaks[, 1] <- as.character(peaks[, 1])
    peaks[, 2:4] <- sapply(peaks[, 2:4], as.numeric)
    
    # Calculate max coordinate for each chromosome and the max peak signal
    cyto_max <- max(cyto[["end"]])
    max_score <- max(peaks[[4]])
    
    # Draw ideogram, interaction curves and peaks signals for each chromosome
    for (i in seq.int(length(chrom_no))) { 
      # Get chromosome name
      #if (grepl("X", chrom_no[i])) {
      #  chrom <- paste("Chr", "X", sep = "")
      #} else if (grepl("Y", chrom_no[i])) {
      #  chrom <- paste("Chr", "Y", sep = "")
      #} else if (grepl("C", chrom_no[i])) {# update by sun
      #  chrom <- paste("chr", "C", sep = "")
      #} else {
      #  chrom <- paste("Chr", i, sep = "")
      #}
      
      chrom <- chrom_no[i]
      
      # Filter data for ideogram
      chrom_cyto <- cyto[which(cyto[[1]] == chrom), ]
      chrom_cyto <- chrom_cyto[order(chrom_cyto[[2]]), ]
      chrom_max <- max(chrom_cyto[["end"]])
      chrom_centro <- centro[which(centro[["chrom"]] == chrom), ]
      
      # Filter data for peaks
      chrom_peaks <- peaks[which(peaks[[1]] == chrom), ]
      chrom_peaks <- chrom_peaks[order(chrom_peaks[[2]]), ]
      #chrom_peaks[[4]] <- chrom_peaks[[4]] / max_score * 10
      #chrom_peaks[[4]] <- chrom_peaks[[4]] / max_score 
      chrom_peaks[[4]][chrom_peaks[[4]] > 1] <- 1
      # chrom_peaks[[4]][chrom_peaks[[4]] > quantile(chrom_peaks[[4]], 0.95)] <- max_score
      
      
      # Filter data for inetraction curves
      chrom_interaction <- interaction[which(interaction[[1]] == chrom), ]
      chrom_interaction <- chrom_interaction[, c(2, 3, 5, 6, 8)]
      #print(head(chrom_interaction))
      
      # Plot chromosome name
      vp_chr = viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 1)
      grid.text(chrom, x = unit(0.5, "npc"), y = unit(0, "npc"),
                just = c("right", "bottom"), gp = gpar(col = "#323745", fontsize = 15) , vp = vp_chr) 
      
      # Plot chromosome ideogram
      vp = viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 2, name = "vp", xscale = c(0, cyto_max))
      
      
      # Cytoband information
      #   chrom_cyto <- rectGrob(x = unit(chrom_cyto[["start"]], "native"), y = unit(0.07, "npc"), width = unit((chrom_cyto[["end"]] - chrom_cyto[["start"]]) * 0.95, "native"), height = unit(0.86, "npc"),
      #                          just = c("left", "bottom"), gp = gpar(fill = chrom_cyto[["col"]], col = NA), vp = vp)            
      #   grid.draw(chrom_cyto)
      #   
      #    
      #   if (nrow(chrom_centro) == 0) {
      #       chrom_rect <- roundrectGrob(x = unit(-0.0015, "npc"), y = unit(0, "npc"), width = unit(chrom_max, "native") + unit(0.0015, "npc"),
      #                                   r = unit(1.3, "mm"),
      #                                   just = c("left", "bottom"), gp = gpar(fill = NA, lex = 2.5, col = "#9F9F9F"), vp = vp)   
      #       grid.draw(chrom_rect)
      #   } else {
      #       # Left part
      #       chrom_left <- chrom_centro[1, ]
      #       chrom_left$start <- 0
      #       chrom_rect_left <- roundrectGrob(x = unit(-0.0015, "npc"), y = unit(0, "npc"), width = unit(chrom_left[["end"]] - chrom_left[["start"]], "native") + unit(0.0015, "npc"),
      #                                        r = unit(1.3, "mm"),
      #                                        just = c("left", "bottom"), gp = gpar(fill = NA, lex = 2.5, col = "#9F9F9F"), vp = vp)   
      #       grid.draw(chrom_rect_left)
      #       
      #       # Right part
      #       chrom_right <- chrom_centro[2, ]
      #       chrom_right$end <- chrom_max
      #       chrom_rect_right <- roundrectGrob(x = unit(chrom_right[["start"]], "native"), y = unit(0, "npc"), width = unit((chrom_right[["end"]] - chrom_right[["start"]]) * 1.05, "native"),
      #                                         r = unit(1.3, "mm"),
      #                                         just = c("left", "bottom"), gp = gpar(fill = NA, lex = 2.5, col = "#9F9F9F"), vp = vp)            
      #       grid.draw(chrom_rect_right)
      #   }           
      
      
      grid.rect(x = unit(chrom_cyto[["start"]], "native"), y = unit(1, "npc"), width = unit(chrom_cyto[["end"]] - chrom_cyto[["start"]], "native"), height = unit(1, "npc"),
                just = c("left", "top"), gp = gpar(fill = chrom_cyto[["col"]], col = NA, col = "#9F9F9F"), vp = vp)    
      
      
      if (nrow(chrom_centro) == 0) {
        grid.rect(x = unit(0, "npc"), y = unit(1, "npc"), width = unit(chrom_max, "native"), height = unit(1, "npc"),
                  just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)   
        
      } else {
        # Left part
        chrom_left <- chrom_centro[1, ]
        chrom_left$start <- 0
        grid.rect(x = unit(0, "npc"), y = unit(1, "npc"), width = unit(chrom_left[["end"]] - chrom_left[["start"]], "native"), height = unit(1, "npc"),
                  just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)   
        
        
        # Right part
        chrom_right <- chrom_centro[2, ]
        chrom_right$end <- chrom_max
        grid.rect(x = unit(chrom_right[["start"]], "native"), y = unit(1, "npc"), width = unit(chrom_right[["end"]] - chrom_right[["start"]], "native"), height = unit(1, "npc"),
                  just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)            
        
      }
      
      
      # Plot interaction curves
      #             draw_curves <- function(vector) {
      #                 x_start = (vector[1] + vector[2]) / 2
      #                 x_end = (vector[3] + vector[4]) / 2
      #                 angle <- seq(0, pi, length=50)
      #                 polygonGrob(x = unit(seq(x_start, x_end, length=50), "native"),
      #                             y = unit(-sin(angle) , "npc"),
      #                             gp = gpar(fill = "#E92F4A", col = NA, alpha = 0.3, lwd = 0.001)
      #                 )
      #             } 
      
      draw_Bezier_Curve <- function(vector) {
        x_start = (as.numeric(vector[1]) + as.numeric(vector[2])) / 2
        x_end = (as.numeric(vector[3]) + as.numeric(vector[4])) / 2
        x1 <- x_start; y1 <- -0.2;
        x4 <- x_end; y4 <- -0.2;
        if ((x4 - x1) <= 100000) {
          x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 50; y2 <- -3;
          x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 50; y3 <- -3;
        } else if((x4 - x1) <= 1000000) {
          x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 5; y2 <- -3;
          x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 5; y3 <- -3;  
        } else {
          x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 0.5; y2 <- -3;
          x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 0.5; y3 <- -3;
        }
        bezierGrob(x = unit(c(x1, x2, x3, x4), "native"), y = c(y1, y2, y3, y4),
                   gp = gpar(col = as.character(vector[5]), fill = as.character(vector[5]), alpha = 0.5, lwd = 0.5))
      }
      
      curve <- apply(chrom_interaction, 1, draw_Bezier_Curve)
      curves <- gList()
      for (k in seq.int(curve)) {
        curves[[k]] <- gList(curve[[k]])
      }
      vp_interactions <- viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 2, name = "vp_interactions", xscale = c(0, cyto_max), just = c("left", "top"))
      pushViewport(vp_interactions)
      grid.draw(curves)
      upViewport()
      
      #  Plot peaks
      vp_peaks <- viewport(layout.pos.row = 3 * i - 2, layout.pos.col = 2, name = "vp_peaks", xscale = c(0, cyto_max), yscale = c(0, max_score), just = c("left", "bottom"))
      if (nrow(chrom_peaks) == 0) {
        next
      } else {
        grid.rect(x = unit(chrom_peaks[[2]], "native"), y = unit(0.09, "npc"), width = unit(chrom_peaks[[3]] - chrom_peaks[[2]], "native"), height = unit(chrom_peaks[[4]], "npc"),
                  just = c("left", "bottom"), gp = gpar(fill = "#77C4D2", col = "#77C4D2"), vp = vp_peaks)
      }
      
      
    }
  }
}