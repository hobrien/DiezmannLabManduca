library(scales)
library(RColorBrewer)
library(grid)

sd_theme <- function() {
  theme_bw() +
    theme(plot.title=element_text(size=14)) +
    theme(legend.position="none") +
    theme(plot.margin = unit(c(1,1,.6,1), "cm")) +
    theme(strip.text.x = element_text(size = 10, face='bold'))  +
    theme(strip.text.y = element_text(size = 10, face='bold'))  +
    theme(panel.background=element_blank()) +
    theme(plot.background=element_blank()) +
    theme(legend.background = element_blank()) +
    theme(legend.key = element_blank()) 
}  

pdf_theme <- function() {
  theme_bw() +
    theme(plot.title=element_text(size=14)) +
    theme(legend.position="none") +
    theme(plot.margin = unit(c(1,1,.6,1), "cm")) +
    theme(strip.text.x = element_text(size = 10, face='bold'))  +
    theme(strip.text.y = element_text(size = 10, face='bold'))  +
    theme(panel.background=element_blank()) +
    theme(plot.background=element_blank()) +
    theme(legend.background = element_blank()) +
    theme(legend.key = element_blank()) 
}  

html_theme <- function() {
  theme_bw() +
    theme(plot.title=element_text(size=14)) +
    theme(legend.position="none") +
    theme(plot.margin = unit(c(1,1,.6,1), "cm")) +
    theme(strip.text.x = element_text(size = 10, face='bold'))  +
    theme(strip.text.y = element_text(size = 10, face='bold'))  +
    theme(panel.background=element_blank()) +
    theme(plot.background=element_blank()) +
    theme(legend.background = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(legend.text=element_text(size=20), legend.title=element_text(size=22)) +
    theme(axis.text.x = element_text(size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(size=20,angle=0,hjust=1,vjust=0,face="plain"),  
          axis.title.x = element_text(margin=margin(30,0,0,0), size=22,angle=0,hjust=.5,vjust=0,face="bold"),
          axis.title.y = element_text(margin=margin(0,30,0,0), size=22,angle=90,hjust=.5,vjust=.5,face="bold")) +
    theme(legend.position='top') 
    
}  


fte_theme <- function() {
    
    # Generate the colors for the chart procedurally with RColorBrewer
    palette <- brewer.pal("Greys", n=9)
    color.background = palette[2]
    color.grid.major = palette[3]
    color.axis.text = palette[6]
    color.axis.title = palette[7]
    color.title = palette[9]
    
    # Begin construction of chart
    theme_bw(base_size=12) +
        
        # Set the entire chart region to a light gray color
        theme(panel.background=element_rect(fill=color.background, color=color.grid.major)) +
        #theme(plot.background=element_rect(fill=color.background, color=color.grid.major)) +
        theme(panel.border=element_rect(color=color.grid.major)) +
        
        # Format the grid
        theme(panel.grid.major=element_line(color="white",size=.45)) +
        theme(panel.grid.minor=element_line(color="white",size=.25)) +
        theme(axis.ticks=element_blank()) +
        
        # Format the legend, but hide by default
        theme(legend.position="none") +
        #theme(legend.background = element_rect(fill=color.background)) +
        theme(legend.key = element_rect(fill=color.background, colour=color.background)) +
        theme(legend.text = element_text(size=9,color=color.axis.title)) +
        theme(legend.title = element_blank()) +
        
        # Set title and axis labels, and format these and tick marks
        theme(plot.title=element_text(color=color.title, size=14, vjust=1.25)) +
        theme(axis.text.x=element_text(size=12,color=color.axis.text)) +
        theme(axis.text.y=element_text(size=12,color=color.axis.text)) +
        theme(axis.title.x=element_text(size=14,color=color.axis.title, vjust=0)) +
        theme(axis.title.y=element_text(size=14,color=color.axis.title, vjust=1.25)) +
        
        # Plot margins
        theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}
reverselog_trans <- function(base = exp(1)) {
     trans <- function(x) -log(x, base)
     inv <- function(x) base^(-x)
     trans_new(paste0("reverselog-", format(base)), trans, inv, 
               log_breaks(base = base), 
               domain = c(1e-100, Inf))
}

