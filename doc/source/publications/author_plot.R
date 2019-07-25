# Code to make a network plot of ASPECT authors
# Jane Carlen, February 2019
# Modify paths and files for inclusion in ASPECT repo and website.
# LJ Hwang May 2019

# Modified for Rayleigh
# LJ Hwang July 2019

# ------------------------------------------------------------------------------------------
# 0. Setup  ----

library(readxl)
library(bibtex)
library(stringr)
library(igraph)
library(ggraph)
library(dplyr)
library(grDevices)
library(grid)

# This is now set to use the file all.bib.
# Files should be downloaded from the rayleigh repository /doc/source/publications
# CONCATENATE all the publicationYYYY.bib files
#       cat *.bib >>all.bib
# Note the files have headers but it seems read.bib just ignores these
#
# Adjustments maybe needed for inconsistent author naming
#
# SET YOUR DIRECTORY
#setwd("~/Documents/CIG/2019_Citation/Rayleigh")
setwd("~/Repos/Rayleigh/doc/source/publications")
#UPDATE TO CURRENT DATE
my_grob = grobTree(textGrob("Last Updated 25 July 2019", x=0.05,  y=0.05, hjust=0,
                            gp=gpar(col="black", fontsize=8, fontface="italic")))
#
# 1.0 Clean author names  ----

# read in bib file from working directory
# we are keeping the aspect variable name out of laziness
aspect.bib = read.bib("all.bib", encoding = "latin1") #had to add school placeholder 'test' for some articles so they'd load
aspect.authors.list = lapply(aspect.bib, function(x) as.character(x$author)); length(unique(unlist(aspect.authors.list))) 
#
# two authors get grouped when we use last name only
# From checking excelt, S. Zhang != N. Zhang,  but ch C. O\\textquoterightNeill == C. J. O\\textquoterightNeill ?
#
# Check for inconsistent authors names and add to list
aspect.authors.list = sapply(aspect.authors.list, function(x) {

  x[str_detect(pattern = "Rene Gassm\303\266ller", x)] = "R Gassm\303\266ller"; x
  x[str_detect(pattern = "Wolfgang Bangerth", x)] = "W Bangerth"; x 
  x[str_detect(pattern = "N. Featherstone", x)] = "N.A. Featherstone"; x 
  x[str_detect(pattern = "N. A. Featherstone", x)] = "N.A. Featherstone"; x 
  x[str_detect(pattern = "B. W. Hindman", x)] = "B.W. Hindman"; x 
  x[str_detect(pattern = "B. Buffett", x)] = "B.A. Buffett"; x 
  x[str_detect(pattern = "M. Miesch", x)] = "M. S. Miesch"; x 
  
}) 
sort(table(unlist(aspect.authors.list)))

# 2.0  Convert to co-author network  ----
author.pairs = lapply(aspect.authors.list, function(x) { if(length(x) > 1) {t(combn(x, 2))} else {NULL}})
author.pairs = data.frame(do.call("rbind", author.pairs))
author.pairs$count = 1
author.pairs = aggregate(author.pairs$count, list(author.pairs$X1, author.pairs$X2), sum)
names(author.pairs)[3] = "co_authorships"
# write out file in case you need to debug
# write.csv(author.pairs, "author_pairs.csv", row.names = F)
aspect.net = graph_from_data_frame(author.pairs, directed = F, vertices = as.vector(unique(unlist(aspect.authors.list))))
plot(aspect.net)

# 3.0  Make plot and add colors ----

# Special authors
V(aspect.net)$name_last = sapply(str_split(V(aspect.net)$name, " "), last)
red = c("Featherstone") #Original Lead Developers
#Added Principal Developers. Keeping this in for now but we are not using here 
#brown = c("Glerum", "Fraters", "Austermann", "Naliboff") 
brown = c(" ") 
color.aspect = c("Other", "Principal Developer")[1 + rowSums(sapply(red, grepl, vertex_attr(aspect.net, "name")))]
color.aspect[V(aspect.net)$name_last%in% brown] = "Principal Developer"
#red.rgb = "#990033"
brown.rgb = "#ff5050"
V(aspect.net)$Author_type =  color.aspect
V(aspect.net)$Papers = table(unlist(aspect.authors.list))[V(aspect.net)$name]
E(aspect.net)$co_authorships

author.plot = ggraph(aspect.net, layout = "kk") + 
  annotation_custom(my_grob) +
  geom_edge_link(aes(width = co_authorships), alpha = 1, color = "gray") +
  # for some reason edge width not working like it doesn in the plots below
  scale_edge_width(name = "Number of co-authorships", range=c(.5,3), breaks = c(1,2,3)) +
  geom_node_point(aes(fill = Author_type, size = Papers), shape = 21, stroke = .5, alpha = .8) + 
  scale_fill_manual(values =  c("white", brown.rgb),
                    guide = guide_legend(title = "Author type", override.aes = list(size = 5), order = 1)) + 
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        title = element_text(size = 16)
        ) +
  ggtitle("Rayleigh Co-author Relationships")
author.plot

ggsave(plot = author.plot, "author_plot_no_labels.png", height = 6, width  = 10)

# Add labels
author.plot = author.plot + geom_node_text(aes(label = name), repel = TRUE, segment.alpha = .5, fontface = "bold")
ggsave(plot = author.plot, "author_plot.png", height = 6, width  = 10)



