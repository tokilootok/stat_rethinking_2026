library(igraph)
library(rethinking)


blank(bg="black")

g <- watts.strogatz.game(1, 100, 1, 0.5, loops = FALSE, multiple = FALSE)

plot(g, vertex.label= NA, edge.arrow.size=0.05, 
    xlab = "" , vertex.color="white" , 
    edge.color="white", bg="black" , vertex.frame.color="white" , edge.width=3 , vertex.size=5)

