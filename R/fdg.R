# force directed graph on all pairs association matrix
fdg <- function(dataSet, dataName=NULL, method="A", cutoff=0.1, dim=2){
  
  d <- dataSet[complete.cases(dataSet),]
  d <- d[,sapply(d,is.numeric)]
  
  # set up adjacency matrix
  if (method == "A") {
    adj <- tap(d)
    adj <- adj + diag(nrow(adj))
    
  } else {
    adj <- cor(d)^2
    adj <- adj - diag(nrow(adj))
  }
  
  # apply cutoff
  adj[adj<cutoff] <- 0
  
  # generate graph
  gr <- graph.adjacency(as.matrix(adj), weighted=TRUE,  mode="upper")
  
  # color the largest cliques green
  V(gr)$color <- "LightBlue"
  V(gr)[unlist(largest.cliques(gr))]$color <- "LightGreen"
  
  if (dim == 2) {
    V(gr)$label <- V(gr)$name
    E(gr)$label <- round(100*E(gr)$weight,0)
  } else {
    V(gr)$label <- paste("   ",V(gr)$name,sep=" ")
    V(gr)$label.color <- "white"
  }
  
  # layout <- layout.fruchterman.reingold(gr,dim=dim,coolexp=1)
  layout <- layout.kamada.kawai(gr,dim=dim,coolexp=0.99)
  main <- "Force Directed Graph\n"
  if(!( is.null(dataName))) main <- paste(main, "name ~",dataName,",")
  main <- paste(main,"attraction ~",method,",")
  main <- paste(main,"cutoff ~",100*cutoff,"%")
  
  if(dim == 2) {
    plot(gr, layout=layout, main=main)
  } else {
    rglplot(gr, layout=layout)
  }
}
