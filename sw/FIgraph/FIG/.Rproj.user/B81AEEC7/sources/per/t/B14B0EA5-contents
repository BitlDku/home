#' Draw feature interaction graph
#'
#' This function Draw feature interaction graph.
#'
#' @import igraph
#' @param th Threshold for showing edge. If th is NULL, th is mean of interaction values of all pair of teatures
#' @param gtype Graph style. 'S' (default) is for static graph and 'I' is for interactive graph.
#' @param FIobj Result object from FItable function.
#' @param task Type of prediction task, 'regression' (default), 'Classification'
#' @param class.name It is for classification task. You can draw feature interaction graph for specific class.
#' @param show.edge.weight Decide if edge weight value is displayed or not.
#' @param seed According to seed value, The shape of graph is changed.
#' @return Feature interaction graph
#' @examples
#' library("FIG")
#' # for regression
#' data("Boston", package = "MASS")
#' model <- lm(medv ~ ., data = Boston)
#' FIobj1 <- FItable(model, train=Boston, target.name="medv", grid=50,
#'                  task="regression", interaction_type="OH2")
#' print(FIobj1)
#' FIgraph(gtype="S", FIobj=FIobj1, task="regression", seed=101)
#'
#' # for classification
#' library(e1071)
#' model2 <- svm(Species~., data=iris)
#' FIobj2 <- FItable(model2, train=iris, target.name="Species", grid=50,
#'                  task="classification", interaction_type="OH2", all.class=F)
#' print(FIobj2)
#'  FIgraph(th=NULL, gtype="I", FIobj=FIobj2, task="classification", class.name='setosa',
#'          show.edge.weight=T, seed=102)
#'
#' # Please ignore warning error when graw graph.
#'
#' @export
FIgraph <- function(th=NULL, gtype="S", FIobj, task="regression", class.name="all",
                    show.edge.weight=TRUE, seed=100) {
  # check error
  c.name <- unique(FIobj$Fint$class)
  if (length(c.name) > 1 & class.name == "all") {
    print("Error! You should give 'class.name'")
    return(NULL)
  }
  if (!(class.name %in% c.name) & !is.null(c.name)) {
    cat("Error! 'class.name' is one of [", c.name, "] \n")
    return(NULL)
  }

  # Get Feature interation  & feature importance
  result <- FIobj$Fint   # feature interaction table
  myimp  <- FIobj$Fimp   # feature importance table

  # filter low interaction
  if(is.null(th)) {
    th <- mean(result$weight)
    if (task=='classification') th <- mean(result$weight[result$class==class.name])
  }
  result <- result[result$weight>=th,]
  cat("th: ", th, "\n")

  ## Draw interaction Graph #############################################


  if (task=="regression") {
    edges <- result
    nodes <- data.frame(myimp)
  } else {
    if (is.null(class.name)) class.name = "all"
    edges <- result[result$class==class.name,-1]
    nodes <- data.frame(myimp[myimp$class==class.name,-1])
  }

  net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)

  # node size
  impsize <- (V(net)$importance - min(V(net)$importance)) /
    (max(V(net)$importance)-min(V(net)$importance))
  impsize <- ((impsize/1.2)+0.2) *50
  V(net)$size <- impsize # V(net)$importance*9  # node size

  # node color
  pal1 <- heat.colors(nrow(nodes), alpha=1)
  idx <- (V(net)$importance - min(V(net)$importance)) /
    (max(V(net)$importance)-min(V(net)$importance))
  idx <- nrow(nodes)-as.integer(idx*(nrow(nodes)-1)+1)+1
  V(net)$color <- pal1[idx]


  E(net)$label <- edges$weight     # edge weight
  #E(net)$label.cex=1
  #E(net)$label.font=2

  # edge width
  ewidth <- (E(net)$weight - min(E(net)$weight)) /
    (max(E(net)$weight)-min(E(net)$weight))
  ewidth <- (ewidth + 0.2)*7
  E(net)$width <- ewidth # edges$weight*20  # edge width

  # edge color
  E(net)$color <- "gray78"
  E(net)$color[E(net)$width < 0] <- "orange"

  # Draw graph

  title <- ""
  if (task=="classification") title <- paste0("class: ", class.name)

  if(!show.edge.weight) {
    elabel <- NA
  } else {
    elabel <- E(net)$value
  }

  set.seed(seed)
  if (gtype=="S") {  # static graph
    plot(net, edge.arrow.size=.4, vertex.label=V(net)$nodes, edge.label=elabel,
         main=title)
  } else {            # interactive graph
    id <- tkplot(net, edge.label=elabel)
    #tk_close(id, window.close = T)
  }

}

