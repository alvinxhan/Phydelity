# make data tables of phylogenies from simData
prepare_data_tables = function(vars,names) {
  library(data.table)
  library(phangorn)
  library(igraph)
  res = list()
  for(v in 1:length(vars)) {
    t = vars[[v]]
    n = length(t)
    for(i in 1:n) {
      n2=length(t[[i]])
      benchdata = data.table(id=1:n2,tree=list(),stepsize=rep(0.1,n2))
      for(j in 1:n2) {
        tmp = phyloObjFromSpreadMat(t[[i]][[j]]$transmissionFrame, dropUndiag = TRUE)
        tmp = reorder.phylo(tmp)
        x = floor(log10(min(tmp$edge.length)))
        benchdata$tree[[j]] = tmp
        benchdata$stepsize[j] = 10^x
      }
      
      setkey(benchdata,"id")
      tfile=paste0(names[v],'_IS_',i,'_data_table.RData')
      save(benchdata, file=tfile)
    }   
  }
}

####################################################################################################################
# Functions provided by Luc Villandre

phyloObjFromSpreadMat <- function(spreadMat, dropUndiag = FALSE) {
 
  simPhylo <- createPhyloTree(spreadMat) ## simPhylo has 4 elements, 
  #"infections" is the infection frame used to build the phylogeny, it's spreadMatScaled without the uninfected nodes, basically. 
  #"reverse" is an environment allowing us to traverse the tree either prospectively or retrospectively, I don't remember, 
  #"phylo" lists all parent nodes and their children nodes, parent uses the phylogeny label, infection should tell me which node in the population was infected, thus creating the split. 
  #Finally "branchSizes" tells you the supporting branchSizes for all internal nodes. External branch sizes are set later.
  phyloObj <- convertIntoPhylo(simPhylo)
  if (dropUndiag) {
    undiagNodes <- as.character(subset(spreadMat, subset = !is.na(infectionTime) & is.na(diagTime))$label)
    phyloObj <- drop.tip(phyloObj, tip = undiagNodes)
  } else{}
  phyloObj
}

.infectionIndex <- function(infectionTable) {
  #infectOnly <- subset(infectionTable, subset = infected == 1)
  infectOnly <- subset(infectionTable, subset = !is.na(infectionTime)) ## The root, whether it is inside or outside the sample, is infected at time 0.
  infectOnly <- infectOnly[order(infectOnly$infectionTime),]
  
  # Hash table to lookup new index given original
  labels <- infectOnly$label
  lookup <- new.env()
  lapply(seq_along(labels), FUN = function(i) {
    lookup[[as.character(labels[[i]])]] <- i
  })
  lookup[["0"]] <- 0 # Root is its own parent
  
  nodeList <- mget(x = as.character(infectOnly$parentNode), envir = lookup)
  infectOnly$parentNode <- do.call("c", nodeList)
  list(infections = infectOnly, reverse = lookup)
}

.createLeafNodes <- function(infections) {
  size <- 2 * nrow(infections) - 1
  data.frame(parent = rep(NA, size), c1 = NA, c2 = NA, infection = NA)
}

.generatePhylogeny <- function(phylo, infections) {
  nleaf <- nrow(infections)
  
  # Map of infection participants to phylo IDs
  map <- 1:nleaf
  
  ##nextRow <- nleaf + 1 # First non-leaf
  nextRow <- nleaf*2 - 1
  lapply(nleaf:2, FUN = function(tgt) {
    src <- infections$parentNode[tgt]
    mtgt <- map[tgt]
    msrc <- map[src]
    
    # Create a new node, update parent links
    phylo[nextRow,] <- c(parent = NA, c1 = mtgt, c2 = msrc,
                         infection = infections$label[tgt])
    phylo$parent[c(mtgt, msrc)] <- nextRow
    phylo <<- phylo
    
    # Update map to point src to parent
    map[src] <- nextRow
    map <<- map
    ## nextRow <<- nextRow + 1
    nextRow <<- nextRow - 1
  })
  
  phylo
}

.supBranchSizes <- function(leaves, phylo, spreadMat) {
  children <- (leaves + 1):(nrow(phylo) - 1)
  parents <- phylo$parent[children]
  mat <- matrix(spreadMat$infectionTime[phylo$infection[c(children, parents)]],
                ncol = 2)
  data.frame(phylo = children, size = mat[,1] - mat[,2])
}

createPhyloTree <- function(spreadMat) {
  index <- .infectionIndex(spreadMat)
  phylo <- .createLeafNodes(index$infections)
  phylo <- .generatePhylogeny(phylo, index$infections)
  branchSizes <- .supBranchSizes(nrow(index$infections), phylo, spreadMat)
  
  index$phylo <- phylo
  index$branchSizes <- branchSizes
  index
}

## The next fct takes the phylogeny outputted by createPhyloTree and turns it into a phylo type object, which will allow it to be used in phangorn
## This function should be modified to scrap the tip corresponding to the artificial root, if it exists.

convertIntoPhylo <- function(createPhyloOut, tipSupLength, defaultLength = 0.01) {
  excludeRoot <- subset(createPhyloOut$phylo, !is.na(parent)) ## the root is excluded
  edge <- cbind(excludeRoot$parent, as.numeric(rownames(excludeRoot)))
  
  environment(.getEdgeLength) <- environment()
  
  edge.length <- .getEdgeLength()
 
  ##tip.label <- paste("s", seq(nrow(createPhyloOut$infections)), sep ="")
  Nnode <- nrow(createPhyloOut$infections) - 1 ## The tree is rooted.
  tipLabelTemp <- paste("s", seq(Nnode+1), sep = "")
  
  phyloObj <- list(edge = edge, edge.length = edge.length, tip.label = tipLabelTemp, Nnode = Nnode)
  class(phyloObj) <- "phylo"
  
  tipLabels <- as.character(tipNodesId(createPhyloOut, phyloObj))
  
  phyloObj$tip.label <- tipLabels
  ##phyloObj <- .fixPhylo(phyloObj)
  ### We could remove branches for undiagnosed nodes right away with drop.tip, but should we? It would be more computationally efficient... We won't do it for now: We'll still simulate the sequences and drop them later, once their DNA configuration is simulated.
  ## This way, we can compare the results of clustering when all sequences are know with those when parts of the network are invisible.
  return(phyloObj)
}

.getEdgeLength <- function() {
  phyloExcludeTips <- subset(createPhyloOut$phylo, !is.na(c1)) ## The tips are excluded, as they do not have children.
  vectorAssign <- Vectorize(assign, vectorize.args = c("x", "value"))
  vectorGet <- Vectorize(get, vectorize.args = "x")
    
  vectorAssign(x = as.character(c(phyloExcludeTips$c1, phyloExcludeTips$c2)), value = rep(defaultLength, nrow(phyloExcludeTips)*2), envir = environment()) ## A bunch of objects named after the children labels are created and given value defaultLength
  vectorAssign(x = as.character(createPhyloOut$branchSizes$phylo), value = createPhyloOut$branchSizes$size, envir = environment()) ## Are all the internal branch lengths updated, leaving all the external branch lengths at defaultLength? Yes.
  ## This next bit is tricky. Containers for branch lengths have been created and initialized at defaultLength. Now, all containers for leaves (nodes 1 to num. of infections), must be given the right length corresponding to diagnostic time minus the last time they transmitted HIV.
  ## The "reverse" environment in phyloExcludeTips associates labels of nodes in the population (they're the object names in the environment) with leaf labels in the phylogeny (I just checked it, it works, as its elements match those in the label column of createPhyloOut$infections).
  ## So, for each leaf, the algorithm must check in createPhyloOut$infections when the node corresponding to it, that we find in the reverse environment, was diagnosed. It must then find the last transmission attributable to it. If it never transmitted HIV, we must find the infection time instead.
  diagTimes <- createPhyloOut$infections$diagTime
  diagTimes[is.na(diagTimes)] <- max(diagTimes, na.rm = TRUE) ## A trick to make the code lighter. All undiagnosed nodes are assumed diagnosed when the last diagnosis occurs, which will produce correct branch lengths, as the time of the last diagnosis is the clock time when the simulations stopped.
  currentEnvir <- environment()
  ## Now we find the latest infection times for each parent node.
  by(createPhyloOut$infections, INDICES = createPhyloOut$infections$parentNode, FUN = function(x){
    parentNode <- createPhyloOut$infections$label[x$parentNode[1]] ## [1] simply because x$parentNode is all identical values
    if (length(parentNode) == 0) return(NULL) ## If this is for the root of the epidemic (parent is labelled 0), do nothing.                                                                                    
    diagTime <- diagTimes[createPhyloOut$infections$label == parentNode] ## This should have length 1, as labels in the label column are unique.
    lastInfectionTime <- max(x$infectionTime)
    branchLength <- diagTime - lastInfectionTime ## This should always be greater than 0, as nodes should shut down once they are diagnosed.
    parentLeafLabInPhylo <- get(as.character(parentNode), envir = createPhyloOut$reverse)
    assign(as.character(parentLeafLabInPhylo), value = branchLength, pos = currentEnvir)    
  })
      
   ## The branch lengths here correspond to diagnostic time - last infection time for diagnosed nodes or clock time - last infection time for undiagnosed nodes.
  
  ## Infected nodes not listed in lastInfectionTimeList did not ever transmit. Branch lengths for them should be diagnostic time - infection time. This is what the next chunk of the code sets.
  
  parentPopLabels <- createPhyloOut$infections$label[createPhyloOut$infections$parentNode]
  nodesNoTransmit <- setdiff(createPhyloOut$infections$label, parentPopLabels) ## This is potentially empty, which becomes integer(0). It's empty when all nodes are responsible for at least one additional infection.
  if (length(nodesNoTransmit) > 0) { ## nodesNoTransmit is not empty
    phyloLeafIndex <- mget(as.character(nodesNoTransmit), envir = createPhyloOut$reverse)    
    branchLengthsNonTransmit <- (diagTimes - createPhyloOut$infections$infectionTime)[match(nodesNoTransmit, createPhyloOut$infections$label)] ## This vector is diagnostic time - infection time for diagnosed infections and clock time - infection time for undiagnosed infections. When time to diagnostic is fixed, some values are naturally repeated.    
    vectorAssign(x = as.character(do.call("c", phyloLeafIndex)), value = branchLengthsNonTransmit, envir = environment()) ## This resets the values of elements representing branch lengths for leaf nodes for non-transmitters in the environment.    
  } else{} 
  
  ##branchSizes is a matrix with two columns, one called phylo, the other called size. phylo is a child index. So, size is a supporting branch length?
  edge.length <- vectorGet(x = as.character(edge[,2]), envir = environment())
  names(edge.length) <- NULL ## edge.length gives supporting branch lengths.
  return(edge.length) ## I checked in one iteration: there are no negative external branch lengths and all of them are inferior or equal to the fixed diagnostic time I set (0.6 * 1/100, 1/100 being a time scaling factor to make clustering easier)
  ## The maximum of edge.length is 0.01, the initial value. Does this make sense?
  ## It does if it's for the root. Is the last labelled node in the phylogeny the root?
  ## It would make sense, given that leaves are numbered 1:number of infections, and that phylo labels are updated everytime two branches merge.
  ## Look at .fixPhylo: The final phylo node in edge.length is indeed the root. .fixPhylo would relabel it number of infections + 1.
}

## The root must be numbered (# of tips + 1)
## This function is actually never used. Why?
.fixPhylo <- function(phyloObj) {
  numTips <- length(phyloObj$tip.label)
  rootNum <- max(phyloObj$edge[,1])
  rootPos <- which(phyloObj$edge[,1] == rootNum)
  replaceFirstCol <- which(phyloObj$edge[,1] == (numTips + 1))
  ## Second column is always ordered, index is numTips + 1
  phyloObj$edge[,1] <- replace(phyloObj$edge[,1], rootPos, numTips + 1)
  phyloObj$edge[,1] <- replace(phyloObj$edge[,1], replaceFirstCol, rootNum)
  phyloObj$edge[numTips+1,2] <- rootNum
  
  return(phyloObj)
}

tipNodesId <- function(createPhyloOut, phyloObj) {
  
  .getNodePair <- function(node) {
    parent <- Ancestors(phyloObj, node = node, type = "parent")
    firstNode <- createPhyloOut$phylo$infection[parent]
    parentNodeIndex <- createPhyloOut$reverse[[as.character(firstNode)]]
    infectionsIndex <- createPhyloOut$infections$parentNode[parentNodeIndex]
    secondNode <- createPhyloOut$infections$label[infectionsIndex]
    contactNodes <- c(firstNode,secondNode)
    contactNodes
  }
  nleaf <- phyloObj$Nnode+1
  tipNodeVec <- integer(nleaf)
  lapply(seq(nleaf), FUN = function(x) {
     
    if (!(tipNodeVec[[x]] == 0)) { ## The leaf is already assigned a node.
      return(NULL)
    } else{}
    
    leafParentConfig <- .getNodePair(x) ## This gives the two nodes involved in the branching event leading to leaf x. One is the newly infected node, the other is the infector.
    ## Now we check if this ancestor has internal nodes as descendants.
    ## 0 should never be in there
    if (0 %in% leafParentConfig) stop("0 should not be in leafParentConfig")
    parentDescent <- Children(x = phyloObj, node = Ancestors(phyloObj, node = x, type = "parent"))
    if (!any(parentDescent > nleaf)) {
      tipNodeVec[parentDescent] <<- leafParentConfig
    } else {
      internalNodeIndex <- which(parentDescent > nleaf)
      ## We go down one node from the child internal node
      oneChildInternalNode <- Children(x = phyloObj, node = parentDescent[internalNodeIndex])[1]
      internalNodeConfig <- .getNodePair(oneChildInternalNode)
      
      leafNodeIndex <- which(!(leafParentConfig %in% internalNodeConfig))
      
      tipNodeVec[x] <<- leafParentConfig[leafNodeIndex]      
    }
    return(NULL)
  })
  return(tipNodeVec)
}