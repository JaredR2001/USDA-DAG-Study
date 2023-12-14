library(lme4)
library(tidyverse)
library(bnlearn)
library(Rgraphviz)
library(tools)
library(furrr)
library(penalized)

# Variables to set are all here
###############################
useMarkers = TRUE
displaySaveData = FALSE

# Input data
InputTierDataFileName = "tiers_list1.csv"
InputTraitDataFileName <- "data_Flood.csv" 
InputMarkerDataFileName <- "data_markers.csv" 

# Output files
outTiffFile <- "m.Flood_wMarkers_ls1.tiff"
nodeListFile <- "Flooded_wMarkers_nodeList.csv"
droppedMarkerFile <- "droppedMarkers.csv"

markerCorrThresh <- 0.9999 #only perfect correlations removed from markers
alphaval <- 0.01
useRidge <- FALSE
kfold <- 10
kfoldruns <- 10
graphLayout <- "dot" #other options: neato, fdp, twopi, circo
graphMinThresh <- 0.5
#graphThresh <- 0.6
graphThresh <- 0.75
markerInteract <- FALSE #allow marker to marker edges
saveDAGfile <- "saveDAG.rds" #filename of saved run
relearnDAG <- TRUE #set to false to skip and use a saved run
###############################


no_cores <- availableCores() - 1
plan(multisession, workers = no_cores) #generic plan for use with furrr future_map

# future::plan( #custom plan for nested future_map
#   list(
#     future::tweak(
#       future::multisession, 
#       workers = 2), 
#     future::tweak(
#       future::multisession,
#       workers = 5)
#   )
# )

InputTraitDataFileName.md5 <- md5sum(InputTraitDataFileName)
print(InputTraitDataFileName.md5) #put in README metadata file
InputMarkerDataFileName.md5 <- md5sum(InputMarkerDataFileName)
print(InputMarkerDataFileName.md5) #put in README metadata file
InputTierDataFileName.md5 <- md5sum(InputTierDataFileName)
print(InputTierDataFileName.md5) #put in README metadata file

# Data preprocessing. Drop redundant or highly correlated traits. Set variable type (numeric, factor, ect.).
inputtraitdata <- read.csv(InputTraitDataFileName) %>%
  mutate_at("Genotype", str_replace, "TIL:", "") %>%
  select(-DaysToScdryT_100Is100, -DHD.j, -Node_BottomT, -Diff_LfAge_minus_N1T, -TN_sixwk, -RB_per_StemNu, -SB_per_StemNu, -RBbySB_ratio, -Node_UppermostT, -NeverGot_ScdryT, -Ht_Lftip_cm, -TotalRtLngth_cm, -ThckRtNu_sum, -Sum_ThickRoot_branches, -B_Proportion_LatB_over_total_Length, -TotRt_SA_cm2, -TotRt_Vol_cm3, -B_LatRtSA, -B_ThickRtSA, -B_ProportionSA_from_Lats, -B_LatRt_vol, -B_ThickRt_vol, -B_Proportion_Vol_is_Lat, -C_thinRtLngthSum_cm, -C_Proportion_Thin_over_total_Length, -C_thinRtSA, -GrowthRt, -TTill, -PctVegTill, -ISdWtPrPan, -SdDenPrPan, -SdCtPerPlt, -PltSdPerPan, -TYldPerPan, -TYldPerTill, -GrainLWR, -PanStress, -LeafStress) %>%
  rename(NodeFirstT = Node_FirstT) %>%
  rename(Ht_collar = Ht_collar_cm) %>%
  rename(LeafLength = LfLngth_TipHt_Minus_CollarHt) %>%
  rename(LatRtAvLngth = LatRtAvLngth_RtSystm_mm) %>%
  rename(LatRtDensity_Max = LatRtDensity_Max_Subsection) %>%
  rename(CrownRts = Nu_Post.embryonic_CrownRts) %>%
  rename(Upper_RtBranches = Nu_Seed.proximal_RtBranches) %>%
  rename(Mid_RtBranches = Nu_Mid.rt_Branches) %>%
  rename(Lower_RtBranches = Nu_TerminalRt_Branches) %>%
  rename(LatRtLngth_sum = B_LatLngthSum_cm) %>%
  rename(CoarseRtLngth_sum = B_ThickLngth_bySubtr_cm) %>%
  rename(AvgRtDiam_RtSystm = TotRt_AvgDiam_cm) %>%
  mutate(across(As:PctBrnChalk, as.numeric)) %>%
  mutate(across(As:PctBrnChalk, ~scale(., center = FALSE, scale = TRUE) %>% #scale to get Z-scores
                  as.numeric %>%
                  round(4))
  ) 

markerdata <- read.csv(InputMarkerDataFileName)

# Store the original marker dataset
original_markers <- names(markerdata)[-1]

#removes correlated markers 
markerdedup <- markerdata[,-1] %>%
  mutate_all(as.numeric) %>%
  dedup(threshold = markerCorrThresh, debug = TRUE)
markerdata <- cbind(Genotype=markerdata[,1],markerdedup)

#merge traits and markers
inputdata <- inputtraitdata %>% 
  inner_join(markerdata, by="Genotype") %>%
  select(-Genotype)

#get tier lists
tierdata <- read.csv(InputTierDataFileName)
tiercolors <- tierdata[1,] #get tier colors
tierdata <- tierdata[-1, , drop = FALSE] #remove row with colors

createTiers <- function(tierColNo, tdata) {
  tlist <- tdata[,tierColNo]
  tlist <- names(inputdata)[names(inputdata) %in% tlist] #eliminates traits not found in data file
  tlist <- tlist[nzchar(tlist)] #removes empty ""
  return(tlist)
}
traittiers <- 1:ncol(tierdata) %>% set_names(names(tierdata)) %>% 
  map(~ createTiers(.x,tierdata))

Markers <- names(markerdata[,-1])
traits <- unlist(traittiers) %>% unname() %>% unique() #get list of traits from blacklist

# Only use complete entries - use only complete cases
partial = inputdata[!complete.cases(inputdata), ]
inputdata = inputdata[complete.cases(inputdata), ]

fit.the.model = function(data, alpha) {
  
  if(useMarkers == TRUE) {
    learndata <- data[, c(Markers,traits)]
  } else {
    learndata <- data[, c(traits)]
  }
  
  learn.fun <- function(trait) {
    cat(paste0("fitting trait ",trait,"\n"))
    learntrait <-learn.nbr(learndata, node = trait, debug = FALSE, 
                           method = "si.hiton.pc", test = "cor", alpha = alpha)
    cat(paste0("fit trait ",trait,"\n"))
    return(learntrait)
  }
  
  cpc <- traits %>% set_names(traits) %>% future_map(~ learn.fun(.x) , .options = furrr_options(seed = TRUE))
  
  cat("Learn done\n")
  
  marker_Tier <- Markers[Markers %in% unlist(cpc)]
  
  nodes = unique(c(traits, unlist(cpc))) #keeping all traits but only markers that are parents of traits
  
  #make blacklist from tiers
  blacklisted = tiers2blacklist(traittiers)
  
  combined_blacklist <- tiers2blacklist(traittiers)
  
  if(useMarkers == TRUE) {
    blacklisted = tiers2blacklist(c(marker_Tier,traittiers))
    markerblacklist <- set2blacklist(marker_Tier)
    combined_blacklist <- rbind(markerblacklist, blacklisted)
  }
  
  cat("run hc\n")
  if (useMarkers == TRUE & markerInteract == TRUE) {
    bn = hc(data[, nodes], blacklist = blacklisted)  #allows marker to marker prediction
  } else {
    bn = hc(data[, nodes], blacklist = combined_blacklist) #does not allow markers to predict other markers
    #bn = tabu(data[, nodes], blacklist = combined_blacklist) #does not allow markers to predict other markers
  }
  return(bn)
}

xval.the.model = function(data, k, alpha, ridge) {
  
  n = nrow(data)
  data <- data[sample(n),] #shuffle rows 
  maxrows <- floor(n/k)*k # max rows to keep divisible by k
  data <- data[1:maxrows,]
  
  predcor = numeric(length(traits))
  names(predcor) = traits
  postcor = numeric(length(traits))
  names(postcor) = traits
  
  kcv = split(sample(maxrows), seq_len(k))
  
  predict.fun <- function(test) {
    pred = matrix(0, nrow = length(test), ncol = length(traits))
    colnames(pred) = traits
    
    post = matrix(0, nrow = length(test), ncol = length(traits))
    colnames(post) = traits
    
    cat("* beginning cross-validation fold.\n")
    
    dtraining = data[-test,]
    dtest = data[test,]
    
    cat("Fit the model\n")
    model <- fit.the.model(dtraining, alpha = alpha)
    cat("bn.fit\n")
    fitted = bn.fit(model, dtraining[, nodes(model)])
    
    cat("ridge\n")
    if (ridge) {
      
      for (no in nodes(fitted)) {
        node.parents = parents(fitted, no)
        if (length(node.parents) < 3)
          next
        
        opt.lambda = optL2(
          response = dtraining[, no],
          penalized = dtraining[, node.parents],
          model = "linear",
          trace = FALSE,
          minlambda2 = 10e-5,
          maxlambda = 500
        )$lambda
        fitted[[no]] = penalized(
          response = dtraining[, no],
          penalized = dtraining[, node.parents],
          model = "linear",
          trace = FALSE,
          lambda1 = 0,
          lambda2 = opt.lambda
        )
      }
    }
    
    dtest = dtest[, nodes(model)]
    cat("  > model has", length(nodes(model)), "nodes.\n")
    
    for (t in traits) {
      pred[, t] = predict(fitted, node = t, data = dtest[, nodes(model)])
    }
    if(useMarkers == TRUE) {
      for (i in seq(nrow(dtest))) {
        post[i, traits] = colMeans(cpdist(fitted, nodes = traits,
                                          evidence = as.list(dtest[i, names(dtest) %in% Markers]),
                                          method = "lw", n = 1000))
      }
    }
    
    return(list(model = fitted, pred = pred, post = post))
    
  }
  
  predicted <- kcv %>% future_map(~ predict.fun(.x), .options = furrr_options(seed = TRUE))
  
  posterior = do.call(rbind, lapply(predicted, `[[`, "post"))
  causal = do.call(rbind, lapply(predicted, `[[`, "pred"))
  
  cat("* overall cross-validated correlations:\n")
  for (t in traits) {
    predcor[t] = cor(causal[, t], data[unlist(kcv), t])
    cat("  > PREDCOR(", t, "):", predcor[t], "\n")
    postcor[t] = cor(posterior[, t], data[unlist(kcv), t])
    cat("  > POSTCOR(", t, "):", postcor[t], "\n")
  }
  
  
  return(list(predicted = causal, posterior = posterior,
              observed = data[unlist(kcv), t], predcor = predcor, postcor = postcor,
              models = lapply(predicted, `[[`, "model")))
  
}

if (relearnDAG == TRUE) {
  pr001 <- 1:kfoldruns %>% future_map(~ xval.the.model(inputdata, k = kfold, alpha = alphaval, ridge = useRidge), .options = furrr_options(seed = TRUE))
  saveRDS(pr001,saveDAGfile)
} else {
  pr001 <- readRDS(saveDAGfile)
}

pred.summary = sapply(pr001, `[[`, "predcor")
pred.means <- rowMeans(pred.summary)
print(pred.means)
pred.all.means <- mean(pred.means)

post.summary = sapply(pr001, `[[`, "postcor")
post.means <- rowMeans(post.summary)
print(post.means)
post.all.means <- mean(post.means)

cat(paste0("K-fold ", kfold, " predict ", round(pred.all.means,4), " posterior ", round(post.all.means,4), "\n"))

arclist = list()

for (i in seq_along(pr001)) {
  run = pr001[[i]]$models
  for (j in seq_along(run))
    arclist[[length(arclist) + 1]] = arcs(run[[j]])
}

nodes = unique(unlist(arclist))
strength = custom.strength(arclist, nodes = nodes)
averaged = averaged.network(strength, threshold = graphMinThresh)
relevant.nodes = nodes(averaged)[sapply(nodes, degree, object = averaged) > 0] #get nodes with at least one edge
averaged2 = bnlearn::subgraph(averaged, relevant.nodes)
strength2 = strength[(strength$from %in% relevant.nodes) &
                       (strength$to %in% relevant.nodes), ]
attr(strength2, "nodes") = relevant.nodes

gR = strength.plot(averaged2, strength2, shape = "ellipse", layout = graphLayout, threshold = graphThresh, groups=traittiers, render=FALSE)

tiff(outTiffFile, height = 35, width = 35, units = 'cm', compression = "lzw", res = 1080)

# Find correlations between variables
trait_correlations <- cor(inputdata[, traits], use = "pairwise.complete.obs")
marker_correlations <- cor(inputdata[, traits], inputdata[ , markerdata], method = "spearman")

#palette <- colorRampPalette(c("ivory", "ivory4"))(100)
negative_palette <- colorRampPalette(c("#FFBEB2", "#AE123A"))(100)
positive_palette <- colorRampPalette(c("#B9DDF1", "#2A5783"))(100)

# Set the number of colors based number of observed strengths
arc_colors <- rep(NA, nrow(strength2))
# Determine colors of arcs as it corresponds to the calculated correlations
for(i in 1:nrow(strength2)) {
  arc_from <- strength2$from[i]
  arc_to <- strength2$to[i]
  
  if (arc_from %in% traits && arc_to %in% traits){
    correlation_value <- trait_correlations[arc_from, arc_to]
    #normalized_value <- (correlation_value +1) / 2
    if(correlation_value > 0){
      arc_colors[i] <- positive_palette[ceiling(correlation_value * 100)]
    } else if (correlation_value < 0){
      arc_colors[i] <- negative_palette[ceiling(-correlation_value * 100)]
    } else {
      arc_colors[i] <- "white"
    }
    
    #arc_colors[i] <- palette[ceiling(normalized_value * 100)]
  } else {
    if (arc_from %in% Markers || arc_to %in% Markers){
      arc_colors[i] <- "darkgreen"
    }
  }
}

# Add the colors to the DAG
for(i in 1:nrow(strength2)) {
  arc_from <- strength2$from[i]
  arc_to <- strength2$to[i]
  arc_name <- as.character(interaction(arc_from, arc_to, sep = "~"))
  edgeRenderInfo(gR)$col[arc_name] <- arc_colors[i]
}

nodeRenderInfo(gR)$fill = "lightblue"
nodeRenderInfo(gR)$col = "darkblue"
nodeRenderInfo(gR)$fill[Markers] = "orange" # For markers 
nodeRenderInfo(gR)$col[traits] = "black"

for (i in 1:length(tiercolors)) {
  nodeRenderInfo(gR)$fill[traittiers[[i]]] = tiercolors[i]
}
# nodeRenderInfo(gR)$width <- 40
# nodeRenderInfo(gR)$height <- 40
# nodeRenderInfo(gR)$fontsize <- 25

a = arcs(bnlearn::subgraph(averaged2, relevant.nodes))
#a = as.character(interaction(a[, "from"], a[, "to"], sep = "~"))

#edgeRenderInfo(gR)$col = "grey"
#edgeRenderInfo(gR)$col[a] = "darkgreen"
renderGraph(gR)
dev.off()

# To see the black list that was applied
#write.csv(blacklisted, "C:/Users/jared.richardson/Desktop/USDASummer/R_Scripts/Bayesian Models/fourStudy_blacklisttable_Flood_09-12.csv")

# Store and save relevant nodes, including markers 
arcs <- averaged2[["arcs"]]
write.csv(arcs, nodeListFile, row.names = FALSE)

## Lets try to store dropped or unused markers (unused as a way to classify markers dropped AFTER deup)
if (useMarkers == TRUE) {
  # Place to store dropped markers
  unused_markers <- c()

  # Separate dataframe for "unused" markers
  unused_markers <- setdiff(original_markers, averaged2$arcs) # We'll just compare markers from the original set to the nodes present in the DAG 
  unused_markers_df <- data.frame(UnusedMarkers = unused_markers)

  write.csv(unused_markers_df, "unusedMarkers.csv", row.names = FALSE)
}

# To be added...First attempt at marginal graphs
# library(gRain)
# fitted_network <- pr001[[1]]$models[[1]]
# graphviz.chart(fitted_network, node = "LfNu")
# tiff(outTiffFile, height = 50, width = 50, units = 'cm', compression = 'lzw', res = 1080)

# To view the save file. Will be a .rds file ext.
if (displaySaveData == TRUE) {
  filename <- file.choose()
  saveData <- readRDS(filename)
  saveData
}

