#data is  a dataframe
#site is a column of in the dataframe
#species is a column of species names in the dataframe

get.site.by.species <- function(data, sites, species, output = c("obs.count", "presence.absence")){
 
  output <- match.arg(output)
   
  sites <- enquo(sites)
  species <- enquo(species)
  
  #consider changing sep = ".", may cause issues, prefer sep = "_"
  
  # we do this to get a variable to count occurrences in each cell
  data <- data %>%  
    select(!!sites, !!species) %>%
    rename( sites = !!sites,
            species = !!species)%>%
    mutate(sites.species = as.factor(paste(sites, species, sep = ".")),
           sites = as.factor(sites),
           species = as.factor(species))
  
  
  levels.species <- levels(data$species)
  num.species <- length(levels.species)
  levels.sites <- levels(data$sites)
  num.sites <- length(levels.sites)
  
  
  levels.sites.species <- levels(data$sites.species)
  
  #here we get a dataframe with factorial combinations of species and sites, to get 0 occupancy
  species.vec <- rep(levels.species, num.sites)
  
  sites.vec <- rep(levels.sites, each = num.species)
  
  levels.frame <- data.frame(cbind(species.vec, sites.vec))
  
  levels.frame$site.species.factorial <- paste(levels.frame$sites.vec, levels.frame$species.vec, sep = ".")
  
  
  
  #here we add factor levels so that we can count zero occurrences later
  data$sites.species <- factor(data$sites.species, levels = levels.frame$site.species.factorial)
  
  #here we count occurrence of each species in each square, including zeros
  
  #this only works for presence data
 # species.sum <- data %>%
 #   count(sites.species, .drop = F)
  
  species.sum <- data %>%
    group_by(site.species) %>%
    summarise(sites.species, .drop = F
  
  #To understand why we need "\\." instead of merely "." we need to understand regular expressions in r
  #apparently "." has a meaning and we need to escape that by using "\\."
  species.sum <- species.sum %>% 
    mutate(sites.species = as.character(sites.species)) %>% 
    separate(col = sites.species, into = c("sites", "species"), sep = "\\.", remove = T, extra = "merge")
  
  
  #here we get number of presences for each species in each cell
  occurrence.matrix <- species.sum %>%
    pivot_wider(names_from = sites, values_from = n)%>%
    data.frame
  
  
  #make species column into row names
  rownames(occurrence.matrix) <- occurrence.matrix[,1]
  
  #remove species column
  
  occurrence.matrix$species <- NULL
  
  if(output == "presence.absence") {
  
  #convert to presence absence
  occurrence.matrix[occurrence.matrix[, ] > 0] <- 1
  

    
 
  }
  occurrence.matrix <- occurrence.matrix %>%
     t()%>%
       data.frame() 
  
  #   interum <- list(levels.species, num.species, levels.sites, num.sites, levels.sites.species)
  return(occurrence.matrix)
  
  
}

