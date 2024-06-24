#functions for butterfly msgdm project



# A function to add numbers of records, species and sites at each step to a table
# Used in Process_CSV_Data.Rmd and process_temporal_data.Rmd
add_prisma <- function(tibble = prisma, recs = records, step_name = "?"){
  records <- recs
  prisma_temp <- tibble(
    "Step" = step_name,
    "Records" = nrow(records),
    "Sites in records" = length(unique(records$site)),
    "Species in records" = length(unique(records$species.no)),
    "Sites in sites" = nrow(sites),
    "Species in species" = nrow(species)) 
  new_prisma <- rbind(prisma, prisma_temp)
}

# A function for filtering data: this assumes that the environment contains objects 
# called records, sites, species as prepared above. The argument is the dataset that 
# will already have been filtered by another criterion. The other 
# two are then filtered to contain records still relevant
filter_ukbms <- function(filter_by){
  
  if(filter_by == "records") {
    species <- filter (species, species.no %in% records$species.no)   
    sites <-  filter(sites, site %in% records$site)  
    
  } else if  ((filter_by == "sites")) {
    records <- filter(records, site %in% sites$site)
    species <- filter (species, species.no %in% records$species.no) 
    
  } else if (filter_by == "species"){
    records <- filter(records, species.no %in% species$species.no)
    sites <-  filter(sites, site %in% records$site)
  }
  
  filtered_data <- list(records, sites, species)
  names(filtered_data) <- c("records", "sites", "species")
  return(filtered_data)
}

# A function to filter the rest of the output dataset after the site by species data is filtered/subsetted
# Used in butterflies_combine_env_data.Rmd
# Used in Process_CSV_Data.Rmd and process_temporal_data.Rmd

filter_dataset <- function(site.by.species, site.by.xy, site.by.env, species, sites){
  site.by.xy <- filter(site.by.xy, rownames(site.by.xy) %in% rownames(site.by.species))
  site.by.env <- filter(site.by.env, rownames(site.by.env) %in% rownames(site.by.species))
  sites <- filter(sites, rownames(sites)  %in% rownames(site.by.species)) 
  species <- filter(species, species.no %in% names(site.by.species))
  
  dataset <- list(site.by.species, site.by.xy, site.by.env, species, sites)
  names(dataset) <- c("site.by.species", "site.by.xy", "site.by.env", "species", "sites")
  
  return(dataset)
}