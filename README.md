# Readme for Butterfly MSGDM project


Template for file documentation

X)Filename
	•What this does: 
	•Inputs
	•Outputs


1)Process_CSV_Data.Rmd
	•What this does: 
		•Data cleaning and subsetting(filtering) for UKBMS records, sites and 
		species data.
		•Adds Easting, Northing; & 1km grid-reference to data
		•Obtains site by species matrices.
		•Data exploration:
			•Site localities: plot on basemap.
			•Survey effort.
			•Survey coverage.
	•Inputs
		•UKBMS records, sites and species datasets from csv files, 
		including a site data csv supplied by Rob Cooke with additional info.
		•Basemaps processed in another project. •"./Data./Maps./gb_multipolygon_simplified.rds"
			•"./Data./Maps./gb_multipolygon_incl_ork_shet_simplified.rds"
	•Outputs
		•Processed Data: a list "./Data/Processed_Data/Spatial/spatial_ukbms_obs_sites_spp.rds" 
		containing: 
			•A site by species matrix (presence absence at the 
			moment).
			•Modified versions of the sites and species datasets.
		•Descriptive data: tables in prisma style, detailing filtering criteria 
		and data kept/omitted:
			•"Output/Spatial/filtering//filtering_summary.csv" with all steps
			and including records, sites and species kept at each step.
			•"Output/Spatial/filtering//sites_filtering_summary.csv" only with 
			sites kept after data filtered to 2015-2019 period, for inclusion in
			methods.
			
2) butterflies_combine_env_data.Rmd
	•What this does:
		•Data synthsis and processing: 
			•Environmental variables:
				•Extracts data for environmental variables for sites included 
				by the other criteria into site by environment dataset.
				•Log transforms certain environmental variables	.
				•Makes correlation tables.
			•Sites, species, and site by species data:
				•subsets to only sites with data for all environmental records.
				•Produces input for zeta diversity, MSGDM and GAM analyses 
				including site by species matrix, site by xy, site.by.env,
				sites, species dataframes and places them in a list.
				•subsets the list containing data for all species into: wider 
				countryside species; habitat specialists INCLUDING sites with 
				no habitat specialists; habitat specialists EXCLUDING sites 
				with no habitat specialists.
		•Data exploration: Environmental variables:
			•Maps
			•Histograms of environmental variables for all grid squares in UK 
			versus our sites to eyeball	bias in site selection wrt environmental 
			variables
	•Inputs
		•Data
			•Butterfly Data
				•List of datasets from Process_CSV_Data.Rmd. "./Data/Processed_Data/Spatial/spatial_ukbms_obs_sites_spp.rds"
		•Environmental data maps processed in other projects
			•human population density: "./Data/Masked_Cropped_Environmental_Data/site_by_humans.rds"
			•ph: "./Data/Masked_Cropped_Environmental_Data/site_by_ph.rds"
			•tree density: "./Data/Input_Environmental_Data/trees_1km_df.rds"
			•nitrogen deposition: "./Data/Masked_Cropped_Environmental_Data/site_by_nitrogen.rds"
			•pesticide.risk: "./Data/Masked_Cropped_Environmental_Data/site_by_pesticide_risk.rds"
			•topographic wetness: "./Data/Masked_Cropped_Environmental_Data/site_by_wetness.rds"
			•mean annual temperature: "./Data/Masked_Cropped_Environmental_Data/site_by_temp.rds"
			•mean annual precipitation: "./Data/Masked_Cropped_Environmental_Data/site_by_prec.rds"
		•A basemap processed in another project. "./Data./Maps./gb_multipolygon_simplified.rds"
	•Outputs
		• basic site locations map "./Output/Spatial/env_var_exploration/1122_sites_maps.pdf"
		•Histogram plots
			•"./Output/Spatial/env_var_exploration/env_var_histograms_1.pdf"
			•"./Output/Spatial/env_var_exploration/env_var_histograms_2.pdf"
			•"./Output/Spatial/env_var_exploration/env_var_histograms_3.pdf"
		•Processed Data: lists including site by species matrix, site by xy, site.by.env,
				sites, species dataframes for
			•all butterfly species "./Data/Processed_Data/Spatial/all_species.rds"
			•widercountryside species "./Data/Processed_Data/Spatial/wc.rds"
			•habitat specialists INCLUDING sites with no habitat specialists "./Data/Processed_Data/Spatial/hs.rds"
			•habitat specialists EXCLUDING sites with no habitat specialists. "./Data/Processed_Data/Spatial/hs_no_zeros.rds"

## Needs tidying and annotating, and eventually this can be pasted into butterflies_combine_env_data.Rmd
2a) butterfly_maps.Rmd		
	•What this does: generates maps of the distribution of sites in our data,
		including showing density of sites. Various attempts were made, and I am
		going with the hexplots
	•Inputs
		•Basemap "./Data./Maps./gb_multipolygon_simplified.rds"
		•Data "./Data/Processed_Data/Spatial/all_species.rds"
	•Outputs
		• hexmap of site density"Output/Spatial/Figures/Hex_dens_transects.png/"
		
		
3) run_zeta_msgdm.Rmd call Scripts to run analyses



combine environmental data and do zeta diversity analyses
see also https://zsmith27.github.io/rmarkdown_crash-course/lesson-7-parameterized-reports.html
Order to run:


