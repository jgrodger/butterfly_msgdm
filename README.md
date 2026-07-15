# Readme for Butterfly MSGDM project

## Note, updating in progress 24/01/2025

Template for file documentation

X)Filename
	•What this does: 
	•Inputs
	•Outputs
----------------------------------------------------------------------------------------------------------
Order in which to run scripts: 

The first two scripts must be run in order

1_Process_CSV_Data.Rmd
2_butterflies_combine_env_data.Rmd

Then the following can be run in any order

3_butterflies_maps.Rmd	
4_Occupancy_analyses.Rmd
5_run_zeta_msgdm.Rmd
	•Note, this file calls other files to run analyses and generate paramterised reports

The following scripts are run via 5_run_zeta_msgdm.Rmd

6_Zeta_declines_decays_analyses.Rmd
	This is run by 5_run_zeta_msgdm.Rmd
7_Zeta_msgdm_change_seed_analyses.Rmd

Scripts

R_Scripts
butterfly_richness_turnover_functions.R
	•What this does:
		•loads functions used by other scripts, especially functions that are used by more than one script

1) 1_Process_CSV_Data.Rmd
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


Note
	•Some plots may not print		
2) 2_butterflies_combine_env_data.Rmd
	•What this does:
		•Data synthesis and processing: 
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
			•Plots richness against environmental variables
			•Calculates correlation between environmental variables
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
		•Processed Data: lists including site by species matrix, site by xy, 
		site.by.env, sites, species dataframes for
			•all butterfly species "./Data/Processed_Data/Spatial/all_species.rds"
			•widercountryside species "./Data/Processed_Data/Spatial/wc.rds"
			•habitat specialists INCLUDING sites with no habitat specialists "./Data/Processed_Data/Spatial/hs.rds"
			•habitat specialists EXCLUDING sites with no habitat specialists. "./Data/Processed_Data/Spatial/hs_no_zeros.rds"
		•Slimmed down site by env to just map environmental data•"./Data/Processed_Data/Spatial/map_all_env_vars.rds"))
# corrections made 02/06/2026 and 06/06/2026. I had previously inadvertently edited in a way that prevented filtering

3) 3_butterflies_maps.Rmd		
	•What this does: generates maps of the distribution of sites in our data,
		including showing density of sites and mean species richness for each butterfly group and makes maps of environmental variables
	•Inputs
		•Basemap "./Data./Maps./gb_multipolygon_simplified.rds"
		•Data "./Data/Processed_Data/Spatial/all_species.rds"
	•Outputs
		• Hexmap of site density "Output/Spatial/Figures/Hex_dens_transects.png/"
		• Fig. 1 multipanel hexmaps of site density and species richness "Output/Spatial/Figures/Hex_dens_transects.png/"
		• Fig. S1.1 multipanel maps of environmental variables "Output/Spatial/Figures/all_env_vars.png/"
		
4) 4_Occupancy_analyses.Rmd
	•What this does: 
		•Analyses and plots species contributions to variance of raw zeta 
		diversity for different zeta orders, grouping by strategy (habitat specialist/habitat generalist).
	•Inputs
	Data: all species data "./Data/Processed_Data/Spatial/all_species.rds"
	•Outputs
		•Fig S1.1 Species occupancy table "Output/Spatial/zeta_declines_decays/all_species/butterfly_occupancy.csv"
		•Table S1.1 Rank occupancy plot grouped by (habitat specialist/habitat generalist). "Output/Spatial/zeta_declines_decays/all_species/rank_occupancy.png"
		•Fig S 4.1 Plot of contributions of different butterfly species to variance in zeta diversity for interpretation of msgdm results "Output/Spatial/zeta_declines_decays/all_species/contribution_zeta_variance.png"

 https://zsmith27.github.io/rmarkdown_crash-course/lesson-7-parameterized-reports.html		

# Needs updating with current script names (updated partially 15/05/2026)
5) 5_run_zeta_msgdm.Rmd
	•What this does: runs r markdown files and renders parameterized reports 
	for data analyses. 
	Analyses done for the first three or all four groups below
		•all butterfly species: all_species.
		•wider countryside species: wc.
		•habitat specialists INCLUDING sites with no habitat specialists: hs.
		•habitat specialists EXCLUDING sites with no habitat specialists: 
		hs_no_zeros.
	The following analyses are done NEEDS UPDATING
		•Zeta diversity decline and decay. "Scripts/Rmd_Scripts/Zeta_declines_decays_analyses.Rmd"
		•MSGDMs 
			•9_process_msgdm_results.RMD to get figures and tables
			"Scripts/Rmd_Scripts/Zeta_msgdm_analyses.Rmd"
			•note: sensitivity and robustness not yet done. They need to be 
			added from butterflies_msgdm_prelim_22_Nov_2023.Rmd
		•GAMs "Scripts/Rmd_Scripts/GAM_analyses.Rmd"
		•GAMMs "Scripts/Rmd_Scripts/GAMM_analyses.Rmd"
	•Inputs
		•Functions.
			•"./Scripts/R_scripts/get_k_fold_performance.R"
			•"./Scripts/R_scripts/Zeta.order.mc.sam.R"
			•"./Scripts/R_scripts/get_validation_different_seeds_and_samps.R"
	•Outputs
		•markdown reports


6) 6_Zeta_declines_decays_analyses.Rmd
	•What this does
		•Creates a parameterised report
		•Called by run_zeta_msgdm.Rmd
		
		•Zeta diversity: Calculates values.
		•Zeta decline: Compares fit of exponential versus power-law decline in 
		zeta diversity with zeta order, for different ranges of zeta 
		order.
		•Zeta decay: Assesses and plots zeta decay for different orders of zeta.
	•Inputs: datasets produced by butterflies_combine_env_data.Rmd
			•all butterfly species "./Data/Processed_Data/Spatial/all_species.rds"
			•widercountryside species "./Data/Processed_Data/Spatial/wc.rds"
			•habitat specialists EXCLUDING sites with no habitat specialists. "./Data/Processed_Data/Spatial/hs_no_zeros.rds"
	•Outputs: For each dataset 
		•Zeta values: values up to zeta = 250 in table "~zeta_decline_exact.csv")
		•Zeta decline, for orders 1 to 10, 1 to 20 and 1 to 50
			•Zeta decline plots
				•"~zeta_decline_to_10.png"
				•"~zeta_decline_to_20.png"
				•"~zeta_decline_to_250.png"
			•Tables comparing aic and r-squared of exponential verus power law 
			decline
				•"~zeta_declines_performance_10.csv"
				•"~zeta_declines_performance_20.csv"
				•"~zeta_declines_performance_250.csv"
			•Nestedness test
				•Boxplots of Simpson, Sorensen and Jaccard dissimilarity "~zeta_boxplot.png"
				•Table testing difference between Simpson and Sorensen dissimilarity  "~nestedness_anova_table.rds"
		•Zeta decays
			•Output of Zeta.ddecays "~zeta_distance_decay.rds"
			•Zeta distance decay plot "~zeta_distance_decays.png"
		Markdown: Zeta_declines_decays_analyses.html
			
7) 7_Zeta_msgdm_change_seed_analyses.Rmd 
	•What this does: 
		•Creates a parameterised report
		•Called by run_zeta_msgdm.Rmd
		•Seed can be set as a paramter
		
	•Analyses 		
		•Pearson and Spearman correlations of predictors
		•MSGDM analyses
		•Calculates deviance explained by MSGDM analyses from deviance and null 
		deviance (note that this should be treated with caution in a non-linear 
		model and other performance measures should routinely be used in
		addition or instead.
		•Ispline plots
		•Deviance explained pie plots
	•Inputs: datasets produced by 2_butterflies_combine_env_data.Rmd
			•all butterfly species "./Data/Processed_Data/Spatial/all_species.rds"
			•widercountryside species "./Data/Processed_Data/Spatial/wc.rds"
			•habitat specialists EXCLUDING sites with no habitat specialists. "./Data/Processed_Data/Spatial/hs_no_zeros.rds"
	•Outputs: For each dataset 
		•ispline plots "~ispline_plot_[dataset]_[normalize_msgdm]_zeta_order[i].png"
		•pie plots "~msgdm_var_pie_[dataset]_[normalize_msgdm]_zeta_order[i].png"
		•msgdm outputs in list "[dataset]_[normalize_msgdm]/zeta_msgdms.rds" 
		• NOW THIS SHOULD ALSO MAKE VARIANCE PIE PLOTS, CHECK WORKING (LAST CHUNK ADDED NOT YET RUN)

8) 8_MSGDM_multipanel_plot.RMD
	•What this does: Makes figure 3, a multipanel plot of MSGDM splines
	•Inputs: files from 7_Zeta_msgdm_change_seed_analyses.Rmd Output
		•./Output/Spatial/msgdms/all_species_Simpson_4//figs_tables//msgdm_single_splines.rds
		•./Output/Spatial/msgdms/wc_Simpson_4//figs_tables//msgdm_single_splines.rds
		•./Output/Spatial/msgdms/hs_no_zeros_Simpson_4//figs_tables//msgdm_single_splines.rds
	•Outputs
		•Fig. 3
			•msgdm_multiplot.png

 9) 9_process_msgdm_results.RMD
	•Run by: 5_run_zeta_msgdm.Rmd
	•What this does: makes forest plots for msgdms
	•Inputs: msgdm outputs  from  7_Zeta_msgdm_analyses.Rmd in list "[dataset]_[normalize_msgdm]/zeta_msgdms.rds" 
	•Outputs: 
		•Single spline msgdm plots "[dataset]_[normalize_msgdm]_zeta_order[i]/msgdm_single_splines.rds"
		•Multipanel plot of single ispline plots (currently not written to output)
		•ispline plots for each variable grouped by zeta order
		•		
		•Fig. S4 
			•Forest plots coloured by predictor, including associated p values "[dataset]_[normalize_msgdm]_zeta_order[i]/forest_plot.png"
		•Variance (Deviance) explained "[dataset]_[normalize_msgdm]_zeta_order[i]/deviance_table.csv"
		•Variance partitionining "[dataset]_[normalize_msgdm]_zeta_order[i]/variance_partitioning_table.csv"
# run from 	 5_run_zeta_msgdm.Rmd 01/06/2026	


		
11) 10_GAM_compare_family_analyses.Rmd
	•Run by: 5_run_zeta_msgdm.Rmd
	•What this does
		•GAM Analyses with different distribution assumptions
			•Residual Plots (not output elsewhere)
			•Basis dimension checks (not output elsewhere)
			•Smooth plots (not output elsewhere)
			•Moran's I test (not output elsewhere)
			•AIC 	
	•Inputs: 
		•Datasets produced by 2_butterflies_combine_env_data.Rmd
		•Basemap "./Data./Maps./gb_multipolygon_simplified.rds"
	•Outputs:
		•AIC 	"./Data/GAM/[dataset]/family_aic_table.csv" 
# run from 	 5_run_zeta_msgdm.Rmd 02/06/2026	


	
11) 11_GAM_comare_spline_analyses.Rmd
	•Run by: 5_run_zeta_msgdm.Rmd
	•What this does
			•GAM analyses, using family with lowest AIC and comparing spline 
			methods
				•Compares AIC.
				•Plots residuals
			•For the GAM with best spline method
				•Plots residuals on a map
				•Calculates concurvity of smoothed terms
				•approximate test null hypothesis that each smooth is actually 
				zero
	•Inputs: 
		•Datasets produced by 2_butterflies_combine_env_data.Rmd
		•Basemap "./Data./Maps./gb_multipolygon_simplified.rds"
	•Outputs: 
		•GAM_compare_spline_analyses.csv: Report to inspect and choose best distribution (family argument) to use.
		•spline_aic_table.csv: AIC table "/Output/Spatial/GAM/GAM_[dataset]/spline_aic_table.csv"
		• variograms
# run from 	 5_run_zeta_msgdm.Rmd 02/06/2026	
		
	
12) 12_GAMM_compare_corstructs_analyses
	•Standalone, run this script directly (NOT run by 5_run_zeta_msgdm.Rmd)
	•What this does
		•Finds best correlation structure (spatial autocorrelation) for each dataset
		•For each dataset
			•Prepares input for GAMMs
			•Runs GAMMs with different correlation structures
			•Writes results to file
	•Inputs: datasets produced by 2_butterflies_combine_env_data.Rmd
	•Outputs: 
		For each dataset, gamm input
			•"./Data/Processed_Data/Spatial/all_species_gamm_input.rds"
			•"./Data/Processed_Data/Spatial/hs_gamm_input.rds"
			•"./Data/Processed_Data/Spatial/wc_gamm_input.rds"
		For each dataset, list with all GAMM outputs 
			•"/Output/Spatial/GAM/GAMM_all_species/gamm_cor_struct_compare.rds"
			•"/Output/Spatial/GAM/GAMM_wc/gamm_cor_struct_compare.rds"
			•"/Output/Spatial/GAM/GAMM_hs/gamm_corExp.rds"
# Have not rerun GAMMS again in 2026, this requires significant computer time	



13) 13_process_GAMM_analyses
	•Run by: 5_run_zeta_msgdm.Rmd
	•What this does: 
		•Performs analyses and makes plots to choose best correlation structure
			•AIC table
			•k-index table and residual plots from gam.check
			•Custom residual plots because the ones from gam.check do not take into account family of distribution, so are not reliable
			•Moran's tests
			•Concurvity
			•GAM summaries
			•Smooth plots with partial residuals
		•Inputs
			•GAMM outputs from script 12_GAMM_compare_corstructs_analyses.Rmd
	•Outputs
		•Paramaterised reports for each dataset
			•Output/Spatial/GAM/GAMM_all_species/gamm_cor_struct_compare/13_process_GAMM_analyses.html
			Output/Spatial/GAM/GAMM_wc/gamm_cor_struct_compare/13_process_GAMM_analyses.html
			Output/Spatial/GAM/GAMM_hs/gamm_cor_struct_compare/13_process_GAMM_analyses.html


# Below still needs annotation and updating in this readme and the powerpoint	


14_GAMM_model_selection.Rmd
	•What this does: 
		•Runs GAMMs sequentially dropping ns variablles
	•Inputs
	•Outputs

15_process_GAMMS_single_result_output
	•Run by: 5_run_zeta_msgdm.Rmd
	•What this does: 
		•Drops ns terms sequentially from model
	•Inputs
		GAMM output from 14_GAMM_model_selection.Rmd
	•Outputs

GAMM_multipanel_plot.Rmd
	•What this does: 
	•Inputs
	•Outputs


	

GAMM_model_selection.Rmd
	•What this does: 
	•Inputs
	•Outputs
	
#below script no longer used, all this done elsewhere
GAMM_analyses.Rmd
	•Run by: 5_run_zeta_msgdm.Rmd
	•NB, these are very time consuming
	•What this does: Called by run_zeta_msgdm.Rmd;
	•analyses
		•GAMM analyses, using family and spline method from best GAM, compares
		GAMMs with different correlation structures for spatial autocorrelation
			•Compares AIC.
	•For the best GAMM
		•Plots residuals
		•Calculates concurvity of smoothed terms
	•Inputs: datasets produced by 2_butterflies_combine_env_data.Rmd
	•Outputs: For each dataset, list with all GAMM outputs 
	"~[dataset]_gamm_cor_struct_compare.rds"
	
11) butterfly_richness_turnover.Rmd
	•What this does: Generates a report containing figures and tables for main 
	text and supplementary materials
	•Inputs: 
		•Figures and tables generated by analyses files
		•Note, some chunks are manually pasted in from generate_code_chunks_for_supp_report.Rmd
	•Outputs: R markdown converted to html butterfly_msgdm_tables_figures_report.html

12)generate_code_chunks_for_supp_report.Rmd
	•What this does: Runs loops to automatically generated chunks for figures 
	to avoid errors
	•Inputs: NA.
	•Outputs: code chunks manually copied to butterfly_msgdm_tables_figures_report.Rmd

# Template for file documentation

X)Filename
	•What this does: 
	•Inputs
	•Outputs
