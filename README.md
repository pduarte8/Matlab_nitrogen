# Project Title: Dominance of physical mixing over biogeochemical processes in the Central Arctic mesopelagic
The scrips and functions listed and detailed below were used to analyse data from the Arctic Ocean
over 70 degrees northern latitude, imported from the GLODAPv2.2023 global interior ocean biogeochemical dataset
(https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0283442).

## Table of Contents
- [Introduction](#introduction)
- [Scripts Overview](#scripts-overview)
- [Usage](#usage)
- [Requirements](#requirements)
- [License](#license)
---
## Introduction
- This project is linked with the paper with the same tile, where methodological details are given, allowing
  a better understanding of the scripts and functions listed below and providing their statistical criteria. 
  Therefore, reading the Methods section of the paper is important to fully understand the explanation given below.  
- The main propose of this project was to evaluate the relative importance of physical mixing, 
  and biogeochemical processes, in determining nitrate concentrations in the central Arctic mesopelagic
- Scripts were implemented to: 
  1) retrieve data from GLODAPv2.2023, retaining only the part flagged
     as good quality data and separated by longitude: 
     - European (Atlantic) Arctic	45°W → 60°E		
     - Siberian Arctic	60°E → 170°E		
     - Pacific Arctic	170°E → 150°W
     - Canadian Arctic and Greenland	150°W → 45°W
     These sectors were further divided in two latitude bands: 70-80 and 80-90N
  2) select vertical profiles separately for each sector and latitude band, with data within the mesopelagic layer
     and in Atlantic Water (AW) or AW derived water masses, such as Modified Atlantic Water;
  3) analyse the relationships nitrate versus salinity and, AOU versus salinity, quantifying
     the fractions with significant linear positive and negative relationships 
     (conservative behavior - no significant biogeochemical sinks), the fractions with significant 
     quadratic “source or sink-type” – parabola opens downward (∩), the vertex is the highest point (a maximum) or upward (∪), 
     the vertex is the lowest point (a minimum), respectively, and the fraction with no significant
     linear or quadratic regressions.
  4) carry out statistical analysis to compare sectors.
- The work flow starts with exporting Glodap data for each to netCDF files, separately, 
  for each of the sectors listed above and using Ocean Data View software (see details below). Then, Script 1 is used to produce
  matlab files for further processing. Script 3 is then used to classify water masses and add them as a new variable to
  the matlabfiles. Scripts 3-5 analyse nitrate versus salinity and apparent oxygen utilization versus nitrate profiles,
  within Atlantic Water or its modified forms and at depths > 100 m. It classifies profiles as linear or quadratic according
  to criteria detailed in the paper. Then results are saved to a file for further processing, done with the Scripts 6-15. 
   

Ocean Data View details:
Version:	5.8.0 - 64 bit (Windows x64)
TEOS-10 Version:	3.2 - Apr/2023. Threadsafe C++ implementation by Reiner Schlitzer, AWI
User Directory:	C:\Users\pedro\Documents\ODV\
Browser:	C:\Program Files\Internet Explorer\IEXPLORE.EXE
Viewer:	notepad
User Name:	pedro
Host Name:	SARSTANGEN
Settings File:	C:\Users\pedro\.odv_settings_w64
Installation Path:	C:\Program Files\Ocean Data View\
Temp Directory:	C:\Users\pedro\Documents\ODV\.temp\ODV_pedro\2026-03-16T10-08-26\

---
## Scripts Overview
Here is a list of the independent scripts in this project, along with their descriptions:
### Script 1: **[ReadVars_from_netCDF_file.m](ReadVars_from_netCDF_file.m)**
- **Purpose**: Reading variables from netCDF files containing GLODAPv2.2023 data extracted for
               each sector. Selecting only good-quality data and storing the data in files
               with two main recrods  - one for each latitude range. The first lines
               of this script should be commented/uncommented depending on the 
               sector and latitude range being processed.
- **Inputs**: the following netCDF files (these input files are at the beginning of the
              script and should be commented/uncommented depending on which is to be used
              - only one file each time): 
                data_from_GLODAPv2.2023_Canadian_sector_70-80N.nc
                data_from_GLODAPv2.2023_Canadian_sector_80-90N.nc
                data_from_GLODAPv2.2023_European_sector_70-80N.nc
                data_from_GLODAPv2.2023_European_sector_80-90N.nc
                data_from_GLODAPv2.2023_Pacific_sector_70-80N_170-180E.nc
                data_from_GLODAPv2.2023_Pacific_sector_80-90N_170-180E.nc
                data_from_GLODAPv2.2023_Pacific_sector_70-80N_150-180W.nc
                data_from_GLODAPv2.2023_Pacific_sector_80-90N_150-180W.nc
                data_from_GLODAPv2.2023_Siberian_sector_70-80N.nc
                data_from_GLODAPv2.2023_Siberian_sector_80-90N.nc.
- **Outputs**: The outputs are matlab files containing selected data per sectors and 
               latitude band (one record per latitude band): 
               Barents_and_Arctic_European_Sector.m
               Canadian_Sector.mat
               Pacific_Sector_E.mat
               Pacific_Sector_W.mat
               Pacific_Sector.mat (this file resulted from merging the former two and it was used for further processing)
               Siberian_Sector.mat. 
- **How to Run**: Just typing the script name and making sure that the input files are accessible.

### Script 2: **[WaterMassClassification.m](WaterMassClassification.m)**
- **Purpose**: Reading hydrographic data from the matlab files listed above,
               classifying the water masses from which the results were obtained and
               adding the water masses as a new variable (type string). This was done
               to facilitate data filtering for specific water masses.
- **Inputs**: One of the matlabfiles listed above, commenting/uncommenting the
              relevant lines at the top of the script.
- **Outputs**: The outputs the same matlab files liosted above but containing a new variable
               with the water mass classification. 
- **How to Run**: Just typing the script name and making sure that the input files are accessible.

### Script 3: **[Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m] (Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m)**
- **Purpose**: Computes linear and curvilinear regressions between nitrate and salinity absolute
               and between apparent oxygen utilization and nitrate for each section, using 
               diagnose_linearity_adaptive_restrictive.m and diagnose_quadratic_structure.m. functions (see below).
- **Inputs**: All the matlab files listed above.
- **Outputs**: One matlab file with detailed results of the linear and quadratic regression analysis, per sector and laritude band
               This file lists results without any significant linear or quadratic regressions, results with significant linear
               regressions, where quadratic fits did not improve significantly the match between predictions and observations,
               and results with significant quadratic regressions. All regression parameters are also available. 
- **How to Run**: Just typing the script name and making sure that the input files are accessible. In this script
                  there is the following line: k = find(WM(:,j) >= 30 & WM(:,j) <= 50 & depth(:,j) >= UpperDepthLimit);
                  This is to select the available resuls for Atlantic Water types and obtained at 
                  depths > 100 m. This lines can be edited to add other selections. For example,
                  just a given period of time, based on UNIXSECS (referred to 1970) in the input matlab files.

### Script 4: **[Nitrate_vs_Salinity_vs_AOU_patterns_whole_data_v5.m] (Nitrate_vs_Salinity_vs_AOU_patterns_whole_data_v5.m)**
- **Purpose**: The same as the previous script but here data per sector and latitude band is agregatted before linear and quadratic
               regressions are calculated.
- **Inputs**: All the matlab files listed above.
- **Outputs**: One matlab file with detailed results of the linear and quadratic regression analysis per sector and latitude band.
- **How to Run**: Just typing the script name and making sure that the input files are accessible. In this script
                  there is the following line: k = find(WM(:,j) >= 30 & WM(:,j) <= 50 & depth(:,j) >= UpperDepthLimit);
                  This is to select the available resuls for Atlantic Water types and obtained at 
                  depths > 100 m. This lines can be edited to add other selections. For example,
                  just a given period of time, based on UNIXSECS (referred to 1970) in the input matlab files.

### Script 5: **[Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_with_mask_v1.m] (Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_with_mask_v1.m)**
- **Purpose**: Computes linear and curvilinear regressions between nitrate and salinity absolute
               and between apparent oxygen utilization and nitrate for each section and for 
               sea ice-covered and open water areas using 
               diagnose_linearity_adaptive_restrictive.m and diagnose_quadratic_structure.m functions (see below)
               and sea-ice masks for the period 1980-2025.
- **Inputs**: All the matlab files listed above and the sea ice masks that may be downloaded from: 
              https://doi.org/10.21334/NPOLAR.2026.807E1A75.
- **Outputs**: One matlab file with detailed results of the linear and quadratic regression analysis for each sector and 
               ice-covered and open-water areas. 
- **How to Run**: Just typing the script name and making sure that the input files are accessible. In this script
                  there is the following line: k = find(WM(:,j) >= 30 & WM(:,j) <= 50 & depth(:,j) >= UpperDepthLimit);
                  This is to select the available resuls for Atlantic Water types and obtained at 
                  depths > 100 m. This lines can be edited to add other selections. For example,
                  just a given period of time, based on UNIXSECS (referred to 1970) in the input matlab files.

### Script 6: **[summarize_linear_and_quadratic_NO3_SA_profiles.m] (summarize_linear_and_quadratic_NO3_SA_profiles.m)**
- **Purpose**: Summarizes the results of script Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m (see above).
- **Inputs**: The matlab file resulting from Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m.
- **Outputs**: A bar graph showing the absolute frequencies of different types of significant regressions.
               or their absence, for each sector and latitude band, between nitrate and salinity absolute. 
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Script 7: **[summarize_linear_and_quadratic_NO3_SA_profiles_rel_frequencies.m] (summarize_linear_and_quadratic_NO3_SA_profiles_rel_frequencies.m)**
- **Purpose**: Summarizes the results of script Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m (see above).
- **Inputs**: The matlab file resulting from Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m.
- **Outputs**: A bar graph showing the relative frequencies of different types of significant regressions.
               or their absence, for each sector and latitude band, between nitrate and salinity absolute. 
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Script 8: **[summarize_linear_and_quadratic_AOU_nitrate_profiles.m] (summarize_linear_and_quadratic_AOU_nitrate_profiles.m)**
- **Purpose**: Summarizes the results of script Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m (see above).
- **Inputs**: The matlab file resulting from Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m.
- **Outputs**: A bar graph showing the absolute frequencies of different types of significant regressions
               or their absence, for each sector and latitude band, between apparent oxygen utilization and nitrate. 
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Script 9: **[summarize_linear_and_quadratic_AOU_NO3_profiles_rel_frequencies.m] (summarize_linear_and_quadratic_AOU_NO3_profiles_rel_frequenciess.m)**
- **Purpose**: Summarizes the results of script Apparent Oxygen Utilization vs Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m (see above).
- **Inputs**: The matlab file resulting from Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v5.m.
- **Outputs**: A bar graph showing the relative frequencies of different types of significant regressions
               or their absence, for each sector and latitude band, between apparent oxygen utilization and nitrate. 
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Script 10: **[SlopeBoxPlots_v1.m] (SlopeBoxPlots_v1.m)**
- **Purpose**: Produce box plots for Deming slopes per geographic sector.
- **Inputs**: The matlab file resulting from Script 3.
- **Outputs**: A box plot with one box per sector for the Deming slopes observed at the latitude range 80-90N.
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Script 11: **[SlopeBoxPlots_v2.m] (SlopeBoxPlots_v2.m)**
- **Purpose**: Produce box plots for Deming slopes per geographic sector.
- **Inputs**: The matlab file resulting from Script 5.
- **Outputs**: A box plot with one box per sector for the Deming slopes observed under sea ice.
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Script 12: **[SlopeStatistics_v3.m] (SlopeStatistics_v3.m)**
- **Purpose**: Compute Kruskal Wallis and post-hoc tests to compare Deming slopes across geographic sectors.
- **Inputs**: The matlab file resulting from Script 3.
- **Outputs**: Significant tests for the null hypothesis of no slope differences between the various grographic sectors
               at the latitude range 80-90N.
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Script 13: **[SlopeStatistics_v4.m] (SlopeStatistics_v4.m)**
- **Purpose**: Compute Kruskal Wallis and post-hoc tests to compare Deming slopes across geographic sectors.
- **Inputs**: The matlab file resulting from Script 5.
- **Outputs**: Significant tests for the null hypothesis of no slope differences between the various grographic sectors
               in sea ice-covered areas.
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Script 14: **[Two_panel_script_for_linear_regressions.m] (Two_panel_script_for_linear_regressions.m)**
- **Purpose**: Produce a tile plot with a sample of linear nitrate versus salinity profiles per sector
- **Inputs**: The matlab file resulting from Script 3 and the sector files listed under Script 1.
- **Outputs**: Tile plot.
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Script 15: **[Two_panel_script_for_linear_regressions.m] (Two_panel_script_for_linear_regressions.m)**
- **Purpose**: Produce a tile plot with a sample of quadratic nitrate versus salinity profiles per sector
- **Inputs**: The matlab file resulting from Script 3 and the sector files listed under Script 1.
- **Outputs**: Tile plot.
- **How to Run**: Just typing the script name and making sure that the input files are accessible. 
                  Check first lines of the script to decide which input file to use.

### Function 1: **[diagnose_linearity_adaptive_restrictive.m] (diagnose_linearity_adaptive_restrictive.m)**
- **Purpose**: Computes Deming regressions using several criteria to evaluate their significance, which 
               are detailed in the accompanying paper (see also below). 
- **Inputs**: Function arguments which are two vectors with the values of x and y variables.
- **Outputs**: A record detailing the regression results.  
- **How to Run**: Just typing the script name followed by the arguments (Inputs, see above).

### Function 2: **[diagnose_quadratic_structure.m] (diagnose_quadratic_structure.m)**
- **Purpose**: Computes quadratic regressions using several criteria to evaluate their significance, which 
               are detailed in the accompanying paper (see also below). 
- **Inputs**: Function arguments which are two vectors with the values of x and y variables.
- **Outputs**: A record detailing the regression results.  
- **How to Run**: Just typing the script name followed by the arguments (Inputs, see above).





