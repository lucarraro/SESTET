# SESTET
Code and data supporting 
Carraro, L., Toffolon, M., Rinaldo, A., Bertuzzo, E. (2020). SESTET: a spatially explicit stream temperature model based on equilibrium temperature. *Hydrological Processes*, 34(2), 355-369. doi:10.1002/hyp.13591 

## Main files
- `RUN_AM.m` calibrates the temperature model via the Adaptive Metropolis algorithm. The model type can be chosen by appropriately setting variable `ModelType`; the calibration subset can be chosen by appropriately setting variable `SimType`.
- `RUN_PSO.m` calibrates the temperature model via the Particle Swarm Optimization algorithm. The model type can be chosen by appropriately setting variable `ModelType`; the calibration subset can be chosen by appropriately setting variable `SimType`.
- `MAIN.m` analyses model simulations and produces the same figures that are displayed in Carraro et al., (2019). Hydrol. Process. (hereafter: CHP)

Scripts `RUN_AM.m` and `RUN_PSO.m` generate output files whose generic name  is `*ModelType*_*SimType*_*AlgorithmType*.mat`. In order to be read by `MAIN.m`, these need to be copied into folder "results", and the `_*AlgorithmType*` part of their name (namely, `_AM` or `_PSO`) must be deleted.

## Other files
- `results`: contains .mat files produced by `RUN_AM.m` and `RUN_PSO.m`. 
- `utilities`: contains the following datasets needed to run the model:
  - `AirTemp_MeteoSuisse_Data.xlsx`: contains hourly air temperature data for the meteorological stations run by MeteoSuisse. The first column contains date and time, the following columns the air temperature time series in <sup>o</sup>C. Column names are acroynms used by MeteoSuisse to identify stations. The respective coordinates (in Swiss system LV1903) are contained in the script `EvalAirTemp.m`.
  - `OtherStations.xlsx`: as `AirTemp_MeteoSuisse_Data.xlsx`, but contains data from meteorological stations not directly run by MeteoSuisse (but rather from cantonal authorities, private companies etc.). These data were made available by Meteosuisse (IDAWEB portal).
  - `BestModels.mat`: this file is produced by `MAIN.m` and contains information on best simulations for each model type and calibration run. When `MAIN.m` is run, `BestModel.mat` is read if present; if absent, it is generated (this may take some minutes).
  - `data_aT_sT.xls`: time series of air and soil temperatures. for each triplet of columns (corresponding to a given station), the first one indicates date and time, the second hourly air temperature (in <sup>o</sup>C), the third the soil temperature at 1-m depth (in <sup>o</sup>C). Information on stations is embedded in script `EvalSoilTemp.m`. Data extracted from MeteoSuisse.
  - `Q_ZOF.txt`: mean daily discharge time series at station #3 (Zofingen). Data obtained from the Swiss Federal Office for the Environment (FOEN).
  - `stage-discharge.xlsx`: stage-discharge relationship for station #3. Data from FOEN.
  - `TempMeas.mat`: contains mean daily measured stream temperatures (in <sup>o</sup>C). Matrix `TempMeas` has one column per each station (column numbers correspond to station IDs as in Figure 2 of CHP. Vector `TimeMeas` contains corresponding dates in MATLAB format.
  - `DataWigger.mat`: morphological information on the catchment. The content of this dataset is described below.
- `DrawWiggerMap.m`: function that produces thematic maps of the Wigger catchment.
- `EvalAirTemp.m`: script that evaluates air temperature time series for the whole catchment. It is called by `RUN_AM.m`, `RUN_PSO.m`, `MAIN.m`. It cannot be run independently.
- EvalSoilTemp.m: script that evaluates soil temperature time series for the whole catchment. It is called by `RUN_AM.m`, `RUN_PSO.m`, `MAIN.m`. It cannot be run independently.
- HydraulicProperties.m: script that evaluates hydraulic properties for the whole catchment. It is called by `RUN_AM.m`, `RUN_PSO.m`, `MAIN.m`. It cannot be run independently.
- `interp.m`: function that performs fast 1-D interpolation (called by `SESTET_solver.m` and `SESTET_solver_weight.m`).
- `SESTET_solver.m`: function that runs models 'Sestet', 'Local' or 'Flat' on a given calibration set for a given parameter set.
- `SESTET_solver_weight.m`: as `SESTET_solver.m` but also stores the various heat fluxes making up the energy budget. It is considerably slower than `SESTET_solver.m`
- `v2struct.m`: function that packs and unpacks structures. It is used to import data into other functions in a synthetic fashion.

## Content of `DataWigger.m`
These data were originated from a 25-m digital elevation model of the region provided by the Swiss Federal Institute of Topography, analyzed via a TauDEM subroutine (D8 method) in ArcMap and subsequently post-processed in MATLAB.
- `AD`: adjacency matrix at the reach level.
- `AD_pixel`: adjacency matrix at the pixel level.
- `DEM`: digital elevation model of the region. `DEM(i,j)` is the elevation (in m a.s.l.) of the cell whose coordinates are (`XX(i,j)`, `YY(i,j)`).
- `N_reach`: number of reaches into which the river network is partitioned.
- `X`, `Y`: longitudinal/latitudinal coordinates (in Swiss LV1903 system) of all pixels constituting the river network.  
- `XX`, `YY`: longitudinal/latitudinal coordinates (in Swiss LV1903 system) of pixels constituting the DEM.
- `X_centr`, `Y_centr`: longitudinal/latitudinal coordinates (in Swiss LV1903 system) of subcatchment centroids.
- `Xc`, `Yc`: longitudinal/latitudinal coordinates (in Swiss LV1903 system) of pixels constituting the contour of the catchment.
- `area_local`: vector of subcatchment area values (in km<sup>2</sup>).
- `area_upstream`: vector of drainage area values for all subcatchments (in km<sup>2</sup>).
- `down_reach`: vector containing the ID of the reach downstream of the considered reach; if `down_reach(i)==0`, then `i` is the outlet reach.
- `length_reach`: vector containing reach length values (in m).
- `mean_subcatchment_altitude`: vector containing the mean elevation (in m a.s.l.) for each subcatchment.
- `nnodes`: number of pixels constituting the river network.
- `outlet`: ID of outlet pixel.
- `reach`: vector indicating the ID of the reach for each pixel constituting the river network.
- `reach_ID`: vector of ID of reaches where the stream temperature loggers were placed.
- `reach_slope`: vector of mean slope values for each reach.
- `reach_width`: vector of mean width values (in m) for each reach. 
- `subcatch`: matrix (with same dimension as `XX`, `YY`, `DEM`) containing the ID of the subcatchment to which every cell of the DEM belongs. If `isnan(subcatch(i,j))==1`, then (`XX(i,j)`, `YY(i,j)`) represents a pixel outside of the catchment. Note that every subcatchment corresponds to a single reach. 
