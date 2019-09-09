%% SESTET - Spatially-Explicit Stream Temperature model based on Equilibrium Temperature
% by Carraro, L. et al. (2019)

% RUN_PSO.m performs calibration of the SESTET model on the case study
% dataset by means of a Particle Swarming Optimization algorithm

clear all; close all; clc
rng('shuffle')

ModelType='Local'; % possible values: 'Sestet'; 'Local'; 'Flat'; 'Equil'
SimType='One'; % possible values: 'All'; 'One'; 'Two'; 'Three'; 'All_NSA'
show_fig=0;  % set to 1 tp calibrate model for air temperature
calib_aT=0;  % set to 1 to show figures on air/soil temperature fit

filename_short=[ModelType];
filename=[filename_short,'_',SimType,'_PSO.mat'];

%% Define parameters
if strcmp(ModelType,'Sestet') || strcmp(ModelType,'Local')
    N_param=6;
elseif strcmp(ModelType,'Flat')
    N_param=5;
elseif strcmp(ModelType,'Equil')
    N_param=4;
end

ParNames= {
    'a'; % multiplies Ta in Teq
    'b'; % constant term in Teq
    'c'; % multiplies Dl in Teq
    'tau'; % ordinal day of max sin
    'k'; % constant k
    'delta'; % exponent for depth \propto area_upstream^delta [0; 1]
    };
ParMinVal=[0.2; 2; 0; 120; 0; 0];
ParMaxVal=[0.8; 8; 4; 240; 4; 1];

%% Load data
load('utilities\TempMeas.mat') % measured temperatures
load('utilities\Q_ZOF.txt')    % discharge time series
[sd_data,~,~]=xlsread('utilities\stage-discharge.xlsx'); % stage-discharge relationship
atms_tmp=xlsread('utilities\AirTemp_MeteoSuisse_Data.xlsx'); % air temperature data
os_tmp=xlsread('utilities\OtherStations.xlsx'); % other air temperature data
load('utilities\DataWigger.mat') % morphological data for the catchment

% air and soil temperature interpolation; evaluate hydraulic properties
EvalAirTemp;
if strcmp(SimType,'All_NSA')
    AirTemp=AirTemp_NSA';
end
EvalSoilTemp;
HydraulicProperties;

% define calibration set
SubsetAll=find(ismember(TimeAir,TimeMeas));
if strcmp(SimType,'All') || strcmp(SimType,'All_NSA')
    StationsValid=[];  SpanTimeCalib=1:length(TimeMeas);
elseif strcmp(SimType,'One')
    StationsValid=[2 6 9]; SpanTimeCalib=1:round(0.8*length(TimeMeas));
elseif strcmp(SimType,'Two')
    StationsValid=[1 5 8]; SpanTimeCalib=round(0.2*length(TimeMeas)):length(TimeMeas);
elseif strcmp(SimType,'Three')
    StationsValid=[]; SpanTimeCalib=1:round(0.6*length(TimeMeas));
end
SubsetCalib=SubsetAll(SpanTimeCalib);
StationsCalib=setdiff([1:11],StationsValid);
TempCalib=TempMeas(SpanTimeCalib,StationsCalib);
SubsetValid=setdiff(SubsetAll,SubsetCalib);

%% INITIAL PARAMETER VALUES

a = 0.6302;            % slope  for Teq
b = 4.0163;              % [^oC] offset   for Teq
c = 1;
tau=150;
k = 5;     % Ke/rho/Cp [m/d] heat exchange coefficient
delta=0.4;

%%
%*************************************************************************
% PARTICLE SWARM OPTIMIZATION
%*************************************************************************
N_part=16*N_param; N_run=100;
omega=0.825; phi_b=1.19; phi_g=1.19;

for param=1:N_param
    ParStruct.([ParNames{param}])=nan(N_part,N_run);
end
RMSE_cal = Inf(N_part,N_run); RMSE_val = Inf(N_part,N_run);

for param=1:N_param
    ParStruct.([ParNames{param}])(:,1)=rand(N_part,1)*(ParMaxVal(param)-ParMinVal(param))+ParMinVal(param);
end

disp('run = 1');
for ind_p=1:N_part
    for param=1:N_param
        eval([ParNames{param},'=ParStruct.([ParNames{param}])(ind_p,1);']);
    end
    if strcmp(ModelType,'Equil')
        tic
        
        DL=cos(2*pi/365*(tau - [1:length(TimeAir)] - (31+28+31+30+31)));
        Teq=@(T,DL) a*T'+b+c*repmat(DL',1,N_reach);
        
        y=Teq(AirTemp,DL);
    else
        for t=1:length(TimeAir)
            reach_depth(:,t)=(area_upstream(:)/area_upstream(13)).^delta*d_ZOF(t);
            u(:,t)=Q_all(:,t)./reach_depth(:,t)./reach_width(:);
        end
        dDdt=[zeros(N_reach,1) diff(reach_depth,1,2)];
        params_Teq=v2struct(ModelType,a,b,c,tau,k);
        parameters = v2struct(Q_all,dDdt,u,Cp,g,N_reach,beta1,beta2,beta3,length_reach,reach_slope,reach_depth);
        tic
        [t,y,weight_Teq,weight_input,weight_lat,weight_frict,weight_dQdt]=...
            SESTET_solver(parameters,params_Teq,AirTemp,SoilTemp,[1:length(AirTemp)],ones(N_reach,1));
    end
    % Eval RMSE for calibration set
    tmp=TempMeas(SubsetCalib-27,StationsCalib)-y(SubsetCalib,reach_ID(StationsCalib));
    RMSE_cal(ind_p,1)=sqrt(nanmean((tmp(:)).^2));
    % Eval RMSE for validation set
    tmp=TempMeas(SubsetValid-27,StationsCalib)-y(SubsetValid,reach_ID(StationsCalib));
    tmp=tmp(:);
    tmp2=TempMeas(:,StationsValid)-y(SubsetAll,reach_ID(StationsValid));
    tmp2=tmp2(:); tmp3=[tmp; tmp2];
    RMSE_val(ind_p,1)=sqrt(nanmean(tmp3.^2));
    best_p(ind_p,1) = find(RMSE_cal(ind_p,:)==min(RMSE_cal(ind_p,:)),1,'first');
    disp(sprintf('particle: %d  -  time %.2f s  -  RMSE cal = %.2f  -  RMSE val = %.2f',...
        ind_p,toc,RMSE_cal(ind_p,1),RMSE_val(ind_p,1)))
end
save(filename,'ParStruct','RMSE_cal','RMSE_val')

best_g = find(RMSE_cal==min(min(RMSE_cal)),1,'first'); % best swarm position

% initialize velocities
for param=1:N_param
    ParStruct.(['v_',ParNames{param}])=2*rand(N_part,1)*...
        (ParMaxVal(param)-ParMinVal(param))-(ParMaxVal(param)-ParMinVal(param));
end
disp(' ')

for ind_run = 2:N_run
    disp(sprintf('run = %d',ind_run))
    rnd=rand(N_param,2);
    % find parameter values corresponding to best p
    for ind_p = 1:N_part
        for param=1:N_param
            ParStruct.([ParNames{param},'_p'])(ind_p,1)=ParStruct.([ParNames{param}])(ind_p,best_p(ind_p));
        end
    end
    % find parameter values corresponding to best g, calculate velocity and
    % update parameters
    for param=1:N_param
        ParStruct.([ParNames{param},'_g'])=ParStruct.([ParNames{param}])(best_g);
        ParStruct.(['v_',ParNames{param}])=omega* ParStruct.(['v_',ParNames{param}]) +....
            phi_b*rnd(param,1)*(ParStruct.([ParNames{param},'_p'])-ParStruct.([ParNames{param}])(:,ind_run-1)) +...
            phi_g*rnd(param,2)*(ParStruct.([ParNames{param},'_g'])-ParStruct.([ParNames{param}])(:,ind_run-1));
        ParStruct.([ParNames{param}])(:,ind_run)=min(ParMaxVal(param),...
            max(ParMinVal(param),ParStruct.([ParNames{param}])(:,ind_run-1)+...
            ParStruct.(['v_',ParNames{param}])));
    end
    for ind_p=1:N_part
        for param=1:N_param
            eval([ParNames{param},'=ParStruct.([ParNames{param}])(ind_p,ind_run);']);
        end
        if strcmp(ModelType,'Equil')
            tic            
            DL=cos(2*pi/365*(tau - [1:length(TimeAir)] - (31+28+31+30+31)));
            Teq=@(T,DL) a*T'+b+c*repmat(DL',1,N_reach);            
            y=Teq(AirTemp,DL);
        else
            for t=1:length(TimeAir)
                reach_depth(:,t)=(area_upstream(:)/area_upstream(13)).^delta*d_ZOF(t);
                u(:,t)=Q_all(:,t)./reach_depth(:,t)./reach_width(:);
            end
            dDdt=[zeros(N_reach,1) diff(reach_depth,1,2)];
            params_Teq=v2struct(ModelType,a,b,c,tau,k);
            parameters = v2struct(Q_all,dDdt,u,Cp,g,N_reach,beta1,beta2,beta3,length_reach,reach_slope,reach_depth);
            tic
            [t,y,weight_Teq,weight_input,weight_lat,weight_frict,weight_dQdt]=...
                SESTET_solver(parameters,params_Teq,AirTemp,SoilTemp,[1:length(AirTemp)],ones(N_reach,1));
        end
        % Eval RMSE for calibration set
        tmp=TempMeas(SubsetCalib-27,StationsCalib)-y(SubsetCalib,reach_ID(StationsCalib));
        RMSE_cal(ind_p,ind_run)=sqrt(nanmean((tmp(:)).^2));
        % Eval RMSE for validation set
        tmp=TempMeas(SubsetValid-27,StationsCalib)-y(SubsetValid,reach_ID(StationsCalib));
        tmp=tmp(:);
        tmp2=TempMeas(:,StationsValid)-y(SubsetAll,reach_ID(StationsValid));
        tmp2=tmp2(:); tmp3=[tmp; tmp2];
        RMSE_val(ind_p,ind_run)=sqrt(nanmean(tmp3.^2));
        best_p(ind_p,1) = find(RMSE_cal(ind_p,:)==min(RMSE_cal(ind_p,:)),1,'first');
        disp(sprintf('particle: %d  -  time %.2f s  -  RMSE cal = %.2f  -  RMSE val = %.2f',...
            ind_p,toc,RMSE_cal(ind_p,ind_run),RMSE_val(ind_p,ind_run)))
    end
    save(filename,'ParStruct','RMSE_cal','RMSE_val')
    
    best_g = find(RMSE_cal==min(min(RMSE_cal)),1,'first'); % best swarm position
    disp(' ')
end
