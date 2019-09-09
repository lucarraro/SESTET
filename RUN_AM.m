%% SESTET - Spatially-Explicit Stream Temperature model based on Equilibrium Temperature
% by Carraro, L. et al. (2019)

% RUN_AM.m performs calibration of the SESTET model on the case study
% dataset by means of an Adaptive Metropolis algorithm

clear all; close all; clc
rng('shuffle')

ModelType='Sestet'; % possible values: 'Sestet'; 'Local'; 'Flat'; 'Equil'
SimType='Three'; % possible values: 'All'; 'One'; 'Two'; 'Three'
show_fig=0;  % set to 1 tp calibrate model for air temperature
calib_aT=0;  % set to 1 to show figures on air/soil temperature fit

filename_short=[ModelType];
filename=[filename_short,'_',SimType,'_AM.mat'];

%% Define parameters
if strcmp(ModelType,'Sestet') || strcmp(ModelType,'local')
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
ParSigma=[0.03; 0.2; 0.2; 4; 0.2; 0.2]; % initial parameters' standard deviation

%% Load data
load('utilities\TempMeas.mat') % measured temperatures
load('utilities\Q_ZOF.txt')    % discharge time series
[sd_data,~,~]=xlsread('utilities\stage-discharge.xlsx'); % stage-discharge relationship
atms_tmp=xlsread('utilities\AirTemp_MeteoSuisse_Data.xlsx'); % air temperature data
os_tmp=xlsread('utilities\OtherStations.xlsx'); % other air temperature data
load('utilities\DataWigger.mat') % morphological data for the catchment

% air and soil temperature interpolation; evaluate hydraulic properties
EvalAirTemp;
EvalSoilTemp;
HydraulicProperties;

% define calibration set
SubsetAll=find(TimeAir==TimeMeas(1)):length(TimeAir);
if strcmp(SimType,'All')
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

%% INITIALIZE MARKOV CHAIN

% initial parameter values
a = 0.38;            % slope  for Teq
b = 6.55;              % [^oC] offset   for Teq
c = 1.9;
tau=195;
% parameters for heat exchange coefficient
k = 0.5;     % this is ln(k), reparametrization of k=Ke/rho/Cp [m/d] heat exchange coefficient;
delta=0.5; % this is logit of delta (reparametrization); real_delta=exp(delta)/(1+exp(delta));

N_run=20000;  sigmaLoglik=10;
indAcc=1;

q = @(COVmat,N_param) mvnrnd(zeros(1,N_param),COVmat);
COVmat = diag(ParSigma(1:N_param).^2);

for param=1:N_param
    eval(['ParOld(param)=',ParNames{param},';']);
    ParStruct.([ParNames{param}])(:,1)=ParOld(param);
end

ind_p=1;
for param=1:N_param
    eval([ParNames{param},'=ParOld(param);']);
end
K=exp(k); % reparameterization K=exp(k) enforces the heat exchange velocity to be positive
disp(sprintf('a = %.3f  -  b = %.3f  - c = %.3f  -  tau = %.1f  -  k = %.3f  -  delta = %.3f',a,b,c,tau,K,delta))

for t=1:length(TimeAir)
    reach_depth(:,t)=(area_upstream(:)/area_upstream(13)).^delta*d_ZOF(t);
    u(:,t)=Q_all(:,t)./reach_depth(:,t)./reach_width(:);
end
dDdt=[zeros(N_reach,1) diff(reach_depth,1,2)];

if strcmp(ModelType,'Equil')
    tic
    DL=cos(2*pi/365*(tau - [1:length(TimeAir)] - (31+28+31+30+31)));
    Teq=@(T,DL) a*T'+b+c*repmat(DL',1,N_reach);
    y=Teq(AirTemp,DL);
else
    params_Teq=v2struct(ModelType,a,b,c,tau,K);
    parameters = v2struct(Q_all,dDdt,u,Cp,g,N_reach,beta1,beta2,beta3,length_reach,reach_slope,reach_depth);
    tic
    [t,y,weight_Teq,weight_input,weight_lat,weight_frict,weight_dQdt]=...
        SESTET_solver(parameters,params_Teq,AirTemp,SoilTemp,[1:length(AirTemp)],ones(N_reach,1));
end
% Eval RMSE for calibration set
tmp=TempMeas(1:length(SubsetCalib),StationsCalib)-y(SubsetCalib,reach_ID(StationsCalib));
tmp=tmp(:);
RMSE_cal(ind_p,1)=sqrt(nanmean(tmp.^2));
Loglik(ind_p,1)=nansum(log(normpdf(tmp,0,sigmaLoglik)));
Loglik_old=Loglik;

% Eval RMSE for validation set
tmp=TempMeas(length(SubsetCalib)+1:end,StationsCalib)-y(SubsetValid,reach_ID(StationsCalib));
tmp=tmp(:);
tmp2=TempMeas(:,StationsValid)-y(SubsetAll,reach_ID(StationsValid));
tmp2=tmp2(:); tmp3=[tmp; tmp2];
RMSE_val(ind_p,1)=sqrt(nanmean(tmp3.^2));

disp(sprintf('run: %d  -  time %.2f s  -  loglik  =  %.1f  -  RMSE cal = %.2f',...
    ind_p,toc,Loglik(ind_p,1),RMSE_cal(ind_p,1)))

%% ADAPTIVE METROPOLIS ALGORITHM

for ind_p=2:N_run
    if mod(indAcc,10)==0
        par_mat=zeros(indAcc,N_param);
        for param=1:N_param
            par_mat(:,param) = ParStruct.([ParNames{param}])';
        end
        COVmat=(2.38/sqrt(N_param))^2 * cov(par_mat(1:indAcc-1,1:N_param)) + 1e-4 * eye(N_param);
    end
    ParNew = ParOld + q(COVmat,N_param);
    
    % "fold" value of delta if boundaries are exceeded
    % (see Vrugt, J.A., Markov chain Monte Carlo Simulation Using the DREAM
    % Software Package: Theory, Concepts, and MATLAB Implementation,
    % Environmental Modelling and Software, 2016)
    if N_param==6
        if ParNew(6)<0
            ParNew(6)=1+ParNew(6);
        elseif ParNew(6)>1
            ParNew(6)=ParNew(6)-1;
        end
    end
    
    for param=1:N_param
        eval([ParNames{param},'=ParNew(param);']);
    end
    K=exp(k); % reparameterization K=exp(k) enforces the heat exchange velocity to be positive
    disp(sprintf('a = %.3f  -  b = %.3f  - c = %.3f  -  tau = %.1f  -  k = %.3f  -  delta = %.3f',a,b,c,tau,K,delta))
    
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
        params_Teq=v2struct(ModelType,a,b,c,tau,K);
        parameters = v2struct(Q_all,dDdt,u,Cp,g,N_reach,beta1,beta2,beta3,length_reach,reach_slope,reach_depth);
        tic
        [t,y,weight_Teq,weight_input,weight_lat,weight_frict,weight_dQdt]=...
            SESTET_solver(parameters,params_Teq,AirTemp,SoilTemp,[1:length(AirTemp)],ones(N_reach,1));
    end
    % Eval Loglik and RMSE for calibration set
    tmp=TempMeas(SubsetCalib-27,StationsCalib)-y(SubsetCalib,reach_ID(StationsCalib));
    tmp=tmp(:); RMSE_cal_new=sqrt(nanmean(tmp.^2));
    Loglik_new=nansum(log(normpdf(tmp,0,sigmaLoglik)));
    % Eval RMSE for validation set
    tmp=TempMeas(SubsetValid-27,StationsCalib)-y(SubsetValid,reach_ID(StationsCalib));
    tmp=tmp(:);
    tmp2=TempMeas(:,StationsValid)-y(SubsetAll,reach_ID(StationsValid));
    tmp2=tmp2(:); tmp3=[tmp; tmp2];
    RMSE_val_new=sqrt(nanmean(tmp3.^2));
    if Loglik_new > Loglik_old || rand < exp(Loglik_new-Loglik_old)
        indAcc=indAcc+1;
        RMSE_cal(indAcc,1)=RMSE_cal_new; Loglik(indAcc,1)=Loglik_new;
        RMSE_val(indAcc,1)=RMSE_val_new;
        for param=1:N_param
            ParStruct.([ParNames{param}])(indAcc,1)=ParNew(param);
        end
        Loglik_old=Loglik_new; ParOld=ParNew;
        disp(sprintf(' ')); disp(sprintf('ACCEPTED!'))
        Loglik=Loglik;
        save(filename,'ParStruct','Loglik','RMSE_cal','RMSE_val')
    end
    disp(sprintf('run %d  -  acc %d  -  time %.1f s  -  loglik %.1f  -  RMSEcal %.2f  -  RMSEval %.2f',...
        ind_p,indAcc,toc,Loglik_new,RMSE_cal_new,RMSE_val_new))
end

%% SAVE RESULTS
% cut burn-in
for param=1:N_param
    ParStruct.(ParNames{param})(1:50)=[];
end
Loglik(1:50)=[]; RMSE_cal(1:50)=[]; RMSE_val(1:50)=[];

save(filename,'ParStruct','Loglik','RMSE_cal','RMSE_val')
