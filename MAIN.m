clear all; close all; clc

calib_aT=0; % set to 1 tp calibrate model for air temperature
show_fig=0; % set to 1 to show figures on air/soil temperature fit

%% Load data
load('utilities\TempMeas.mat') % measured temperatures
load('utilities\Q_ZOF.txt')    % discharge time series
[sd_data,~,~]=xlsread('utilities\stage-discharge.xlsx'); % stage-discharge relationship
atms_tmp=xlsread('utilities\AirTemp_MeteoSuisse_Data.xlsx'); % air temperature data
os_tmp=xlsread('utilities\OtherStations.xlsx'); % other air temperature data
load('utilities\DataWigger.mat') % morphological data for the catchment

%% Initialization
EvalAirTemp;
EvalSoilTemp;
HydraulicProperties;

%% Define model types and parameters
ParNames= {
    'a'; % multiplies Ta in Teq [0.2; 0.8]
    'b'; % constant term in Teq [2; 8]
    'c'; % multiplies Dl in Teq  [0; 4]
    'tau'; % ordinal day of max sin [120; 240]
    'k'; % constant k [0; 4] 
    'delta'; % exponent for depth \approx area_upstream^delta [0; 1]
    };
ParMinVal=[0.2; 2; 0; 120; 0; 0];
ParMaxVal=[0.8; 8; 4; 240; 4; 1];

ModelType={'Sestet';'Local';'Flat';'Equil'};
SimType={'All';'One';'Two';'Three';'All_NSA'};
N_par_vec=[6 6 5 4];

%% load model runs
for indModel=1:length(ModelType)
    for indSim=1:length(SimType)
        N_param=N_par_vec(indModel);
        eval(['load(''results\',ModelType{indModel},'_',SimType{indSim},'.mat'')'])
        if indSim==1
            eval([ModelType{indModel},'.',SimType{indSim},'.Loglik=Loglik;'])
        end
        eval([ModelType{indModel},'.',SimType{indSim},'.RMSE_cal=RMSE_cal;'])
        eval([ModelType{indModel},'.',SimType{indSim},'.RMSE_val=RMSE_val;'])
        for param=1:N_param
            eval([ModelType{indModel},'.',SimType{indSim},'.',ParNames{param},'=[ParStruct.',ParNames{param},'];']);
        end
        eval([ModelType{indModel},'.',SimType{indSim},...
            '.bestRMSEcal=min(min(',ModelType{indModel},'.',SimType{indSim},'.RMSE_cal));']);
        eval([ModelType{indModel},'.',SimType{indSim},...
            '.indBest=find(',ModelType{indModel},'.',SimType{indSim},'.RMSE_cal==',ModelType{indModel},'.',SimType{indSim},'.bestRMSEcal);']);
        for indPar=1:N_param
            eval([ModelType{indModel},'.',SimType{indSim},...
                '.BEST_',ParNames{indPar},'=',ModelType{indModel},'.',SimType{indSim},'.',ParNames{indPar},'(',ModelType{indModel},'.',SimType{indSim},'.indBest);']);
        end
        if indSim==1
            eval([ModelType{indModel},'.',SimType{indSim},...
                '.Deviance=-2*',ModelType{indModel},'.',SimType{indSim},'.Loglik;'])
            eval([ModelType{indModel},'.',SimType{indSim},...
                '.DIC=mean(',ModelType{indModel},'.',SimType{indSim},'.Deviance)+0.5*var(',ModelType{indModel},'.',SimType{indSim},'.Deviance);'])
        else
            eval([ModelType{indModel},'.',SimType{indSim},...
                '.AIC=2*N_param+N_obs*log(min(',ModelType{indModel},'.',SimType{indSim},'.bestRMSEcal)^2);'])
            eval([ModelType{indModel},'.',SimType{indSim},...
                '.BIC=log(N_obs)*N_param+N_obs*log(min(',ModelType{indModel},'.',SimType{indSim},'.bestRMSEcal)^2);'])
        end
    end
end
clear Loglik ParStruct RMSE_cal RMSE_val

%% run best models
try load('utilities\BestModels.mat')
catch
    for indModel=1:length(ModelType)
        for indSim=1:length(SimType)
            disp(sprintf('%s %s...',ModelType{indModel},SimType{indSim}))            
            %K=1; delta=0.4;
            for i=1:N_par_vec(indModel)
                eval([ParNames{i},'=',ModelType{indModel},'.',SimType{indSim},'.BEST_',ParNames{i},';'])
            end 
            if strcmp(ModelType{indModel},'Equil')
                tic
                DL=cos(2*pi/365*(tau - [1:length(TimeAir)] - (31+28+31+30+31)));
                Teq=@(T,DL) a*T'+b+c*repmat(DL',1,N_reach);
                y=Teq(AirTemp,DL);
            else
                for t=1:length(TimeAir)
                    reach_depth(:,t)=(area_upstream(:)/area_upstream(13)).^delta*d_ZOF(t);
                    u(:,t)=Q_all(:,t)./reach_depth(:,t)./reach_width(:);
                    reach_V(:,t)=reach_depth(:,t).*reach_width.*length_reach;
                end
                dDdt=[zeros(N_reach,1) diff(reach_depth,1,2)];
                str=ModelType{indModel};
                if strcmp(SimType{indSim},'All')
                    k=exp(k);
                end
                params_Teq=v2struct(str,a,b,c,tau,k);
                parameters = v2struct(Q_all,dDdt,u,Cp,g,N_reach,beta1,beta2,beta3,length_reach,reach_slope,reach_depth);
                if strcmp(ModelType{indModel},'Sestet') && strcmp(SimType{indSim},'All')
                    tic
                    [t,y,weight_Teq,weight_input,weight_lat,weight_frict,weight_dQdt,tt]=...
                        SESTET_solver_weight(parameters,params_Teq,AirTemp,SoilTemp,[1:length(AirTemp)],ones(N_reach,1));                      
                    toc
                    for indTime=1:length(TimeAir)
                        tmp=and(tt>indTime-0.5,tt<indTime+0.5);
                        Sestet.All.WeightTeq(:,indTime)=mean(weight_Teq(:,tmp),2);
                        Sestet.All.WeightInput(:,indTime)=mean(weight_input(:,tmp),2);
                        Sestet.All.WeightLat(:,indTime)=mean(weight_lat(:,tmp),2);
                        Sestet.All.WeightFrict(:,indTime)=mean(weight_frict(:,tmp),2);
                        Sestet.All.WeightDQ(:,indTime)=mean(weight_dQdt(:,tmp),2);
                    end
                    clear weight_Teq weight_input weight_lat weight_frict weight_dQdt tt
                elseif strcmp(SimType{indSim},'All_NSA')
                    tic
                    [t,y,weight_Teq,weight_input,weight_lat,weight_frict,weight_dQdt]=...
                        SESTET_solver(parameters,params_Teq,AirTemp_NSA',SoilTemp,[1:length(AirTemp)],ones(N_reach,1));
                    toc
                else
                    tic
                    [t,y,weight_Teq,weight_input,weight_lat,weight_frict,weight_dQdt]=...
                        SESTET_solver(parameters,params_Teq,AirTemp,SoilTemp,[1:length(AirTemp)],ones(N_reach,1));
                    toc
                end
            end
            eval([ModelType{indModel},'.',SimType{indSim},'.V=reach_V;'])
            eval([ModelType{indModel},'.',SimType{indSim},'.BestTemp=y;'])
            tmp=TempMeas-y([1:length(TempMeas)]+27,reach_ID);
            eval([ModelType{indModel},'.',SimType{indSim},'.RMSEstat=sqrt(nanmean((tmp).^2));'])

            SubsetAll=find(TimeAir==TimeMeas(1)):length(TimeAir);
            if (strcmp(SimType{indSim},'All') || strcmp(SimType{indSim},'All_NSA'))
                StationsValid=[];  SpanTimeCalib=1:length(TimeMeas);
            elseif strcmp(SimType{indSim},'One')
                StationsValid=[2 6 9]; SpanTimeCalib=1:round(0.8*length(TimeMeas));
            elseif strcmp(SimType{indSim},'Two')
                StationsValid=[1 5 8]; SpanTimeCalib=round(0.2*length(TimeMeas)):length(TimeMeas);
            elseif strcmp(SimType{indSim},'Three')
                StationsValid=[]; SpanTimeCalib=1:round(0.6*length(TimeMeas));
            end
            SubsetCalib=SubsetAll(SpanTimeCalib);
            StationsCalib=setdiff([1:11],StationsValid);
            TempCalib=TempMeas(SpanTimeCalib,StationsCalib);
            SubsetValid=setdiff(SubsetAll,SubsetCalib);
            
            % calib set
            obs=TempMeas(SubsetCalib-27,StationsCalib);
            tmp=obs-y(SubsetCalib,reach_ID(StationsCalib));
            tmp=tmp(:); obs=obs(:);           
            
            % valid set
            obs=TempMeas(SubsetValid-27,StationsCalib);
            tmp=obs-y(SubsetValid,reach_ID(StationsCalib));
            tmp=tmp(:); obs=obs(:);
            obs2=TempMeas(:,StationsValid);
            tmp2=obs2-y(SubsetAll,reach_ID(StationsValid));
            tmp2=tmp2(:); obs2=obs2(:); tmp3=[tmp; tmp2]; obs3=[obs; obs2];
        end
    end
    save('utilities/BestModels.mat','Sestet','Local','Flat','Equil')
end

%% plot posterior distribution
figure('units','centimeters','position',[0 0 30 20])
for indModel=1:length(ModelType)
    subplot(4,7,(indModel-1)*7+1); hold on
    eval(['histogram(',ModelType{indModel},'.All.RMSE_cal,[0.85:0.3/20:1.15],''normalization'',''probability'')'])
    box off; eval('ylabel(ModelType{indModel})')
    if indModel==1; title('RMSE'); end
    set(gca,'ylim',[0 1],'tickdir','out','xlim',[0.85 1.15],'xtick',[0.85:0.15:1.15])
    for indPar=1:N_par_vec(indModel)
        subplot(4,7,(indModel-1)*7+1+indPar); hold on
        if indPar~=5
            eval(['histogram(',ModelType{indModel},'.All.',ParNames{indPar},',[ParMinVal(indPar):0.05*(ParMaxVal(indPar)-ParMinVal(indPar)):ParMaxVal(indPar)],''normalization'',''probability'')'])
        else
            eval(['histogram(exp(',ModelType{indModel},'.All.',ParNames{indPar},'),[ParMinVal(indPar):0.05*(ParMaxVal(indPar)-ParMinVal(indPar)):ParMaxVal(indPar)],''normalization'',''probability'')'])
        end
        box off; set(gca,'ylim',[0 1],'xlim',[ParMinVal(indPar) ParMaxVal(indPar)],'tickdir','out')
        if indPar==1; set(gca,'xtick',[0.2:0.3:0.8]); end; if indPar==2; set(gca,'xtick',[2:3:8]); end; if indPar==4; set(gca,'xtick',[120:30:240]); end
        for indSim=1:length(SimType)
            if indSim==1 && indPar==5
                eval(['plot([exp(',ModelType{indModel},'.',SimType{indSim},'.BEST_',ParNames{indPar},') exp(',...
                    ModelType{indModel},'.',SimType{indSim},'.BEST_',ParNames{indPar},')],[0 1])'])
            else
                eval(['plot([',ModelType{indModel},'.',SimType{indSim},'.BEST_',ParNames{indPar},' ',...
                    ModelType{indModel},'.',SimType{indSim},'.BEST_',ParNames{indPar},'],[0 1])'])
            end
        end
        if indModel==1; title(ParNames(indPar)); end
    end
end

%% Show time series of heat fluxes
lag=15;
% evaluate equilibrium temperature
DL=cos(2*pi/365*(Sestet.All.BEST_tau - [1:length(Q_outlet)] - (31+28+31+30+31)));
Teq=@(T,DL) Sestet.All.BEST_a*T + Sestet.All.BEST_b + Sestet.All.BEST_c*DL;
EqTemp=Teq(AirTemp,DL);

reach1=reach_ID(1); reach2=reach_ID(11);
figure('units','centimeters','position',[0 0 22 10])
subplot(2,2,1)
plot(TimeAir',movavg(1e-6*1000*Cp*Sestet.All.V(reach1,:)'.*Sestet.All.WeightFrict(reach1,:)'/86400/length_reach(reach1)/reach_width(reach1),'simple',lag),'color',[0.4 0.4 0.4]); hold on; box off; grid on
plot(TimeAir',movavg(1e-6*1000*Cp*Sestet.All.V(reach1,:)'.*Sestet.All.WeightTeq(reach1,:)'/86400/length_reach(reach1)/reach_width(reach1),'simple',lag),'color',[70/255 148/255 209/255]);
plot(TimeAir',movavg(1e-6*1000*Cp*Sestet.All.V(reach1,:)'.*Sestet.All.WeightInput(reach1,:)'/86400/length_reach(reach1)/reach_width(reach1),'simple',lag),'color',[2/255 87/255 59/255]);
plot(TimeAir',movavg(1e-6*1000*Cp*Sestet.All.V(reach1,:)'.*(Sestet.All.WeightLat(reach1,:)'+Sestet.All.WeightDQ(reach1,:)')/86400./length_reach(reach1)/reach_width(reach1),'simple',lag),'color',[240/255 127/255 26/255]);
set(gca,'tickdir','out','xlim',[TimeAir(28) TimeAir(end)],'ylim',[-1.1e-4 1.1e-4],'ytick',[-1e-4:5e-5:1e-4],...
    'xtick',[datenum(2015,01,01) datenum(2016,01,01) datenum(2017,01,01) datenum(2018,01,01)],'xticklabel',[])
ylabel('Heat fluxes [MW m^{-2}]'); title(num2str(reach1))
subplot(2,2,2)
plot(TimeAir',movavg(1e-6*1000*Cp*Sestet.All.V(reach2,:)'.*Sestet.All.WeightFrict(reach2,:)','simple',lag)/86400/length_reach(reach2)/reach_width(reach2),'color',[0.4 0.4 0.4]); hold on; box off; grid on
plot(TimeAir',movavg(1e-6*1000*Cp*Sestet.All.V(reach2,:)'.*Sestet.All.WeightTeq(reach2,:)','simple',lag)/86400/length_reach(reach2)/reach_width(reach2),'color',[70/255 148/255 209/255]);
plot(TimeAir',movavg(1e-6*1000*Cp*Sestet.All.V(reach2,:)'.*Sestet.All.WeightInput(reach2,:)','simple',lag)/86400/length_reach(reach2)/reach_width(reach2),'color',[2/255 87/255 59/255]);
plot(TimeAir',movavg(1e-6*1000*Cp*Sestet.All.V(reach2,:)'.*(Sestet.All.WeightLat(reach2,:)'+Sestet.All.WeightDQ(reach2,:)')/86400/length_reach(reach2)/reach_width(reach2),'simple',lag),'color',[240/255 127/255 26/255]);
set(gca,'tickdir','out','xlim',[TimeAir(28) TimeAir(end)],'ylim',[-1.1e-4 1.1e-4],'ytick',[-1e-4:5e-5:1e-4],...
    'xtick',[datenum(2015,01,01) datenum(2016,01,01) datenum(2017,01,01) datenum(2018,01,01)],'xticklabel',[],'yticklabel',[]); title(num2str(reach2))
subplot(2,2,3)
plot(TimeAir',movavg(SoilTemp(reach1,:)','simple',lag),'color',[0.7 0.2 0]); hold on; box off; grid on
plot(TimeAir',movavg(EqTemp(reach1,:)','simple',lag),'color',[0.2 0.7 0]);
plot(TimeAir',movavg(Sestet.All.BestTemp(:,reach1),'simple',lag),'color',[0 0.2 0.7]);
set(gca,'tickdir','out','xlim',[TimeAir(28) TimeAir(end)],'ylim',[0 20],...
    'xtick',[datenum(2015,01,01) datenum(2016,01,01) datenum(2017,01,01) datenum(2018,01,01)])
datetick('x','mmm-yy','keeplimits','keepticks'); ylabel('Temperature [^oC]');
subplot(2,2,4)
plot(TimeAir',movavg(SoilTemp(reach2,:)','simple',lag),'color',[0.7 0.2 0]); hold on; box off; grid on
plot(TimeAir',movavg(EqTemp(reach2,:)','simple',lag),'color',[0.2 0.7 0]);
plot(TimeAir',movavg(Sestet.All.BestTemp(:,reach2),'simple',lag),'color',[0 0.2 0.7]);
set(gca,'tickdir','out','xlim',[TimeAir(28) TimeAir(end)],'ylim',[0 20],...
    'xtick',[datenum(2015,01,01) datenum(2016,01,01) datenum(2017,01,01) datenum(2018,01,01)],'yticklabel',[])
datetick('x','mmm-yy','keeplimits','keepticks'); 

%% plot temperatures in 2016
figure('units','centimeters','position',[0 0 28 14])
for i=1:11
    subplot(3,4,i)
    plot(TimeMeas,TempMeas(:,i),'.b'); hold on; box off; grid on
    plot(TimeAir,Sestet.All.BestTemp(:,reach_ID(i)),'r'); %plot(TimeAir,Local.BestTemp_from_Sestet(:,reach_ID(i)),'--','color',[0.5 0.5 0.5]);
    plot(TimeMeas,TempMeas(:,i)-Sestet.All.BestTemp(28:end,reach_ID(i)),'k'); hold on; box off
    set(gca,'tickdir','out','xlim',[datenum(2016,01,01) datenum(2017,01,01)],'ylim',[-5 20],'ytick',[-5:5:20],...
        'xtick',[datenum(2016,01,01):366/4:datenum(2017,01,01)])
    datetick('x','mmm','keeplimits','keepticks'); title(['#',num2str(i)])
end

%% fit linear model
for t=1:length(TimeAir)
    [ss]=polyfit(mean_subcatchment_altitude',Sestet.All.BestTemp(t,:),1);
    TempSlope(t,1)=ss(1); TempOffset(t,1)=ss(2);
end
for t=1:length(TimeMeas)
    if not(isnan(TempMeas(t,11)))
        [ss]=polyfit(mean_subcatchment_altitude(reach_ID)',TempMeas(t,:),1);
    else
        [ss]=polyfit(mean_subcatchment_altitude(reach_ID(1:10))',TempMeas(t,1:10),1);
    end
    TempMeasSlope(t,1)=ss(1); TempOffset(t,1)=ss(2);
end
for t=1:length(TimeAir)
    [ss]=polyfit(mean_subcatchment_altitude,AirTemp(:,t),1);
    AirTempSlope(t,1)=ss(1); AirTempOffset(t,1)=ss(2);
end
figure('units','centimeters','position',[0 0 25 15])
DayMin=708+27; DayMax=570+27;%find(TempSlope==max(TempSlope));
subplot(2,3,1);  hold on; plot(TimeAir,AirTempOffset,'color',[0.5 0.5 0.5]); plot(TimeAir,TempOffset,'color',[0 0.5 1]); box off;
plot([TimeAir(DayMin) TimeAir(DayMin)],[-15 30],'k'); plot([TimeAir(DayMax) TimeAir(DayMax)],[-15 30],'k')
set(gca,'tickdir','out','xlim',[TimeAir(28) TimeAir(end)],'ylim',[-15 30],'ytick',[-15:15:30],...
    'xtick',[datenum(2015,01,01) datenum(2016,01,01) datenum(2017,01,01) datenum(2018,01,01)])
datetick('x','mmm-yy','keeplimits','keepticks'); ylabel('Offset [^oC]')
subplot(2,3,2);  hold on; plot(TimeAir,1000*AirTempSlope,'color',[0.5 0.5 0.5]); plot(TimeAir,1000*TempSlope,'color',[0 0.5 1]); box off;
plot([TimeAir(DayMin) TimeAir(DayMin)],[-15 15],'k');  plot([TimeAir(DayMax) TimeAir(DayMax)],[-15 15],'k')
set(gca,'tickdir','out','xlim',[TimeAir(28) TimeAir(end)],'ylim',[-15 15],'ytick',[-15:5:15],...
    'xtick',[datenum(2015,01,01) datenum(2016,01,01) datenum(2017,01,01) datenum(2018,01,01)])
datetick('x','mmm-yy','keeplimits','keepticks'); ylabel('Slope [^oCkm^{-1}]')
subplot(2,3,4); hold on
histogram(AirTempSlope,[-0.015:0.03/30:0.015],'facecolor',[0.5 0.5 0.5],'normalization','probability'); hold on;
histogram(TempSlope,[-0.015:0.03/30:0.015],'facecolor',[0 0.5 1],'normalization','probability');
set(gca,'tickdir','out','ytick',[0:0.1:0.2],'xtick',[-0.015 0 0.015],'xlim',[-0.015 0.015],'ylim',[0 0.3],'ytick',[0:0.1:0.3]);
xlabel('Slope'); ylabel('Frequency [-]')
subplot(2,3,5); hold on
plot(mean_subcatchment_altitude,AirTemp(:,DayMin),'.','color',[0.5 0.5 0.5]);
plot(mean_subcatchment_altitude',Sestet.All.BestTemp(DayMin,:),'.','color',[0 0.5 1]); box off;
plot(StationElev,AirTempMeas(DayMin,:),'o','color',[0.2 0.2 0.2])
plot(mean_subcatchment_altitude(reach_ID),TempMeas(DayMin-27,:),'o','color',[0 0.25 0.75])
%title(datestr(TimeAir(DayMin)))
set(gca,'tickdir','out','xlim',[400 1500],'ylim',[5 20])
ylabel('Temperature [^oC]')
subplot(2,3,6);  hold on
plot(mean_subcatchment_altitude,AirTemp(:,DayMax),'.','color',[0.5 0.5 0.5]);
plot(mean_subcatchment_altitude',Sestet.All.BestTemp(DayMax,:),'.','color',[0 0.5 1]); box off;
plot(StationElev,AirTempMeas(DayMax,:),'o','color',[0.2 0.2 0.2])
plot(mean_subcatchment_altitude(reach_ID),TempMeas(DayMax-27,:),'o','color',[0 0.25 0.75])
% title(datestr(TimeAir(DayMax)))
xlabel('Elevation [m a.s.l.]'); ylabel('Temperature [^oC]')
set(gca,'tickdir','out','xlim',[400 1500],'ylim',[-15 10])

RMSE_day=sqrt(nanmean((TempMeas-Sestet.All.BestTemp(28:end,reach_ID)).^2,2));

[~,months]=datevec(TimeAir);
ThermalInversionWater_summer=sum(and(TempSlope(28:end)'>0,and(months(28:end)>=4,months(28:end)<=9)))/sum(and(months(28:end)>=4,months(28:end)<=9))
ThermalInversionWater_winter=sum(and(TempSlope(28:end)'>0,or(months(28:end)<4,months(28:end)>9)))/sum(or(months(28:end)<4,months(28:end)>9))
ThermalInversionAir_summer=sum(and(AirTempSlope(28:end)'>0,and(months(28:end)>=4,months(28:end)<=9)))/sum(and(months(28:end)>=4,months(28:end)<=9))
ThermalInversionAir_winter=sum(and(AirTempSlope(28:end)'>0,or(months(28:end)<4,months(28:end)>9)))/sum(or(months(28:end)<4,months(28:end)>9))


%% model performance table
for i=1:length(ModelType)
    eval(['DIC_all(i,1)=',ModelType{i},'.All.DIC;']);
    eval(['RMSE_all(i,1)=',ModelType{i},'.All.bestRMSEcal;']);
    
    eval(['AIC_one(i,1)=',ModelType{i},'.One.AIC;']);
    eval(['RMSEcal_one(i,1)=',ModelType{i},'.One.bestRMSEcal;']);
    eval(['RMSEval_one(i,1)=',ModelType{i},'.One.RMSE_val(',ModelType{i},'.One.indBest);']);
    
    eval(['AIC_two(i,1)=',ModelType{i},'.Two.AIC;']);
    eval(['RMSEcal_two(i,1)=',ModelType{i},'.Two.bestRMSEcal;']);
    eval(['RMSEval_two(i,1)=',ModelType{i},'.Two.RMSE_val(',ModelType{i},'.Two.indBest);']);

    eval(['AIC_three(i,1)=',ModelType{i},'.Three.AIC;']);
    eval(['RMSEcal_three(i,1)=',ModelType{i},'.Three.bestRMSEcal;']);
    eval(['RMSEval_three(i,1)=',ModelType{i},'.Three.RMSE_val(',ModelType{i},'.Three.indBest);']);
    
    eval(['AIC_NSA(i,1)=',ModelType{i},'.All_NSA.AIC;']);
    eval(['RMSEcal_NSA(i,1)=',ModelType{i},'.All_NSA.bestRMSEcal;']);
end

DeltaDIC_all=DIC_all-min(DIC_all); DeltaAIC_one=AIC_one-min(AIC_one);
DeltaAIC_two=AIC_two-min(AIC_two); DeltaAIC_three=AIC_three-min(AIC_three); DeltaAIC_NSA=AIC_NSA-min(AIC_NSA);
table(DeltaDIC_all,RMSE_all,...
    DeltaAIC_one,RMSEcal_one,RMSEval_one,...
    DeltaAIC_two,RMSEcal_two,RMSEval_two,...
    DeltaAIC_three,RMSEcal_three,RMSEval_three,...
    DeltaAIC_NSA,RMSEcal_NSA)

%% station-by-station RMSE table
for i=1:length(ModelType)
    for j=1:length(SimType)
        row=(i-1)*length(SimType)+j;
        eval(['RMSE_statbystat(row,:)=',ModelType{i},'.',SimType{j},'.RMSEstat;'])
    end
end
array2table(RMSE_statbystat)

%% best-fit parameters
a_all=zeros(length(ModelType)*length(SimType),1); b_all=a_all; c_all=a_all; tau_all=a_all; k_all=a_all; delta_all=a_all;
for i=1:length(ModelType)
    for j=1:length(SimType)
        row=(i-1)*length(SimType)+j;
        eval(['a_all(row,1)=',ModelType{i},'.',SimType{j},'.BEST_a;'])
        eval(['b_all(row,1)=',ModelType{i},'.',SimType{j},'.BEST_b;'])
        eval(['c_all(row,1)=',ModelType{i},'.',SimType{j},'.BEST_c;'])
        eval(['tau_all(row,1)=',ModelType{i},'.',SimType{j},'.BEST_tau;'])
        if i<4
            if j==1
                eval(['k_all(row,1)=exp(',ModelType{i},'.',SimType{j},'.BEST_k);'])
            else
                eval(['k_all(row,1)=',ModelType{i},'.',SimType{j},'.BEST_k;'])
            end
        else
            k_all(row,1)=NaN;
        end
        if i<3
            eval(['delta_all(row,1)=',ModelType{i},'.',SimType{j},'.BEST_delta;'])
        else
            delta_all(row,1)=NaN;
        end
    end
end

table(a_all,b_all,c_all,tau_all,k_all,delta_all)

%% maps of temperature differences with respect to reach 1

geometry=v2struct(X,Y,Xc,Yc,AD_pixel,nnodes,outlet,area_upstream,reach); % structure containing river network data used to produce maps

% path from highest/farthest headwater to outlet
j=166; path=j;
while j~=0
   j=down_reach(j);
   path=[path j];
end
path(path==0)=[];

SummerFallAvg=mean(Sestet.All.BestTemp(months>6,:));
WinterSpringAvg=mean(Sestet.All.BestTemp(months<7,:));
SummerFallDiff=SummerFallAvg-SummerFallAvg(reach_ID(1));
WinterSpringDiff=WinterSpringAvg-WinterSpringAvg(reach_ID(1));

SummerFallAirAvg=mean(AirTemp(:,months>6)');
WinterSpringAirAvg=mean(AirTemp(:,months<7)');
SummerFallAirDiff=SummerFallAirAvg-SummerFallAirAvg(reach_ID(1));
WinterSpringAirDiff=WinterSpringAirAvg-WinterSpringAirAvg(reach_ID(1));

DrawWiggerMap(SummerFallDiff,5,-5,'summer-fall differences','FR2',geometry,1,1);
DrawWiggerMap(WinterSpringDiff,5,-5,'winter-spring differences','FR2',geometry,1,1);

SummerFallEqAvg=mean(EqTemp(:,months>6)');
WinterSpringEqAvg=mean(EqTemp(:,months<7)');
SummerFallEqDiff=SummerFallAvg-SummerFallEqAvg;
WinterSpringEqDiff=WinterSpringAvg-WinterSpringEqAvg;

DrawWiggerMap(SummerFallEqDiff,5,-5,'summer-fall Eq stream-equilibrium differences','FR2',geometry,1,1);
DrawWiggerMap(WinterSpringEqDiff,5,-5,'winter-spring Eq stream-equilibrium differences','FR2',geometry,1,1);

SummerFallSoilAvg=mean(SoilTemp(:,months>6)');
WinterSpringSoilAvg=mean(SoilTemp(:,months<7)');
SummerFallSoilDiff=SummerFallAvg-SummerFallSoilAvg;
WinterSpringSoilDiff=WinterSpringAvg-WinterSpringSoilAvg;

DrawWiggerMap(SummerFallSoilDiff,5,-5,'summer-fall soil stream-equilibrium differences','FR2',geometry,1,1);
DrawWiggerMap(WinterSpringSoilDiff,5,-5,'winter-spring soil stream-equilibrium differences','FR2',geometry,1,1);


%% reach statistics
figure('units','centimeters','position',[0 0 20 6])
subplot(1,4,1); histogram(reach_width,10,'facecolor',[1 0.5 0],'normalization','probability')
set(gca,'tickdir','out','ylim',[0 0.5],'ytick',[0 0.25 0.5]','xlim',[0 12],'xtick',[0:6:12]); box off; title('Width [m]')
subplot(1,4,2); histogram(length_reach,10,'facecolor',[0.5 0 1],'normalization','probability')
set(gca,'tickdir','out','ylim',[0 0.5],'ytick',[0 0.25 0.5],'xlim',[0 5000],'xtick',[0:2500:5000]); box off; title('Length [m]')
subplot(1,4,3); histogram(area_local,10,'facecolor',[0.5 1 0],'normalization','probability')
set(gca,'tickdir','out','ylim',[0 0.5],'ytick',[0 0.25 0.5],'xlim',[0 10],'xtick',[0:5:10]); box off; title('Local Area [km^2]')
subplot(1,4,4); histogram(mean_subcatchment_altitude,[400:800/11:1200],'facecolor',[0 0.5 1],'normalization','probability')
set(gca,'tickdir','out','ylim',[0 0.5],'ytick',[0 0.25 0.5],'xlim',[400 1200],'xtick',[400:400:1200]); box off; title('Local Elevation [m a.s.l.]')

table([1:11]',reach_width(reach_ID),length_reach(reach_ID),area_local(reach_ID),mean_subcatchment_altitude(reach_ID),...
    'VariableNames',{'Reach','Width','Length','LocalArea','MeanAltitude'})
