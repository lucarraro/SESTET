% coordinates of MeteoSuisse air temperature stations
StationNameMS=['BUS';'EGO';'GOE';'KOP';'MOA';'NAP';'PIL';'WYN'];
StationCoordMS=[648389 248365; 642913 225540; 640417 245937; 612662 218664; 660127 232851; 638132 206078; 661910 203410; 626400 233850]';
StationElevMS=[386; 521; 380; 484; 453; 1403; 2106; 422]';

% coordinates of other air temperature stations
StationCoordOther=[655842 226777; 665478 213381; 665102 215165; 628865 240180; 626820 231515; 622979 225749; 665543 209849; 628675 205310; 619752 243493; 639560 232110; 668855 202428; 648430 230838]';
StationElevOther=[797 485 426 430 457 480 454 1180 560 462 438 486];
StationNamesOther=['NABBRM';'INNEBI';'EMM   ';'NABHAE';'BELAT ';'MMLBB ';'LUZ   ';'BELUA ';'MMMUM ';'INNRED';'MMSTA ';'MMTRG '];

%% read air temperature data
% MeteoSuisse data
% adapt time in Matlab format
yy=floor(atms_tmp(:,1)/1e6); mm=floor((atms_tmp(:,1)-yy*1e6)/1e4); dd=floor((atms_tmp(:,1)-yy*1e6-mm*1e4)/1e2);
hh=floor(atms_tmp(:,1)-yy*1e6-mm*1e4-dd*1e2);
time_h=datenum(yy,mm,dd,hh,zeros(size(dd)),zeros(size(dd)));
% now let's read until Jun 30
%time_h(35065+90*24:end)=[]; % stop reading at 31 Mar 2018
time_h(1:3624)=[]; % start reading at 1 Jun 2014
%TempMeas_h=tmp(3625:(35064+90*24),2:end);
TempMeas_h=atms_tmp(3625:end,2:end);

TimeAirMS=time_h(1):floor(time_h(end));
for i=1:length(TimeAirMS)
    AirTempMeasMS(i,:)=mean(TempMeas_h((i-1)*24+1:24*i,:));
end
clear TempMeas_h atms_tmp

% other stations data
% adapt time in Matlab format
yy=floor(os_tmp(:,1)/1e6); mm=floor((os_tmp(:,1)-yy*1e6)/1e4); dd=floor((os_tmp(:,1)-yy*1e6-mm*1e4)/1e2);
hh=floor(os_tmp(:,1)-yy*1e6-mm*1e4-dd*1e2);
time_h=datenum(yy,mm,dd,hh,zeros(size(dd)),zeros(size(dd)));
% now let's read until Jun 30
%time_h(35065+90*24:end)=[]; % stop reading at 31 Mar 2018
time_h(1:3624)=[]; % start reading at 1 Jun 2014
%TempMeas_h=tmp(3625:(35064+90*24),2:end);
TempMeas_h=os_tmp(3625:end,2:end);

TimeAirOther=time_h(1):floor(time_h(end));
for i=1:length(TimeAirOther)
    AirTempMeasOther(i,:)=mean(TempMeas_h((i-1)*24+1:24*i,:));
end
TimeAir=TimeAirMS;
clear TempMeas_h os_tmp TimeAirOther

AirTempMeas=[AirTempMeasMS AirTempMeasOther];
efficiency=(length(TimeAir)-sum(isnan(AirTempMeas)))/length(TimeAir);

StationCoord=[StationCoordMS StationCoordOther];
StationElev=[StationElevMS StationElevOther];

%% calibration
if calib_aT
    dim=11; RMSE=zeros(dim); RMSE_k=zeros(dim,dim,length(StationElev)-1);
    a_vec=linspace(0,10,dim); b_vec=linspace(0,10,dim);
    for ind_a=1:dim
        for ind_b=1:dim
            a=a_vec(ind_a); b=b_vec(ind_b);
            for k=1:length(StationElev)
                TempNew=AirTempMeas; TempNew(:,k)=[];
                ElevNew=StationElev; ElevNew(k)=[];
                CoordNew=StationCoord; CoordNew(:,k)=[];
                Distance=0.001*sqrt((CoordNew(2,:)-StationCoord(2,k)).^2 + (CoordNew(1,:)-StationCoord(1,k)).^2); %[km]
                DeltaZ=abs(ElevNew-StationElev(k)); %[m]
                for t=1:length(TempNew)
                    TempMod_k(t,1)=nansum(TempNew(t,:)./(Distance.^a.*DeltaZ.^b))/sum(not(isnan(TempNew(t,:)))./(Distance.^a.*DeltaZ.^b));
                end
                if k~=7
                    k2=k; if k2>7; k2=k-1; end
                    RMSE_k(ind_a,ind_b,k2)=sqrt(nanmean((TempMod_k-AirTempMeas(:,k)).^2));
                end
                RMSE(ind_a,ind_b)=sum(squeeze(RMSE_k(ind_a,ind_b,:)).*efficiency([1:6 8:end])')/sum(efficiency([1:6 8:end]));
            end
        end
    end
    [ind_a,ind_b]=find(RMSE==min(min(RMSE)));
    dist_exp=a_vec(ind_a);
    alt_exp=b_vec(ind_b);
else
    dist_exp=2;
    alt_exp=2;
end

%% apply spatial interpolation
RMSE_k=[];
for k=1:length(StationElev)
    TempNew=AirTempMeas; TempNew(:,k)=[];
    ElevNew=StationElev; ElevNew(k)=[];
    CoordNew=StationCoord; CoordNew(:,k)=[];
    Distance=0.001*sqrt((CoordNew(2,:)-StationCoord(2,k)).^2 + (CoordNew(1,:)-StationCoord(1,k)).^2); %[km]
    DeltaZ=abs(ElevNew-StationElev(k)); %[m]
    for t=1:length(TempNew)
        TempMod_k(t,k)=nansum(TempNew(t,:)./(Distance.^dist_exp.*DeltaZ.^alt_exp))/sum(not(isnan(TempNew(t,:)))./(Distance.^dist_exp.*DeltaZ.^alt_exp));
    end
    RMSE_k(k)=sqrt(nanmean((TempMod_k(:,k)-AirTempMeas(:,k)).^2));
end

for i=1:N_reach
    Distance=0.001*sqrt((StationCoord(2,:)-Y_centr(i)).^2 + (StationCoord(1,:)-X_centr(i)).^2); %[km]
    DeltaZ=abs(StationElev-mean_subcatchment_altitude(i)); %[m]
    % if there is a NaN, weights change!
    MeanWeight(i,:)=1./(Distance.^dist_exp.*DeltaZ.^alt_exp)/sum(1./(Distance.^dist_exp.*DeltaZ.^alt_exp));
    for t=1:length(AirTempMeas)
        AirTemp(i,t)=nansum(AirTempMeas(t,:)./(Distance.^dist_exp.*DeltaZ.^alt_exp))/sum(not(isnan(AirTempMeas(t,:)))./(Distance.^dist_exp.*DeltaZ.^alt_exp));
    end
end

%% Test NSA model
R=8.3143; m=0.02897; Cp=1006; P0=100000;

DeltaA=zeros(length(TimeAir),length(StationElev)-1);
lambda=-0.0065; Tb=295;
for k=1:length(StationElev)
    TempNew=AirTempMeas; TempNew(:,k)=[];
    ElevNew=StationElev; ElevNew(k)=[];
    CoordNew=StationCoord; CoordNew(:,k)=[];
    
    Distance=0.001*sqrt((CoordNew(2,:)-StationCoord(2,k)).^2 + (CoordNew(1,:)-StationCoord(1,k)).^2); %[km]
    
    Pz=P0* (Tb./(Tb+lambda*ElevNew(:))).^(m*9.806/lambda/R);
    weights=(P0./Pz).^(R/m/Cp);
    for kk=1:length(StationElev)-1
        DeltaA(:,kk)=weights(kk)*(TempNew(:,kk)+273.15); % potential temperatures
    end
    for t=1:length(TempNew)
        DeltaA_k(t,:)=nansum(DeltaA(t,:)./(Distance.^2))/sum(not(isnan(DeltaA(t,:)))./(Distance.^2));
    end
    Pz_k=P0* (Tb/(Tb+lambda*StationElev(k)))^(m*9.806/lambda/R); weights_k=(P0/Pz_k)^(R/m/Cp);
    TempMod_NSA_k(:,k)=DeltaA_k/weights_k-273.15;
    RMSE_NSA_k(k)=sqrt(nanmean((TempMod_NSA_k(:,k)-AirTempMeas(:,k)).^2));
end

DeltaA=[]; DeltaA_k=[];
% Apply NSA
for k=1:N_reach
    Distance=0.001*sqrt((StationCoord(2,:)-Y_centr(k)).^2 + (StationCoord(1,:)-X_centr(k)).^2);
    Pz=P0* (Tb./(Tb+lambda*StationElev)).^(m*9.806/lambda/R);
    weights=(P0./Pz).^(R/m/Cp);
    for kk=1:length(StationElev)
        DeltaA(:,kk)=weights(kk)*(AirTempMeas(:,kk)+273.15); % potential temperatures
    end
    for t=1:length(AirTemp)
        DeltaA_k(t,:)=nansum(DeltaA(t,:)./(Distance.^2))/sum(not(isnan(DeltaA(t,:)))./(Distance.^2));
    end
    Pz_k=P0* (Tb/(Tb+lambda*mean_subcatchment_altitude(k)))^(m*9.806/lambda/R); weights_k=(P0/Pz_k)^(R/m/Cp);
    AirTemp_NSA(:,k)=DeltaA_k/weights_k-273.15;
end

%% Draw fig. S2
bin=[]; quants=[]; midpoints=[];
for i=1:length(StationElev)
    ind_subset=find(not(isnan(AirTempMeas(:,i))));
    bin(i,:)=quantile(AirTempMeas(ind_subset,i),[0:0.05:1]);
    midpoints(i,:)=[diff(bin(i,:))/2+bin(i,1:end-1)];
    for j=1:20
        low=bin(i,j); up=bin(i,j+1);
        tmp=find(and(AirTempMeas(:,i)>=low,AirTempMeas(:,i)<up));
        subset=TempMod_k(tmp,i);
        subsetNSA=TempMod_NSA_k(tmp,i);
        tmp=quantile(subset,[0.025 0.25 0.75 0.975]);
        tmpNSA=quantile(subsetNSA,[0.025 0.25 0.75 0.975]);
        for k=1:4; quants(i,j,k)=tmp(k); quantsNSA(i,j,k)=tmpNSA(k); end
    end
end

[~,ind_air]=sort(StationCoord(2,:),'descend');
if show_fig
    figure('units','centimeters','position',[0 0 18 22])
    for ind=1:20
        i=ind_air(ind);
        subplot(5,4,ind);
        hold on;
        patch([midpoints(i,:) flip(midpoints(i,:))],[quantsNSA(i,:,1) flip(quantsNSA(i,:,4))],[0.75 0 0],'edgecolor','none')
        patch([midpoints(i,:) flip(midpoints(i,:))],[quantsNSA(i,:,2) flip(quantsNSA(i,:,3))],[0.5 0 0],'edgecolor','none')
        patch([midpoints(i,:) flip(midpoints(i,:))],[quants(i,:,1) flip(quants(i,:,4))],[0.75 0.75 0.75],'edgecolor','none')
        patch([midpoints(i,:) flip(midpoints(i,:))],[quants(i,:,2) flip(quants(i,:,3))],[0.5 0.5 0.5],'edgecolor','none')
        hold on; plot([-15 30],[-15 30],'k'); title([num2str(ind),'  RMSE = ',num2str(RMSE_k(i))])
        set(gca,'tickdir','out','xtick',[-15:15:30],'ytick',[-15:15:30])
        axis([-15 30 -15 30]); box off
    end
end
