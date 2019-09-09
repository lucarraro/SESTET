function [t,y,weight_Teq,weight_input,weight_lat,weight_frict,weight_dQdt,Teq_day]=SESTET_solver(parameters,params_Teq,aT,sT,tspan,y0)
%
[Q,dDdt,u,Cp,g,N_reach,beta1,beta2,beta3,length_reach,reach_slope,reach_depth]=v2struct(parameters);
[type_model,a,b,c,tau,Kconst]=v2struct(params_Teq);
ind_reach=1:N_reach;

Teq_day=zeros(N_reach,length(Q));

DL=cos(2*pi/365*(tau - [1:length(Q)] - (31+28+31+30+31)));
Teq=@(T,DL) a*T+b+c*DL;

% choose linear or constant heat exchange coefficient
 K=@(T,Q) Kconst;

%%%%%%options=odeset('NonNegative',[1:N_reach]);
if  strcmp(type_model,'Flat')
    [t,y]=ode45(@(t,y)eqs_45(t,y),tspan,y0);   
else
    [t,y]=ode15s(@(t,y)eqs_15(t,y),tspan,y0);   
end

    function dy=eqs_45(t,y)
        dy=zeros(N_reach,1);
        t1=floor(t-tspan(1))+1; t2=ceil(t-tspan(1))+1;
        if t2-t1<1
            airT=aT(ind_reach,t1);
            dayL=DL(t1);
            qQ=Q(:,t1);
            dD=reach_depth(:,t1);
        else
            airT=interp(aT(ind_reach,t1),aT(ind_reach,t2),t1,t2,t);
            dayL=interp(DL(t1),DL(t2),t1,t2,t);
            qQ=interp(Q(:,t1),Q(:,t2),t1,t2,t);
            dD=interp(reach_depth(:,t1),reach_depth(:,t2),t1,t2,t);
        end
        Teq_day=Teq(airT,dayL);
        % calculate Teq contribution
        if strcmp(type_model,'Flat')
            weight_Teq = (Teq_day-y(ind_reach))...
                .*K(airT,qQ);
        else
            weight_Teq = (Teq_day-y(ind_reach))...
                .*K(airT,qQ)./dD;
        end
        % calculate other contributions
        weight_input = 0;  weight_lat = 0;  weight_frict = 0;  weight_dQdt = 0;
        dy(ind_reach) = weight_Teq + weight_input + weight_lat + weight_frict + weight_dQdt; %
    end

    function dy=eqs_15(t,y)
        dy=zeros(N_reach,1);
        t1=floor(t-tspan(1))+1; t2=ceil(t-tspan(1))+1;
        if t2-t1<1
            airT=aT(ind_reach,t1);
            dayL=DL(t1);
            qQ=Q(:,t1);
            dD=reach_depth(:,t1);
            uU=u(:,t1);
            ddDdt=dDdt(:,t1);
            soilT=sT(ind_reach,t1);
        else
            airT=interp(aT(ind_reach,t1),aT(ind_reach,t2),t1,t2,t);
            dayL=interp(DL(t1),DL(t2),t1,t2,t);
            qQ=interp(Q(:,t1),Q(:,t2),t1,t2,t);
            dD=interp(reach_depth(:,t1),reach_depth(:,t2),t1,t2,t);
            uU=interp(u(:,t1),u(:,t2),t1,t2,t);
            ddDdt=interp(dDdt(:,t1),dDdt(:,t2),t1,t2,t);
            soilT=interp(sT(ind_reach,t1),sT(ind_reach,t2),t1,t2,t);
        end
        % calculate Teq contribution
        Teq_day=Teq(airT,dayL);
        weight_Teq = (Teq_day-y(ind_reach))...
            .*K(y(ind_reach),qQ)./dD;
        % calculate other contributions
        if strcmp(type_model,'Local')
        weight_input = 0;  weight_lat = 0;  weight_frict = 0;  weight_dQdt = 0;
        else
        weight_input = 86400*uU.*diag(beta1*(repmat(y,1,N_reach)-repmat(y,1,N_reach)'))./length_reach(ind_reach);
        weight_lat = 86400*uU.*beta2(ind_reach).*(soilT - y(ind_reach))./length_reach(ind_reach);
        weight_frict = 86400*uU.*beta3(ind_reach).*reach_slope(ind_reach)/g/Cp;
        weight_dQdt = ddDdt./dD.*(soilT - y(ind_reach));
        end
        dy(ind_reach) = weight_Teq + weight_input + weight_lat + weight_frict + weight_dQdt; %
    end

end
