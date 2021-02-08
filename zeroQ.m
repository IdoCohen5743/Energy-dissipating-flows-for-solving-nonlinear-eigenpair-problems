function [uoutput...%,dtc,dut
    ] = zeroQ(c,u,dTheta,lapi)
%  dt = 0.0001;

uoutput = u;



%% moment zero
Qu = -c*u+(u.^2)/2;
% if(abs(sum(Qu))<10)
%     return;
% end

%% first moment
%ut = 0;


sum_Qu = sum(Qu);
Tu =u*lapi;
dQu = -c+u;
mE = -sum_Qu*dQu;
mJ = -Tu;

%     ut = mEk + mJk;
%     dt = (dTheta^2)*min(norm(mE),norm(mJ))/4;
ut = mE - sum(mE.*mJ)*mJ/norm(mJ)^2;

if(abs(sum(dQu.*ut))<1e-13)
    return;
end
   
dt = -sum(Qu)/sum(dQu.*ut);
if(isnan(dt))
    return;
end

%     dt = 1/(norm(dQu)^2-sum(dQu.*mJ)^2/norm(mJ)^2);
%dtc = dt;
%     dut = ut(1:1:end-1)-ut(2:1:end);
%dut = dut/dTheta;


deltaU = dt*(ut);
if(~isnan(deltaU))
else
    oooo = 0;
end
u = u + deltaU;

uoutput = u;
