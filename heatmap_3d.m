function hmap = heatmap_3d(L,N,ts_traj,smc_lsites_traj,smc_rsites_traj,avgn)
% Convert 1d loop configurations to 3d polymer distribution probabilities
% using the effective genomic distance and Gaussian chain.
% Predict time-averaged contact probabilites and generate the contact map.

Cor_NP=zeros(L,L);
Cor_NPa=Cor_NP;
for i = 1:avgn
    tm = ts_traj(end-1)/avgn*i;
    mid=find(ts_traj>=tm,1,'first');
    for j=1:L-1
        for k=j+1:L
            PT1=j;
            PT2=k;            
            l_foot = smc_lsites_traj(mid,:);
            r_foot = smc_rsites_traj(mid,:);            
            Cor_NP(j,k) = N_eff_TwoPoints(l_foot,r_foot,PT1,PT2,N);
        end
    end    
    Cor_M = Cor_NP+Cor_NP';
    Cor_M(isinf(Cor_M)) = 1e-6;
    Cor_M(Cor_M==0) = 1e-6;
    a = 3/2./(Cor_M).^2;
    Cor_M2 = erf(sqrt(a))-2*sqrt(a/pi).*exp(-a);    
    if isnan(Cor_M2)
        disp(i)
    end
    Cor_NPa = Cor_M2+Cor_NPa;
end
Cor_NPa = Cor_NPa/avgn;
Cor_NPa(logical(eye(L)))=0;
hmap=Cor_NPa;
end