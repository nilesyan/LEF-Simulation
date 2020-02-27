function [MSD1,MSD2]=PolymerSim(l_sites_traj,r_sites_traj,ts_traj,L)
%  The chromatin is modeled as Rouse model chain.
%  Both the polymers without loop and with loops are simulated in 2D.
%  The loop evolution is taken from the loop extrusion simulation as input.
%  By computing the eigenmodes of the polymer, the exact positions of the 
%  beads are simulated.
%  The MSDs of different beads in the chain are then calculated.


rng(1) % generate random seeding
% Eliminate non-steady state
steady = 1.5e4;% Estimated start of steady state
l_sites_traj(1:steady,:)=[];% Left LEF trajactory
r_sites_traj(1:steady,:)=[];% Right LEF trajactory
ts_traj=ts_traj-ts_traj(steady);% Time
ts_traj(1:steady)=[];

% Parameters
zeta=9.47e-6; % friction coefficient
det=1; % time interval for MSD calculation
ted=2e4; % # of loop extrusion events to be included
ts=[ts_traj(1:ted)',ts_traj(1):det:ts_traj(ted)]'; % time series for simulation
ts(1:ted,2)=1:ted;
ts=sortrows(ts);
fst=find(ts(:,2)>0,1,'first');
for i=fst+1:size(ts,1)
    if ts(i,2)==0
        ts(i,2)=ts(i-1,2);
    end
end
ts(1,:)=[];
timeM=ts_traj(1):det:ts_traj(ted); % time series for MSD calculation

for loopin=1:2 % 1: with loops, 2: without loop
    for dmrouse=1:2 % dimension
        kb=1.38e-23;
        Temp=300;
        
        b=2e-7; % averaged bead-bead separation
        lk=138e-9; % estimated kuhn length
        Nk=b/lk;        
        NN=L; % number of beads
                
        k=2*kb*Temp/Nk/lk^2; % spring between beads
        k2=k; % additional spring due to SMC loops
        
        for indts=1:size(ts,1)
            np=indts+1;
            if indts==1
                dt=ts(1,1);
            else
                dt=ts(indts,1)-ts(indts-1,1);
            end
            
            VMu=ones(NN,1);
            VMu=VMu*(1/zeta); 
            MMu=diag(VMu,0); % diagonal friction matrix
            
            VKappa=ones(NN,1);
            VKappa=VKappa*(2*k);
            VKappaoff=ones(NN-1,1)*(-k);
            KKappa=diag(VKappa,0)+diag(VKappaoff,1)+diag(VKappaoff,-1); % spring matrix
            KKappa(1,1)=k; % free-end polymer
            KKappa(end,end)=k; % free-end polymer
            
            if loopin==1 % when loop is considered, spring matrix adds additional off-diagonal terms
                for kapi=1:N                    
                    if l_sites_traj(ts(indts,2),kapi)>0
                        KKappa(l_sites_traj(ts(indts,2),kapi),r_sites_traj(ts(indts,2),kapi))=...
                            KKappa(l_sites_traj(ts(indts,2),kapi),r_sites_traj(ts(indts,2),kapi))-1*k2;
                        KKappa(r_sites_traj(ts(indts,2),kapi),l_sites_traj(ts(indts,2),kapi))=...
                            KKappa(r_sites_traj(ts(indts,2),kapi),l_sites_traj(ts(indts,2),kapi))-1*k2;
                        KKappa(l_sites_traj(ts(indts,2),kapi),l_sites_traj(ts(indts,2),kapi))=...
                            KKappa(l_sites_traj(ts(indts,2),kapi),l_sites_traj(ts(indts,2),kapi))+1*k2;
                        KKappa(r_sites_traj(ts(indts,2),kapi),r_sites_traj(ts(indts,2),kapi))=...
                            KKappa(r_sites_traj(ts(indts,2),kapi),r_sites_traj(ts(indts,2),kapi))+1*k2;
                    end                   
                end
            end
            
            MMukappa=MMu*KKappa;
            [vecs,Lambda]=eig(MMukappa); 
            lambda=diag(Lambda);
            lambda(1)=0; % for free-end polymer, the 1st eigenvalue is 0
            
            for i=1:NN
                if vecs(1,i)<0
                    YY(:,i)=-vecs(:,i);
                else
                    YY(:,i)=vecs(:,i);
                end
                ss(i)=exp(-lambda(i)*dt);
                tt(i)=sqrt(1-ss(i)^2);
            end
                       
            Lcheck=YY'*KKappa*YY;             
            for i=1:NN
                uu(i)=sqrt(kb*Temp./Lcheck(i,i)); 
                UT(i)=uu(i)*tt(i);
            end
                        
            SS=diag(ss,0);
            
            % start with equilibrium position
            if np==2
                P1(1)=0;
                for iniPW=2:NN
                    P1(iniPW)=randn(1)*uu(iniPW);
                end
                PhyWalk(:,1)=YY*P1';
            end
            
            % physical positions of all the beads
            PhyWalk(:,np)=YY*SS*(inv(YY)*PhyWalk(:,np-1))+YY*diag(UT'*randn(1,NN));
                      
        end
        if loopin==1 && dmrouse==1
            dimPhyWalk11=PhyWalk;
        elseif loopin==1 && dmrouse==2
            dimPhyWalk12=PhyWalk;
        elseif loopin==2 && dmrouse==1
            dimPhyWalk21=PhyWalk;
        elseif loopin==2 && dmrouse==2
            dimPhyWalk22=PhyWalk;
        end
    end
    
end

[~,loc]=ismember(timeM,ts(:,1));
x1=dimPhyWalk11(:,loc);
y1=dimPhyWalk12(:,loc);
x2=dimPhyWalk21(:,loc);
y2=dimPhyWalk22(:,loc);

bn=[1:10:600]; % calculate MSD for 60 beads
for jj=1:length(bn)
    j=bn(jj);
    for i=1:length(timeM)-1
        MSD1(jj,i)=mean((x1(j,1:end-i)'-x1(j,1+i:end)').^2+(y1(j,1:end-i)'-y1(j,1+i:end)').^2,1);
    end
end

for jj=1:length(bn)
    j=bn(jj);
    for i=1:length(timeM)-1
        MSD2(jj,i)=mean((x2(j,1:end-i)'-x2(j,1+i:end)').^2+(y2(j,1:end-i)'-y2(j,1+i:end)').^2,1);
    end
end


end