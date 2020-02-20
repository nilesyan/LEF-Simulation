function []=Main(tra)

load('zfp508_ES.mat') 
% load the genomic data
% The genomic data is converted from the gene positions in the library
% FST: Forward DNA strain
% RST: Reversed DNA srain
% The first two columns are the positions of gene promoters and terminators
% (coverted with the desired resolution e.g. 10kb in mouse and 0.5kb in S. pombe)
% The 3rd and the 4th columns are not used in the current model.
% "zfp508_ES" is an example of a 6Mb region in the mouse genome.

% Parameter Setting Section
rng(tra) % generate random seeding
L = 600; % # of lattice
N = 30;  % # of SMC pairs
R_EXTEND = 3e-3; % SMC translocation velocity
R_SMC = 1.22e-6*ones(N,1); % SMC fall-off rate (i.e. 1/lifetime)
SMC_TIME = 2; % SMC residence time (time in solution)
REBINDING_TIME = 1.7; % RNA polymerase binding time 
%                       (i.e. SMC dwell time at promoters after initial binding)

INIT_SMCL_SITES = -1 * ones(N,1); % initiate the LEF positions
INIT_SMCR_SITES = -1 * ones(N,1); % -1 represent "unbound".
ACTIVATION_TIMES = zeros(N,1);
T_MAX = 5e2; % control parameter for events per step
N_SNAPSHOTS = 1e6; % number of simulation steps
BE = []; % boundary elements
BE_perms = []; % bermeability of boundary elements
DR(:,1) = [FST(:,1);RST(:,2)]; % FST: Forward strain; RST: Reversed strain 
DR(:,2) = [FST(:,4);RST(:,4)]; % RNA level (not used in the simple sliding model)
DR2(:,1) =  [FST(:,1);RST(:,2)];
DR2 = unique(DR2);
PERMS = []; % lattice permeability
verbose = 1;

% Program Initiation Section
VELS = ones(2*N,1)*R_EXTEND;
LIFESSMCS = zeros(N,1);
REBINDING_TIMES = zeros(N,1);
SMC_TIMES = zeros(N,1);

for i = 1:N
    if length(REBINDING_TIME) > 1
        REBINDING_TIMES(i) = REBINDING_TIME(i);
    else
        REBINDING_TIMES(i) = REBINDING_TIME;
    end
end

for i=1:N
    if length(R_SMC) > 1
        LIFESSMCS(i) = 1/R_SMC(i);
    else
        LIFESSMCS(i) = 1/R_SMC;
    end
    
    if length(SMC_TIME) > 1
        SMC_TIMES(i) = SMC_TIME(i);
    else
        SMC_TIMES(i) = SMC_TIME;
    end
end

INIT_SMCS = -1 * ones(2*N,1);

for i = 1:N
    INIT_SMCS(i) = INIT_SMCL_SITES(i);
    INIT_SMCS(N+i) = INIT_SMCR_SITES(i);
end

if ~isempty(PERMS) && length(PERMS) ~= L+1
    msgID = 'PERMS:Length';
    msg = 'The length of the provided array of permeabilities should be L+1.';
    Exception = MException(msgID,msg);
    throw(Exception)    
end

LEFSYSTEM = LEFSystem(L, N, VELS, LIFESSMCS, REBINDING_TIMES, ...
    SMC_TIMES, INIT_SMCS, PERMS, BE, BE_perms, DR, DR2, FST, RST);
LEFSYSTEM.time = 0;

%lef_sites_traj = zeros(N_SNAPSHOTS, N);
smc_lsites_traj = zeros(N_SNAPSHOTS, N);
smc_rsites_traj = zeros(N_SNAPSHOTS, N);
ts_traj = zeros(N_SNAPSHOTS, 1);

prev_snapshot_t = 0;
snapshot_idx = 1;

evheap = Event_heap();

% Move LEFs onto the lattice at the corresponding activations times.
% If the positions were predefined, initialize the fall-off time in the
% standard way.
for i = 1:LEFSYSTEM.N
    if INIT_SMCS(i) == -1
        evheap.add_event(i + 3 * LEFSYSTEM.N, ACTIVATION_TIMES(i));
    else
        regenerate_all_smc_events(LEFSYSTEM, evheap, i)
    end
end

% Simulation Section 
% - refer to do_event.m for event descriptions
while snapshot_idx <= N_SNAPSHOTS
    LEFEvent = evheap.pop_event();
    LEFSYSTEM.time = LEFEvent.time;
    event_idx = LEFEvent.event_idx;
    
    LEFStatus = do_event(LEFSYSTEM, evheap, event_idx);
    
    if LEFStatus == 0
        disp('an assertion failed somewhere')
        return
    end
    
    if LEFSYSTEM.time > prev_snapshot_t + T_MAX / N_SNAPSHOTS
        prev_snapshot_t = LEFSYSTEM.time;
        % lef_sites_traj(snapshot_idx,1:N) = LEFSYSTEM.locs(1:N);
        smc_lsites_traj(snapshot_idx,1:N) = LEFSYSTEM.smcs(1:N);
        smc_rsites_traj(snapshot_idx,1:N) = LEFSYSTEM.smcs(N+1:end);
        ts_traj(snapshot_idx) = LEFSYSTEM.time;
        
        snapshot_idx = snapshot_idx + 1;
        if verbose && mod(snapshot_idx,1e4) == 0
            disp([snapshot_idx/N_SNAPSHOTS])
        end
    end
end

% Hi-C Map Generation Section
intv = 1e3; % # of averaged maps
hmap = heatmap_3d(L,N,ts_traj,smc_lsites_traj,smc_rsites_traj,intv);

end
