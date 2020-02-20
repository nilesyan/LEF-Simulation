classdef LEFSystem < handle
    properties
        L
        N
        time
        vels
        rebinding_times
        perms
        smclattice
        smcs
        lifesmcs
        smc_times
        smc_dir
        DR
        DR2
        FST
        RST
        geneon                        
    end
    methods
        function obj = LEFSystem(L,N,vels,lifesmcs,rebinding_times,...
                smc_times,init_smcs,perms,BE,BE_perms,DR,DR2,FST,RST)
            
            % the occupancy of each site of the system. If -1, the site is unoccupied.
            obj.L = L;
            obj.N = N;
            obj.smclattice = -1 * ones(L,1);
            obj.smcs = -1 * ones(2*N,1);
            obj.smc_dir = -1 * ones(N,1); % direction of SMC (1:left, 2:right)
            obj.vels = vels;
            obj.lifesmcs = lifesmcs;
            obj.rebinding_times = rebinding_times;
            obj.smc_times = smc_times;
            obj.DR = DR;
            obj.DR2 = DR2;
            obj.FST = FST;
            obj.RST = RST;                       
            obj.geneon = -1 * ones(N,1); % 1 if SMC is translocating on the chromatin, else 0 or -1
            
            if isempty(BE)
                if isempty(perms)
                    obj.perms = ones(L+1,1);
                else
                    obj.perms=perms;
                end
            else
                if isempty(perms)
                    obj.perms = ones(L+1,1);
                else
                    obj.perms=perms;
                end
                for i = 1:length(BE)
                    obj.perms(BE(i)+1) = BE_perms(i);
                end                
            end            
            
            obj.perms(1)=0;
            obj.perms(end)=0;
            
            % Initialize non-random loops
            for i = 1:2*N
                if init_smcs(i) < 0
                    continue
                end
                obj.smcs(i) = init_smcs(i);
            end
        end
        
        function r = make_step(obj,leg_idx,direction)
            % The variable `direction` can only take values +1 or -1.
            new_pos = obj.locs(leg_idx) + direction;
            r = obj.move_leg(leg_idx, new_pos);
        end
        
        function r = move_leg(obj,leg_idx,new_pos)
            if new_pos > 0 && obj.lattice(new_pos) > 0
                r = 0;
                return
            end
            prev_pos = obj.locs(leg_idx);
            obj.locs(leg_idx) = new_pos;
            
            if prev_pos > 0
                if obj.lattice(prev_pos) <= 0
                    r = 0;
                    return
                end
                obj.lattice(prev_pos) = -1;
            end
            if new_pos > 0
                obj.lattice(new_pos) = leg_idx;
            end
            r = 1;
        end
        
        function r = make_smcstep(obj,leg_idx,direction)
            % The variable `direction` can only take values +1 or -1.
            new_pos = obj.smcs(leg_idx) + direction;
            r = obj.move_smc(leg_idx, new_pos);
        end
        
        function r = move_smc(obj,leg_idx,new_pos)
            if new_pos > 0 && obj.smclattice(new_pos) > 0
                r = 0;
                return
            end
            prev_pos = obj.smcs(leg_idx);
            obj.smcs(leg_idx) = new_pos;
            
            if prev_pos > 0
                if obj.smclattice(prev_pos) <= 0
                    r = 0;
                    return
                end
                obj.smclattice(prev_pos) = -1;
            end
            if new_pos > 0
                obj.smclattice(new_pos) = leg_idx;
            end
            r = 1;
        end
    end
end