function regenerate_event(LEFSystem, evheap, event_idx)

%     Regenerate an event in an event heap. If the event is currently impossible (e.g. a step
%     onto an occupied site) then the new event is not created, but the existing event is not
%     modified.
%
%     Possible events:
%     1 to N : SMC translocation
%     N+1 to 2N : SMC initialtion at gene promoter
%     2N+1 to 3N : SMC passive unbinding
%     3N+1 to 4N : SMC rebinding to a randomly chosen site


% A step to the left or to the right.
if event_idx <= LEFSystem.N
    dir = LEFSystem.smc_dir(event_idx);
    loop_idx = event_idx;
    if dir ==2
        leg_idx = loop_idx+LEFSystem.N;
    elseif dir == 1
        leg_idx = loop_idx;
    end
    if dir==1
        direction=-1;
    elseif dir==2
        direction=1;
    end
    
    local_vel = 0;
    if LEFSystem.smcs(leg_idx) > 0
        % Local velocity = velocity * permeability
        
        if LEFSystem.smc_dir(loop_idx) == 1
            local_vel = (LEFSystem.perms(LEFSystem.smcs(leg_idx) + (direction+1)/2)...
                * LEFSystem.vels(loop_idx));
        elseif LEFSystem.smc_dir(loop_idx) == 2
            local_vel = (LEFSystem.perms(LEFSystem.smcs(leg_idx) + (direction+1)/2)...
                * LEFSystem.vels(loop_idx));
        elseif LEFSystem.smc_dir(loop_idx) == 0
            local_vel = 0;
        end
        
        if local_vel > 0
            if LEFSystem.smclattice(LEFSystem.smcs(leg_idx) + direction) <=0
                evheap.add_event(event_idx,LEFSystem.time + exprnd(1/local_vel));
            end
        end
    end
    
    %SMC initialtion gene start
elseif event_idx > LEFSystem.N && event_idx <= 2 * LEFSystem.N
    loop_idx = event_idx - LEFSystem.N;
    if LEFSystem.smcs(loop_idx) > 0
        evheap.add_event(event_idx,LEFSystem.time +  exprnd(LEFSystem.rebinding_times(loop_idx)))
    end
    
    %SMC unbinding
elseif event_idx > 2 * LEFSystem.N && event_idx <= 3 * LEFSystem.N
    loop_idx = event_idx - 2 * LEFSystem.N;
    if LEFSystem.smcs(loop_idx) > 0
        evheap.add_event(event_idx,LEFSystem.time + exprnd(LEFSystem.lifesmcs(loop_idx)))
    end
    
    %SMC rebinding
elseif event_idx > 3 * LEFSystem.N && event_idx <= 4 * LEFSystem.N
    loop_idx = event_idx - 3 * LEFSystem.N;
    if LEFSystem.smcs(loop_idx) <= 0
        evheap.add_event(event_idx,LEFSystem.time + exprnd(LEFSystem.smc_times(loop_idx)))
    end
    
end

