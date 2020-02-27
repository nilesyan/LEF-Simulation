function sta = do_event(LEFSystem, evheap, event_idx)

%     Apply an event from a heap on the system and then regenerate it.
%     If the event is currently impossible (e.g. a step onto an occupied site),
%     it is not applied, however, no warning is raised.
%
%     Also, partially checks the system for consistency. Returns 0 if the system
%     is not consistent (a very bad sign), otherwise returns 1 if the event was a step
%     and 2 if the event was rebinding.
%
%     Possible events:
%     1 to N : LEF translocation
%     N+1 to 2N : LEF initialtion at gene promoter
%     2N+1 to 3N : LEF passive unbinding
%     3N+1 to 4N : LEF rebinding to a randomly chosen site

    %LEF translocation
if event_idx <= LEFSystem.N        
    dir = LEFSystem.smc_dir(event_idx);
    loop_idx = event_idx;
    if dir ==2
        leg_idx = loop_idx+LEFSystem.N;
    elseif dir == 1
        leg_idx = loop_idx;
    else
        leg_idx=-1;
        prev_pos = 0;
    end
    if dir==1
        direction=-1;
    elseif dir==2
        direction=1;
    end
    if leg_idx>0
        prev_pos = LEFSystem.smcs(leg_idx);
    end
    % check if the loop was attached to the chromatin
    sta = 1;
    if prev_pos > 0
        % make a step only if there is no boundary and the new position is unoccupied
        if LEFSystem.perms(prev_pos + (direction + 1) / 2) > 0
            if LEFSystem.smclattice(prev_pos+direction) < 0 && LEFSystem.geneon(loop_idx)==1
                sta = sta * LEFSystem.make_smcstep(leg_idx, direction);
                regenerate_event(LEFSystem, evheap, event_idx);
                regenerate_neighbours(LEFSystem, evheap, prev_pos-direction,dir);
            end
        end
    end
    
    %LEF initialtion gene start
elseif event_idx > LEFSystem.N && event_idx <= 2 * LEFSystem.N    
    loop_idx = event_idx - LEFSystem.N;
    sta=2;
    if LEFSystem.smc_dir(loop_idx) == 1
        smcpos = LEFSystem.smcs(loop_idx);
    elseif LEFSystem.smc_dir(loop_idx) == 2
        smcpos = LEFSystem.smcs(loop_idx + LEFSystem.N);
    else
        smcpos = -1;
    end
    
    if smcpos>0
        LEFSystem.geneon(loop_idx) = 1;
        regenerate_event(LEFSystem, evheap, loop_idx);
    end
    
    %LEF unbinding
elseif event_idx > 2 * LEFSystem.N && event_idx <= 3 * LEFSystem.N    
    loop_idx = event_idx - 2 * LEFSystem.N;
    
    sta = 2;
    % check if the loop was attached to the chromatin
    if LEFSystem.smcs(loop_idx) <= 0
        sta = 0;
    end
    
    % save previous positions, but don't update neighbours until the loop
    % has moved
    prev_pos = LEFSystem.smcs(loop_idx);
    prev_pos2 = LEFSystem.smcs(loop_idx+LEFSystem.N);
    
    sta = sta * LEFSystem.move_smc(loop_idx, -1);
    sta = sta * LEFSystem.move_smc(loop_idx+LEFSystem.N, -1);
    
    LEFSystem.smc_dir(loop_idx) = -1;
    LEFSystem.geneon(loop_idx) = -1;
    
    % regenerate events for the loop itself and for its previous neighbours
    regenerate_event(LEFSystem, evheap, loop_idx+3*LEFSystem.N);
    
    % update the neighbours after the loop has moved
    regenerate_neighbours(LEFSystem, evheap, prev_pos-1,2);
    regenerate_neighbours(LEFSystem, evheap, prev_pos+1,1);
    regenerate_neighbours(LEFSystem, evheap, prev_pos2-1,2);
    regenerate_neighbours(LEFSystem, evheap, prev_pos2+1,1);
    
    %LEF rebinding
elseif event_idx > 3 * LEFSystem.N && event_idx <= 4 * LEFSystem.N    
    loop_idx = event_idx - 3 * LEFSystem.N;
    
    sta = 2;
    % check if the loop was not attached to the chromatin
    if LEFSystem.smcs(loop_idx) > 0
        sta = 0;
    end
    
    % find a new position for the LEM (a brute force method,
    % can be improved)
    eptsr = [];
    if ismember(-1,LEFSystem.smclattice(LEFSystem.DR2(:,1)))        
        while 1
            if length(eptsr) == length(LEFSystem.DR2(:,1))
                new_pos = -1;
                break
            end
            p1 = randi([1 length(LEFSystem.DR2(:,1))],1,1);
            if ismember(p1,eptsr)
                continue
            else
                new_pos = LEFSystem.DR2(p1);
                if (ismember(new_pos,LEFSystem.FST(:,1)) && LEFSystem.smclattice(new_pos+1) <= 0 && LEFSystem.smclattice(new_pos) <= 0 ) ||...
                        (ismember(new_pos,LEFSystem.RST(:,2)) && LEFSystem.smclattice(new_pos-1) <= 0 && LEFSystem.smclattice(new_pos) <= 0 )
                    break
                end
            end
            eptsr = [eptsr,p1];
        end
        if new_pos > 0
            if (ismember(new_pos,LEFSystem.FST(:,1)) && LEFSystem.smclattice(new_pos+1) <= 0 && LEFSystem.smclattice(new_pos) <= 0 ) &&...
                    (ismember(new_pos,LEFSystem.RST(:,2)) && LEFSystem.smclattice(new_pos-1) <= 0 && LEFSystem.smclattice(new_pos) <= 0 )
                LEFSystem.smc_dir(loop_idx) = randi(2,1);
            elseif ismember(new_pos,LEFSystem.FST(:,1)) && LEFSystem.smclattice(new_pos+1) <= 0 && LEFSystem.smclattice(new_pos) <= 0
                LEFSystem.smc_dir(loop_idx) = 2;
            elseif ismember(new_pos,LEFSystem.RST(:,2)) && LEFSystem.smclattice(new_pos-1) <= 0 && LEFSystem.smclattice(new_pos) <= 0
                LEFSystem.smc_dir(loop_idx) = 1;
            end
            
            if LEFSystem.smc_dir(loop_idx) == 1
                sta = sta * LEFSystem.move_smc(loop_idx, new_pos-1);
                sta = sta * LEFSystem.move_smc(loop_idx+LEFSystem.N, new_pos);
            elseif LEFSystem.smc_dir(loop_idx) == 2
                sta = sta * LEFSystem.move_smc(loop_idx, new_pos);
                sta = sta * LEFSystem.move_smc(loop_idx+LEFSystem.N, new_pos+1);
            end
            regenerate_event(LEFSystem, evheap, loop_idx+LEFSystem.N);
            regenerate_event(LEFSystem, evheap, loop_idx+2*LEFSystem.N);
        else
            regenerate_event(LEFSystem, evheap, event_idx);
        end
    else
        regenerate_event(LEFSystem, evheap, event_idx);
    end    
    
else
    disp(strcat('event_idx assumed a forbidden value :', num2str(event_idx)))
    sta = 0;
end
end
