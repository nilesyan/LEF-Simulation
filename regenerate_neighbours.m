function regenerate_neighbours(LEFSystem, evheap, pos, dir)

%     Regenerate the motion events for the adjacent loop legs.
%     Use to unblock the previous neighbors and block the new ones.

if pos > 1 && pos < LEFSystem.L
    smcid = LEFSystem.smclattice(pos);
    if smcid>0
        if smcid>LEFSystem.N
            smcid = smcid-LEFSystem.N;
        end
        smc_dir = LEFSystem.smc_dir(smcid);
        % If the SMC direction is the same as the direction of the empty lattice,
        % regenerate SMC translocation event.
        if smc_dir == dir 
            if LEFSystem.geneon(smcid)==1
                regenerate_event(LEFSystem, evheap, smcid)
            end
        end
    end
end
end