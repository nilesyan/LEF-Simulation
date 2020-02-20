function Cor_NP =  N_eff_TwoPoints(l_foot,r_foot,PT1,PT2,N)
% Calculate the effective genomic distance between two points

if PT1<0 || PT2 <0
    Cor_NP = -1e6;
    return
end

for i=1:N
    if l_foot(i) > r_foot(i)
        tp = l_foot(i);
        l_foot(i) = r_foot(i);
        r_foot(i) = tp;
    end
end

pl = FindPairs(l_foot,r_foot,PT1,PT2);
InvCovM = Cal_InvCovM(pl,PT1,PT2);
U0 = InvCovM(1:2,1:2);
W0 = InvCovM(3:end,3:end);
V = InvCovM(1:2,3:end);

U = U0-V*inv(W0)*V';

Cor_NP = 1/U(1,1);

end