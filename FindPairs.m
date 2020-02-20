function pl=FindPairs(l_foot,r_foot,PT1,PT2)

lf=l_foot;
rf=r_foot;
pl=[PT1,PT1;PT2,PT2];
sz=0;
while size(pl,1)~=sz
    Ps = min(min(pl));
    Pl = max(max(pl));
    sz = size(pl,1);
    k=[];
    for i=1:length(lf)
        if (rf(i)>=Ps && lf(i)<=Pl) || (rf(i)>=Pl && lf(i)<=Ps)
            pl=[pl;[lf(i),rf(i)]];
            k=[k,i];
        end
    end
    lf(k)=[];
    rf(k)=[];
end

end
