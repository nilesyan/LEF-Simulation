function InvCovM = Cal_InvCovM(pl,PT1,PT2)

point_sort = [pl(:,1);pl(:,2)];
point_sort = sort(point_sort);
point_sort = unique(point_sort);
pl(:,3) = -1;
pl(1,3) = 1;
pl(2,3) = 2;
sz=size(pl,1)-2;
label = 3;

PT1id=find(pl(3:end,:)==PT1);
PT2id=find(pl(3:end,:)==PT2);
if PT1id>sz
    PT1id=PT1id-sz;
end
pl(PT1id+2,3)=1;
if PT2id>sz
    PT2id=PT2id-sz;
end
pl(PT2id+2,3)=2;

for i=3:size(pl,1)
    if pl(i,3)==-1
        pl(i,3)=label;
        label=label+1;
    end
end

Di = max(pl(:,3));
CovM = zeros(Di,Di);
for i=1:length(point_sort)
    if i==1
        Pa = point_sort(i);
        id1=find(pl(3:end,1:2)==point_sort(i));
        if isempty(id1)
            if point_sort(i)==PT1
                a=1;
            elseif point_sort(i)==PT2
                a=2;
            end
        else
            if id1>sz
                id1=id1-sz;
            end
            a=pl(id1+2,3);
        end
    else
        Pb = point_sort(i);
        id2=find(pl(3:end,1:2)==point_sort(i));
        if isempty(id2)
            if point_sort(i)==PT1
                b=1;
            elseif point_sort(i)==PT2
                b=2;
            end
        else
            if id2>sz
                id2=id2-sz;
            end
            b=pl(id2+2,3);
        end
        if a<b
            if CovM(a,b)<0
                CovM(a,b) = -1/(Pb-Pa)+CovM(a,b);                
            else
                CovM(a,b) = -1/(Pb-Pa);
            end
        elseif a>b
            if CovM(b,a)<0
                CovM(b,a) = -1/(Pb-Pa)+CovM(b,a);                
            else
                CovM(b,a) = -1/(Pb-Pa);
            end
        end
        Pa = Pb;
        a = b;
    end
end

CovMa = CovM+CovM';
InvCovM = zeros(Di,Di);
for i=1:Di
    InvCovM(i,i) = -sum(CovMa(i,:));
end
InvCovM = InvCovM+CovMa;

end