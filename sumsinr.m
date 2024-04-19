function [maxx,tow] = sumsinr(A)
    t=[];
    for i=1:numel(A(:,1))
        t=[t sum(A(i,:))];
    end
    [maxx,tow]=max(t);
end