function inde = sumop(A)
    t=[];
    for i=1:numel(A(:,1))
        t=[t sum(A(i,:))];
    end
    [M,ind]=min(t);
    inde=ind;
end