function [m,to]=sumse(A)
    t=[];
    for i=1:numel(A(:,1))
        t=[t sum(A(i,:))];
    end
    [m,to]=max(t);
end