function [ra,rb,rc,rd]= randselect(A,B,C,D)
    ra=randi(numel(A(:,1)));
    rb=randi(numel(B(:,1)));
    rc=randi(numel(C(:,1)));
    rd=randi(numel(D(:,1)));
end