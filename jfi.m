function qos=jfi(X,ind)
    x=[];
    for i = 1:numel(X(1,:))
        x=[x X(ind,i)];
    end
    x=x./mean(x);
%     nr=sum(x)*sum(x);
%     x2=x.*x;
%     dr=numel(X(1,:))*sum(x2);
    qos=1/(1+var(x)/mean(x)); %Coeffecient of variation = Variance/Mean
end