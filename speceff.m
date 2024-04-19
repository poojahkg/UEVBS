function se=speceff(p1,p2,X)
    hnm=randn(1)+j*randn(1);
    x=[];
    for i = 1:numel(X(:,1))
        x=[x X(i,4)];
    end
    sx=sum(x);
    sx=(sx-p2)*hnm;
    c2= abs((p1*hnm)/(p2*(randn(1)+j*randn(1))+sx+1));
    se=log2(1+c2);
end