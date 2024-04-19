
clear
nPtsPerClust = 50;
nClust  = 2;
totalNumPts = nPtsPerClust*nClust;
m(:,1) = [1 1]';
m(:,2) = [-1 -1]';
m(:,3) = [1 -1]';
var = .6;
bandwidth = 35;
clustMed = [];
%clustCent;
% x = var*randn(2,nPtsPerClust*nClust);
% %*** build the point set
% for i = 1:nClust
%     x(:,1+(i-1)*nPtsPerClust:(i)*nPtsPerClust)       = x(:,1+(i-1)*nPtsPerClust:(i)*nPtsPerClust) + repmat(m(:,i),1,nPtsPerClust);   
% end

num_ue = 100;
uex = 100*rand(num_ue,1)+ 24;
uey = 100*rand(num_ue,1) + 12;
x = [uex,uey]';

instbatt = 4000*rand(1,num_ue) + 500;
decay =   randi(15,1,num_ue)+5;
disratee = instbatt./decay;
X = [uex,uey];
[tp4ms, pc4ms, sizee] = mspc(X,disratee,num_ue);
if(sizee==3)
    [tp4knn,pc4knn] = clustere3(X,disratee,num_ue);
elseif(sizee==4)
    [tp4knn,pc4knn] = clustere(X,disratee,num_ue);
else
    [tp4knn,pc4knn] = clustere5(X,disratee,num_ue);
end

figure

barX = categorical({'KNN','Mean shift'});
barX= reordercats(barX,{'KNN','Mean shift'});
bar(barX, [tp4ms(2) tp4knn(2)],0.4);
ylabel('Sum rate in Bps');
title('kNN vs Meanshift : Sumrate');
figure
barX = categorical({'KNN','Mean shift'});
barX= reordercats(barX,{'KNN','Mean shift'});
bar(barX, [pc4ms pc4knn],0.4);
ylabel('Power consumption in dBm');
title('kNN vs Meanshift : Power Consumption');
% tic
% [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(x,bandwidth);
% toc
% numClust = length(clustMembsCell);
% figure(10),clf,hold on
% cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
% for k = 1:min(numClust,length(cVec))
%     myMembers = clustMembsCell{k};
%     myClustCen = clustCent(:,k);
%     plot(x(1,myMembers),x(2,myMembers),[cVec(k) '.'])
%     plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
% end
% title(['no shifting, numClust:' int2str(numClust)])