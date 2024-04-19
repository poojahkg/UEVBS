%Dynamic Allocation of UE VBS in 5G networks

%The explanation for the many blocks in the code is given directly UNDER
% each block
num_ue = 100;
uex = 100*rand(num_ue,1)+ 24;
uey = 100*rand(num_ue,1) + 12;
X = [uex,uey];
%Initializing the position of UEs as Cartesian pairs,i.e, (x,y) coordinates
figure;
plot(X(:,1),X(:,2),'.','MarkerSize',15);
%Variable X is a matrix containing the details of all UEs. The same name
%shall be used in all the user defined functions to access the position,
%received power from BS and dis rate of any 'i'th UE
hold on;
dt=20;
speed=randi([-10,10],num_ue,1);
angle=randi([10 350],num_ue,1);
%Initializing the angle,speed and time interval in order to faciliate movement
%of UEs. This is the part where we make the selection dynamic.
uex = uex+speed*dt;
uey = uey+cos(angle*pi/180)*dt;
X1 = [uex,uey];
plot(X1(:,1),X1(:,2),'*','MarkerSize',5);

uex = uex+speed*dt;
uey = uey+cos(angle*pi/180)*dt;
X2 = [uex,uey];
plot(X2(:,1),X2(:,2),'o','MarkerSize',5);
hold on;

instbatt = 4000*rand(1,num_ue) + 500;
decay =   randi(15,1,num_ue)+5;
disratee = instbatt./decay;
%Intialising the battery discharge rate at random

[tp2,pc2] =clustere2(X,disratee,num_ue);
% 
[tp3,pc3] = clustere3(X,disratee,num_ue);
% %clustere3(X1,disratee);
% %clustere3(X2,disratee);
% 
[tp4,pc4] = clustere(X,disratee,num_ue);
%clustere(X1,disratee,num_ue);
%clustere(X2,disratee,num_ue);
% 
[tp5,pc5]=clustere5(X,disratee,num_ue);
% %clustere5(X1,disratee);
% %clustere5(X2,disratee);
% 
[tp6,pc6]=clustere6(X,disratee,num_ue);
% 
[tp7,pc7]=clustere7(X,disratee,num_ue);
% 
[tp8,pc8]=clustere8(X,disratee,num_ue);

%Each cluster function returns the average throughput found in the
%particular cluster size. The explanation for the said functions are given
%inside the function blocks

figure;
barX = categorical({'Cluster size 2','Cluster size 3','Cluster size 4','Cluster size 5','Cluster size 6','Cluster size 7','Cluster size 8'});
barX= reordercats(barX,{'Cluster size 2','Cluster size 3','Cluster size 4','Cluster size 5','Cluster size 6','Cluster size 7','Cluster size 8'});
bar(barX,[tp2(2) tp3(2) tp4(2) tp5(2) tp6(2) tp7(2) tp8(2)],0.4);
ylabel('Sum rate in bps');
title('Average Sum rate for varying cluster sizes');

%Semi log graph is used to plot the throughput of 4G and 'No UE VBS'
%variants along with the rest of the cluster sizes as the latter has TP in
%the order of Mbps, whereas UE VBS networks have TP in the order of Gbps

figure;
barX = categorical({'Cluster size 2','Cluster size 3','Cluster size 4','Cluster size 5','Cluster size 6','Cluster size 7','Cluster size 8'});
barX= reordercats(barX,{'Cluster size 2','Cluster size 3','Cluster size 4','Cluster size 5','Cluster size 6','Cluster size 7','Cluster size 8'});
bar(barX,[pc2 pc3 pc4 pc5 pc6 pc7 pc8],0.4);
ylabel('Power Consumption in dBm');
title('Power consumption for varying cluster sizes');
