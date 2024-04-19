num_ue = 200;
uex = 100*rand(num_ue,1)+ 24;
uey = 100*rand(num_ue,1) + 12;
X = [uex,uey];
nPtsPerClust = 10;
nClust  = 2;
totalNumPts = nPtsPerClust*nClust;
m(:,1) = [1 1]';
m(:,2) = [-1 -1]';
m(:,3) = [1 -1]';
var = .6;
bandwidth = 35;
clustMed = [];
x = [uex,uey]';
[clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(x,bandwidth);

%size = max(point2cluster);
size=4;
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
if (size==2)
    [tp2,pc2] =clustere2(X,disratee,num_ue);
    [tp21,pc21] =clustere2(X1,disratee,num_ue);
    [tp22,pc22] =clustere2(X2,disratee,num_ue);
% 
elseif (size==3)
    [tp3,pc3] = clustere3(X,disratee,num_ue);
    [tp31,pc31] =clustere3(X1,disratee,num_ue);
    [tp32,pc32] =clustere3(X2,disratee,num_ue);
% 
elseif(size==4)
    [tp4,pc4] = clustere(X,disratee,num_ue);
    [tp41,pc41] =clustere(X1,disratee,num_ue);
    [tp42,pc42] =clustere(X2,disratee,num_ue);
% 
elseif(size==5)
    [tp5,pc5]=clustere5(X,disratee,num_ue);
    [tp51,pc51] =clustere5(X1,disratee,num_ue);
    [tp52,pc52] =clustere5(X2,disratee,num_ue);
% 
elseif(size==6)
    [tp6,pc6]=clustere6(X,disratee,num_ue);
    [tp61,pc61] =clustere6(X1,disratee,num_ue);
    [tp62,pc62] =clustere6(X2,disratee,num_ue);
% 
elseif(size==7)
    [tp7,pc7]=clustere7(X,disratee,num_ue);
    [tp71,pc71] =clustere7(X1,disratee,num_ue);
    [tp72,pc72] =clustere7(X2,disratee,num_ue);
% 
else
    [tp8,pc8]=clustere8(X,disratee,num_ue);
    [tp81,pc81] =clustere8(X1,disratee,num_ue);
    [tp82,pc82] =clustere8(X2,disratee,num_ue);
end
%Each cluster function returns the average throughput found in the
%particular cluster size. The explanation for the said functions are given
%inside the function blocks
