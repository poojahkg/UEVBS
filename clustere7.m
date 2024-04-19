%The algorithm implemented in all the cluster functions are the same with
%the only difference being the cluster size. Hence, the explanation will be
%given in this function block alone.
function [tp7, powcons7] = cluster7(X,disratee,num_ue)
    disrate=disratee;
    figure;
    alpha=0.5;
    beta=0.5;
    gnbc = [80,60]; %gNB - Base Station
    d=pdist2(gnbc,X,'euclidean')*20;
    clusterX = kmeans(X,7);
    plot(X(clusterX==1,1),X(clusterX==1,2),'r.','MarkerSize',15)
    hold on
    plot(X(clusterX==2,1),X(clusterX==2,2),'g.','MarkerSize',15);
    plot(X(clusterX==3,1),X(clusterX==3,2),'c.','MarkerSize',15);
    plot(X(clusterX==4,1),X(clusterX==4,2),'b.','MarkerSize',15);
    plot(X(clusterX==5,1),X(clusterX==5,2),'m.','MarkerSize',15);
    plot(X(clusterX==6,1),X(clusterX==6,2),'y.','MarkerSize',15);
    plot(X(clusterX==7,1),X(clusterX==7,2),'k.','MarkerSize',15);
    drawnow;
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Clustering using K-Means algorithm');
    %Plotting the clustered points
    DistanceMsr=60;
    LightSpeedC=3e8; 
    Freq=5e9; %FR1
    Freq4g=1.9e9;
    TXAntennaGain=1;%db
    RXAntennaGain=1;%d
    PTx=0.25;%watt
    sigma=6;%Sigma from 6 to 12 %Principles of communication systems simulation with wireless application P.548
    PathLossExponent=4;%Line Of sight
    Wavelength=LightSpeedC/Freq;
    Wavelength4g=LightSpeedC/Freq4g;
    %PTxdBm=10*log10(PTx*1000);
    PTxdBm = 23.98;
    Dref=10;
    %The above parameters are used in the log-normal path loss model
    cone = [];
    ctwo = [];
    cthree = [];
    cfour = [];
    cfive = [];
    csix = [];
    cseven = [];
    %The arrays defined above contain the UE details for each cluster
    
    pr1w =[];
    pr2w = [];
    pr3w = [];
    pr4w = [];
    pr5w = [];
    pr6w = [];
    pr7w = [];
    
    %Arrays for storing received power of each UE from BBS in watt (for
    %evey cluster)
    
    Pr4g=[];
    rstate = randn('state');
    GaussRandom= (randn*0.1+0);
    %Power consumed in each cluster = (No.of UE in each cluster)*(BS power in watts) 
    % - (Pr of each UE from BS) + (No UEs-1)*(Pr of UE VBS) - (Pr of each UE from UE VBS)
    for i = 1:num_ue
        M = Wavelength4g / (4 * pi * d(i));
        Pr0=43 + TXAntennaGain + RXAntennaGain- (20*log10(1/M)); %43dBm is the standard BS power
        Pr4g=[Pr4g Pr0+(10*PathLossExponent* log10(d(i)/Dref))+GaussRandom];
        randn('state', rstate);
        if(clusterX(i)==1)
            M = Wavelength / (4 * pi * d(i));
            Pr0=PTxdBm + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            Pr1=Pr0+(10*PathLossExponent* log10(d(i)/Dref))+GaussRandom; %Pr equation
            pr1w = [pr1w 10^((-1*Pr1-30)/10)];
            randn('state', rstate);
            X(i,3) = disrate(1,i);
            X(i,4) = Pr1;
            X(i,5)=alpha*X(i,3)+beta*X(i,4); %The parameter combining Pr and dis. rate is included in the last column - master parameter
            X(i,6) = d(i);
            cone = [cone;X(i,:)]; %The details get appended wrt cluster number. Same procedure is followed in the entire loop for other clusters.
        elseif (clusterX(i)==2)   
            M = Wavelength / (4 * pi * d(i));
            Pr0=PTxdBm + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            Pr2=Pr0+(10*PathLossExponent* log10(d(i)/Dref))+GaussRandom;
            pr2w = [pr2w 10^((-1*Pr2-30)/10)];
            randn('state', rstate); 
            X(i,3) = disrate(1,i);
            X(i,4) = Pr2;
            X(i,5)=alpha*X(i,3)+beta*X(i,4);
            X(i,6) = d(i);
            ctwo = [ctwo;X(i,:)];
        elseif (clusterX(i)==3)
            M = Wavelength / (4 * pi * d(i));
            Pr0=PTxdBm + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            Pr3=Pr0+(10*PathLossExponent* log10(d(i)/Dref))+GaussRandom;
            pr3w = [pr3w 10^((-1*Pr3-30)/10)];
            randn('state', rstate); 
            X(i,3) = disrate(1,i);
            X(i,4) = Pr3; 
            X(i,5)=alpha*X(i,3)+beta*X(i,4);
            X(i,6) = d(i);
            cthree = [cthree;X(i,:)]; 
        elseif (clusterX(i)==4)
            M = Wavelength / (4 * pi * d(i));
            Pr0=PTxdBm + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            Pr4=Pr0+(10*PathLossExponent* log10(d(i)/Dref))+GaussRandom;
            pr4w = [pr4w 10^((-1*Pr4-30)/10)];
            randn('state', rstate); 
            X(i,3) = disrate(1,i);
            X(i,4) = Pr4;
            X(i,5)=alpha*X(i,3)+beta*X(i,4);
            X(i,6) = d(i);
            cfour = [cfour;X(i,:)]; 
        elseif (clusterX(i)==5)
            M = Wavelength / (4 * pi * d(i));
            Pr0=PTxdBm + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            Pr5=Pr0+(10*PathLossExponent* log10(d(i)/Dref))+GaussRandom;
            pr5w = [pr5w 10^((-1*Pr5-30)/10)];
            randn('state', rstate); 
            X(i,3) = disrate(1,i);
            X(i,4) = Pr5;
            X(i,5)=alpha*X(i,3)+beta*X(i,4);
            X(i,6) = d(i);
            cfive = [cfive;X(i,:)]; 
        elseif (clusterX(i)==6)
            M = Wavelength / (4 * pi * d(i));
            Pr0=PTxdBm + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            Pr6=Pr0+(10*PathLossExponent* log10(d(i)/Dref))+GaussRandom;
            pr6w = [pr6w 10^((-1*Pr6-30)/10)];
            randn('state', rstate); 
            X(i,3) = disrate(1,i);
            X(i,4) = Pr6;
            X(i,5)=alpha*X(i,3)+beta*X(i,4);
            X(i,6) = d(i);
            csix = [csix;X(i,:)]; 
        else
            M = Wavelength / (4 * pi * d(i));
            Pr0=PTxdBm + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            Pr7=Pr0+(10*PathLossExponent* log10(d(i)/Dref))+GaussRandom;
            pr7w = [pr7w 10^((-1*Pr7-30)/10)];
            randn('state', rstate); 
            X(i,3) = disrate(1,i);
            X(i,4) = Pr7;
            X(i,5)=alpha*X(i,3)+beta*X(i,4);
            X(i,6) = d(i);
            cseven = [cseven;X(i,:)]; 
        end
    end
    %Power consumed in each cluster = (No.of UE in each cluster)*(BS power in watts) 
    % - (Pr of each UE from BS) + (No UEs-1)*(Pr of UE VBS) - (Pr of each UE from UE VBS)
    pc1h = numel(cone)*PTx - sum(pr1w);
    pc2h = numel(ctwo)*PTx - sum(pr2w);
    pc3h = numel(cthree)*PTx - sum(pr3w);
    pc4h = numel(cfour)*PTx - sum(pr4w);
    pc5h = numel(cfive)*PTx - sum(pr5w);
    pc6h = numel(csix)*PTx - sum(pr6w);
    pc7h = numel(cseven)*PTx - sum(pr7w);
   
    hold on;
    e1=[];
    e2=[];
    e3=[];
    e4=[];
    e5=[];
    e6=[];
    e7=[];
    %Arrays ei are the ones holding the indices (in X) of the eligible UE VBSs. i
    %ranges from 1 to cluster-size
    ss1=[];
    ss2=[];
    ss3=[];
    ss4=[];
    ss5=[];
    ss6=[];
    ss7=[];
    
    sp1=[];
    sp2=[];
    sp3=[];
    sp4=[];
    sp5=[];
    sp6=[];
    sp7=[];
    
    %Arrays ssi are the ones holding the SINR of all the eligible UE VBSs in each cluster
    ra=randi(numel(cone(:,1)));
    rb=randi(numel(ctwo(:,1)));
    rc=randi(numel(cthree(:,1)));
    rd=randi(numel(cfour(:,1)));
    re=randi(numel(cfive(:,1)));
    rf=randi(numel(csix(:,1)));
    rg=randi(numel(cseven(:,1)));
    
    %Selecting UE VBS at random for comparison
    sr1=[];
    sr2=[];
    sr3=[];
    sr4=[];
    sr5=[];
    sr6=[];
    sr7=[];
    %sri holds the SINR of random UE VBSs
    a1=cone(:,5);
    %ai holds the master parameter of each UE. It's a column vector which
    %can be used to find the size of the cluster as well,i.e, numel(a1).
    [M,I1]=max(a1);
    k=1;
    for i=1:numel(cone(:,1))
        if(a1(i)>mean(a1)) %UEs with master parameters more than the mean values are selected as eligible UE VBSs (eUE)
            k=i;
            e1=[e1 k];
        end
    end
    for i=1:numel(e1)
        for j=1:numel(cone(:,1))
            if(j~=e1(i))
                diss= sqrt((cone(e1(i),1)-cone(j,1))^2-(cone(e1(i),2)-cone(j,2))^2);
                M = Wavelength / (4 * pi * diss);
                Pr0=X(e1(i),4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
                PrX=Pr0+(10*PathLossExponent* log10(cone(e1(i),6)/Dref))+GaussRandom; %Measuring the power received by an UE from an eligible UE VBS
                %s1=sinnr(e1(i),cone(e1(i),1),cone(e1(i),2),cone(j,1),cone(j,2),cone,cone(e1(i),3));
                s1=sinnr(e1(i),cone(e1(i),1),cone(e1(i),2),cone(j,1),cone(j,2),cone,PrX); %Sending index of eUE, its position, position of receiving UE and received power
                ss1=[ss1 s1]; %Compiling the SINR values
                se1=speceff(PrX,cone(e1(i),4),cone);
                sp1=[sp1 se1];
            else
                ss1=[ss1 0]; 
                sp1=[sp1 0];
            end
        end
    %ee1=[ee1;e1];
    end
    for j=1:numel(cone(:,1))
        if(j~=ra)
            diss= sqrt((cone(ra,1)-cone(j,1))^2-(cone(ra,2)-cone(j,2))^2);
            M = Wavelength / (4 * pi * diss);
            Pr0=X(ra,4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            PrX=Pr0+(10*PathLossExponent* log10(cone(ra,6)/Dref))+GaussRandom;
            s1=sinnr(ra,cone(ra,1),cone(ra,2),cone(j,1),cone(j,2),cone,PrX);
            sr1=[sr1 s1];
            %Same procedure as above but for a random UE
        else
            sr1=[sr1 0];
        end
    end
    
    %The above procedure is repeated in each cluster as shown below
    
    a2=ctwo(:,5);
    [M,I2]=max(a2);

    k=1;
    for i=1:numel(ctwo(:,1))
        if(a2(i)>mean(a2))
            k=i;
            e2=[e2 k];
        end
    end
    for i=1:numel(e2)
        for j=1:numel(ctwo(:,1))
            if(j~=e2(i))
                diss= sqrt((ctwo(e2(i),1)-ctwo(j,1))^2-(ctwo(e2(i),2)-ctwo(j,2))^2);
                M = Wavelength / (4 * pi * diss);
                Pr0=X(e2(i),4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
                PrX=Pr0+(10*PathLossExponent* log10(ctwo(e2(i),6)/Dref))+GaussRandom;
                %s2=sinnr(e2(i),ctwo(e2(i),1),ctwo(e2(i),2),ctwo(j,1),ctwo(j,2),ctwo,ctwo(e2(i),3));
                s2=sinnr(e2(i),ctwo(e2(i),1),ctwo(e2(i),2),ctwo(j,1),ctwo(j,2),ctwo,PrX);
                ss2=[ss2 s2];
                se2=speceff(PrX,ctwo(e2(i),4),ctwo);
                sp2=[sp2 se2];
            else
                ss2=[ss2 0]; 
                sp2=[sp2 0];
            end
        end
    %ee2=[ee2;e2];
    end       
    for j=1:numel(ctwo(:,1))
        if(j~=rb)
            diss= sqrt((ctwo(rb,1)-ctwo(j,1))^2-(ctwo(rb,2)-ctwo(j,2))^2);
            M = Wavelength / (4 * pi * diss);
            Pr0=X(rb,4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            PrX=Pr0+(10*PathLossExponent* log10(ctwo(rb,6)/Dref))+GaussRandom;
            s2=sinnr(rb,ctwo(rb,1),ctwo(rb,2),ctwo(j,1),ctwo(j,2),ctwo,PrX);
            sr2=[sr2 s2];
        else
            sr2=[sr2 0];
        end
    end
    a3=cthree(:,5);
    [M,I3]=max(a3);

    k=1;
    for i=1:numel(cthree(:,1))
        if(a3(i)>mean(a3))
            k=i;
            e3=[e3 k];
        end
    end
    for i=1:numel(e3)
        for j=1:numel(cthree(:,1))
            if(j~=e3(i))
                diss= sqrt((cthree(e3(i),1)-cthree(j,1))^2-(cthree(e3(i),2)-cthree(j,2))^2);
                M = Wavelength / (4 * pi * diss);
                Pr0=X(e3(i),4) + TXAntennaGain + RXAntennaGain - (20*log10(1/M));
                PrX=Pr0+(10*PathLossExponent* log10(cthree(e3(i),6)/Dref))+GaussRandom;
                %s3=sinnr(e3(i),cthree(e3(i),1),cthree(e3(i),2),cthree(j,1),cthree(j,2),cthree,cthree(e3(i),3));
                s3=sinnr(e3(i),cthree(e3(i),1),cthree(e3(i),2),cthree(j,1),cthree(j,2),cthree,PrX);
                ss3=[ss3 s3];
                se3=speceff(PrX,cthree(e3(i),4),cthree);
                sp3=[sp3 se3];
            else
                ss3=[ss3 0]; 
                sp3=[sp3 0];
            end
        end
    %ee3=[ee3;e3];
    end    
    for j=1:numel(cthree(:,1))
        if(j~=rc)
            diss= sqrt((cthree(rc,1)-cthree(j,1))^2-(cthree(rc,2)-cthree(j,2))^2);
            M = Wavelength / (4 * pi * diss);
            Pr0=X(rc,4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            PrX=Pr0+(10*PathLossExponent* log10(cthree(rc,6)/Dref))+GaussRandom;
            s3=sinnr(rc,cthree(rc,1),cthree(rc,2),cthree(j,1),cthree(j,2),cthree,PrX);
            sr3=[sr3 s3];
        else
            sr3=[sr3 0];
        end
    end
    a4=cfour(:,5);
    [M,I4]=max(a4);

    k=1;

    for i=1:numel(cfour(:,1))
        if(a4(i)>mean(a4))
            k=i;
            e4=[e4 k];
        end
    end
    for i=1:numel(e4)
        for j=1:numel(cfour(:,1))
            if(j~=e4(i))
                diss= sqrt((cfour(e4(i),1)-cfour(j,1))^2-(cfour(e4(i),2)-cfour(j,2))^2);
                M = Wavelength / (4 * pi * diss);
                Pr0=X(e4(i),4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
                PrX=Pr0+(10*PathLossExponent* log10(cfour(e4(i),6)/Dref))+GaussRandom;
                %s4=sinnr(e4(i),cfour(e4(i),1),cfour(e4(i),2),cfour(j,1),cfour(j,2),cfour,cfour(e4(i),3));
                s4=sinnr(e4(i),cfour(e4(i),1),cfour(e4(i),2),cfour(j,1),cfour(j,2),cfour,PrX);
                ss4=[ss4 s4];
                se4=speceff(PrX,cfour(e4(i),4),cfour);
                sp4=[sp4 se4];
            else
                ss4=[ss4 0]; 
                sp4=[sp4 0];
            end
        end
    %ee4=[ee4;e4];
    end

    for j=1:numel(cfour(:,1))
        if(j~=rd)
            diss= sqrt((cfour(rd,1)-cfour(j,1))^2-(cfour(rd,2)-cfour(j,2))^2);
            M = Wavelength / (4 * pi * diss);
            Pr0=X(rd,4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            PrX=Pr0+(10*PathLossExponent* log10(cfour(rd,6)/Dref))+GaussRandom;
            s4=sinnr(rd,cfour(rd,1),cfour(rd,2),cfour(j,1),cfour(j,2),cfour,PrX);
            sr4=[sr4 s4];
        else
            sr4=[sr4 0];
        end
    end
    
     a5=cfive(:,5);
   [M,I5]=max(a5);
   k=1;
   for i=1:numel(cfive(:,1))
       if(a5(i)>mean(a5))
           k=i;
           e5=[e5 k];
       end
   end
   for i=1:numel(e5)
       for j=1:numel(cfive(:,1))
           if(j~=e5(i))
               diss= sqrt((cfive(e5(i),1)-cfive(j,1))^2-(cfive(e5(i),2)-cfive(j,2))^2);
               M = Wavelength / (4 * pi * diss);
               Pr0=X(e5(i),4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
               PrX=Pr0+(10*PathLossExponent* log10(cfive(i,6)/Dref))+GaussRandom;
               s5=sinnr(e5(i),cfive(e5(i),1),cfive(e5(i),2),cfive(j,1),cfive(j,2),cfive,PrX);
               ss5=[ss5 s5];
               se5=speceff(PrX,cfive(e5(i),4),cfive);
               sp5=[sp5 se5];
            else
                ss5=[ss5 0]; 
                sp5=[sp5 0];
           end
       end
   %ee1=[ee1;e1];
   end
   for j=1:numel(cfive(:,1))
       if(j~=re)
           diss= sqrt((cfive(re,1)-cfive(j,1))^2-(cfive(re,2)-cfive(j,2))^2);
           M = Wavelength / (4 * pi * diss);
           Pr0=X(re,3) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
           PrX=Pr0-(10*PathLossExponent* log10(1/Dref))+GaussRandom;
           s5=sinnr(re,cfive(re,1),cfive(re,2),cfive(j,1),cfive(j,2),cfive,PrX);
           sr5=[sr5 s5];
       else
           sr5=[sr5 0];
       end
   end
   
    a6=csix(:,5);
   [M,I6]=max(a6);
   k=1;
   for i=1:numel(csix(:,1))
       if(a6(i)>mean(a6))
           k=i;
           e6=[e6 k];
       end
   end
   for i=1:numel(e6)
       for j=1:numel(csix(:,1))
           if(j~=e6(i))
               diss= sqrt((csix(e6(i),1)-csix(j,1))^2-(csix(e6(i),2)-csix(j,2))^2);
               M = Wavelength / (4 * pi * diss);
               Pr0=X(e6(i),4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
               PrX=Pr0+(10*PathLossExponent* log10(csix(i,6)/Dref))+GaussRandom;
               s6=sinnr(e6(i),csix(e6(i),1),csix(e6(i),2),csix(j,1),csix(j,2),csix,PrX);
               ss6=[ss6 s6];
               se6=speceff(PrX,csix(e6(i),4),csix);
               sp6=[sp6 se6];
            else
                ss6=[ss6 0]; 
                sp6=[sp6 0];
           end
       end
   %ee1=[ee1;e1];
   end
   for j=1:numel(csix(:,1))
       if(j~=rf)
           diss= sqrt((csix(rf,1)-csix(j,1))^2-(csix(rf,2)-csix(j,2))^2);
           M = Wavelength / (4 * pi * diss);
           Pr0=X(rf,4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
           PrX=Pr0+(10*PathLossExponent* log10(csix(i,6)/Dref))+GaussRandom;
           s6=sinnr(rf,csix(rf,1),csix(rf,2),csix(j,1),csix(j,2),csix,PrX);
           sr6=[sr6 s6];
       else
           sr6=[sr6 0];
       end
   end
   
   a7=cseven(:,5);
   [M,I7]=max(a7);
   k=1;
   for i=1:numel(cseven(:,1))
       if(a7(i)>mean(a7))
           k=i;
           e7=[e7 k];
       end
   end
   for i=1:numel(e7)
       for j=1:numel(cseven(:,1))
           if(j~=e7(i))
               diss= sqrt((cseven(e7(i),1)-cseven(j,1))^2-(cseven(e7(i),2)-cseven(j,2))^2);
               M = Wavelength / (4 * pi * diss);
               Pr0=X(e7(i),4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
               PrX=Pr0+(10*PathLossExponent* log10(cseven(i,6)/Dref))+GaussRandom;
               s7=sinnr(e7(i),cseven(e7(i),1),cseven(e7(i),2),cseven(j,1),cseven(j,2),cseven,PrX);
               ss7=[ss7 s7];
               se7=speceff(PrX,cseven(e7(i),4),cseven);
               sp7=[sp7 se7];
            else
                ss7=[ss7 0]; 
                sp7=[sp7 0];
           end
       end
   %ee1=[ee1;e1];
   end
   for j=1:numel(cseven(:,1))
       if(j~=rg)
           diss= sqrt((cseven(rg,1)-cseven(j,1))^2-(cseven(rg,2)-cseven(j,2))^2);
           M = Wavelength / (4 * pi * diss);
           Pr0=X(rg,4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
           PrX=Pr0+(10*PathLossExponent* log10(cseven(i,6)/Dref))+GaussRandom;
           s7=sinnr(rg,cseven(rg,1),cseven(rg,2),cseven(j,1),cseven(j,2),cseven,PrX);
           sr7=[sr7 s7];
       else
           sr7=[sr7 0];
       end
   end
    ss1=reshape(ss1,numel(e1),numel(a1));
    ss2=reshape(ss2,numel(e2),numel(a2));
    ss3=reshape(ss3,numel(e3),numel(a3));
    ss4=reshape(ss4,numel(e4),numel(a4));
    ss5=reshape(ss5,numel(e5),numel(a5));
    ss6=reshape(ss6,numel(e6),numel(a6));
    ss7=reshape(ss7,numel(e7),numel(a7));
    %ssi is converted into a 2D matrix with each row indicating eUE and
    %columns with SINR values from a particular eUE to a UE in the same
    %cluster.

    [maxx1,tow1]=sumsinr(ss1);
    [maxx2,tow2]=sumsinr(ss2);
    [maxx3,tow3]=sumsinr(ss3);
    [maxx4,tow4]=sumsinr(ss4);
    [maxx5,tow5]=sumsinr(ss5);
    [maxx6,tow6]=sumsinr(ss6);
    [maxx7,tow7]=sumsinr(ss7);
    
    %The maximum mean SINR the eUEs is selected. maxx1 contains SINR value
    %and tow1 contains the index (wrt e1) of eUE which is selected
    
    sp1=reshape(sp1,numel(e1),numel(a1));
    sp2=reshape(sp2,numel(e2),numel(a2));
    sp3=reshape(sp3,numel(e3),numel(a3));
    sp4=reshape(sp4,numel(e4),numel(a4));
    sp5=reshape(sp5,numel(e5),numel(a5));
    sp6=reshape(sp6,numel(e6),numel(a6));
    sp7=reshape(sp7,numel(e7),numel(a7));

    [m1,t1]=sumse(sp1);
    [m2,t2]=sumse(sp2);
    [m3,t3]=sumse(sp3);
    [m4,t4]=sumse(sp4);
    [m5,t5]=sumse(sp5);
    [m6,t6]=sumse(sp6);
    [m7,t7]=sumse(sp7);

    q1=jfi(sp1,t1);
    q2=jfi(sp2,t2);
    q3=jfi(sp3,t3);
    q4=jfi(sp4,t4);
    q5=jfi(sp5,t5);
    q6=jfi(sp6,t6);
    q7=jfi(sp7,t7);
    
    %These are the QoS values of UE VBS in every cluster.
    
    [rt1,tt]=sumsinr(sr1);
    [rt2,tt]=sumsinr(sr2);
    [rt3,tt]=sumsinr(sr3);
    [rt4,tt]=sumsinr(sr4);
    [rt5,tt]=sumsinr(sr5);
    [rt6,tt]=sumsinr(sr6);
    [rt7,tt]=sumsinr(sr7);
    
    %Same procedure for random UE
    
    thwoa=[];
    thwoa=[throughput(rt1,100e7) throughput(rt2,100e7) throughput(rt3,100e7) throughput(rt4,100e7) throughput(rt5,100e7) throughput(rt6,100e7) throughput(rt7,100e7)];
    %Throughput for random selection
    
    %op1=[];
    %op2=[];
    %op3=[];
    %op4=[];
    %for i=1:numel(e1)
     %   opp=outageprob(cone,e1(i),tow1);
      %  op1=[op1 opp];
    %end
    %for i=1:numel(e4)
     %   opp=outageprob(cfour,e4(i),tow4);
      %  op4=[op4 opp];
    %end
    inde1=e1(tow1);
    inde2=e2(tow2);
    inde3=e3(tow3);
    inde4=e4(tow4);
    inde5=e5(tow5);
    inde6=e6(tow6);
    inde7=e7(tow7);
    %inde contains the index (in X) of the selected UE VBS.
    ss4g=[];
    
    for i=1:num_ue
        s14=sinr4g(X(i,1),X(i,2),X,0.25,Pr4g,d(i),num_ue);
        ss4g=[ss4g s14];
    end
    %ss4g is used to contain the SINR value for 4G network
    disp(ss4g);
    maxx4g=mean(ss4g);
    
    th4temp=[throughput(maxx4g,20e7) throughput(maxx4g,20e7) throughput(maxx4g,20e7) throughput(maxx4g,20e7) throughput(maxx4g,20e7) throughput(maxx4g,20e7) throughput(maxx4g,20e7)]; 
    
    %throughput of 4G network
    plot(cone(I1,1),cone(I1,2),'r*','MarkerSize',15,'DisplayName','UEVBS1');
    plot(cone(inde1,1),cone(inde1,2),'kd','MarkerSize',15,'DisplayName','UEVBS1 with SINR');
    plot(ctwo(I2,1),ctwo(I2,2),'g*','MarkerSize',15,'DisplayName','UEVBS2');
    plot(ctwo(inde2,1),ctwo(inde2,2),'kd','MarkerSize',15,'DisplayName','UEVBS2 with SINR');
    plot(cthree(I3,1),cthree(I3,2),'c*','MarkerSize',15,'DisplayName','UEVBS3');
    plot(cthree(inde3,1),cthree(inde3,2),'kd','MarkerSize',15,'DisplayName','UEVBS3 with SINR');
    plot(cfour(inde4,1),cfour(inde4,2),'kd','MarkerSize',15,'DisplayName','UEVBS4 with SINR');
    plot(cfour(I4,1),cfour(I4,2),'b*','MarkerSize',15,'DisplayName','UEVBS4');
    plot(cfive(inde5,1),cfive(inde5,2),'kd','MarkerSize',15,'DisplayName','UEVBS5 with SINR');
    plot(csix(inde6,1),csix(inde6,2),'kd','MarkerSize',15,'DisplayName','UEVBS6 with SINR');
    plot(cseven(inde7,1),cseven(inde7,2),'kd','MarkerSize',15,'DisplayName','UEVBS7 with SINR');
    
    plot(gnbc(1),gnbc(2),'hk','MarkerSize',15,'DisplayName','Base Station');
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    title('Selection of UE-VBS');
    hold off;
    lgd=legend;
    lgd.NumColumns=2;
    
    rp1 = recpow(inde1, cone(inde1,1), cone(inde1,2), cone);
    rp2 = recpow(inde2, ctwo(inde2,1), ctwo(inde2,2), ctwo);
    rp3 = recpow(inde3, cthree(inde3,1), cthree(inde3,2), cthree);
    rp4 = recpow(inde4, cfour(inde4,1), cfour(inde4,2), cfour);
    rp5 = recpow(inde5, cfive(inde5,1), cfive(inde5,2), cfive);
    rp6 = recpow(inde6, csix(inde6,1), csix(inde6,2), csix);
    rp7 = recpow(inde7, cseven(inde7,1), cseven(inde7,2), cseven);
    
    pc1 = abs(pc1h + 10^((-1*cone(inde1,4)-30)/10)*(numel(cone)-1) - rp1);
    pc2 = abs(pc2h + 10^((-1*ctwo(inde2,4)-30)/10)*(numel(ctwo)-1) - rp2);
    pc3 = abs(pc3h + 10^((-1*cthree(inde3,4)-30)/10)*(numel(cthree)-1) - rp3);
    pc4 = abs(pc4h + 10^((-1*cfour(inde4,4)-30)/10)*(numel(cfour)-1) - rp4);
    pc5 = abs(pc5h + 10^((-1*cfive(inde5,4)-30)/10)*(numel(cfive)-1) - rp5);
    pc6 = abs(pc6h + 10^((-1*csix(inde6,4)-30)/10)*(numel(csix)-1) - rp6);
    pc7 = abs(pc7h + 10^((-1*cseven(inde7,4)-30)/10)*(numel(cseven)-1) - rp7);

    pctot = pc1+pc2+pc3+pc4+pc5+pc6+pc7;
    
    powcons7 = 10*log10(pctot*1000);
    
    figure;
    bar([1 2 3 4 5 6 7], [pc1 pc2 pc3 pc4 pc5 pc6 pc7], 0.4);
    ylabel('Power consumed in each cluster in dBm');
    xlabel('Cluster number');
    title('Power consumption');

    figure;
    thh=[];
    thwa=[throughput(maxx1,100e7) throughput(maxx2,100e7) throughput(maxx3,100e7) throughput(maxx4,100e7) throughput(maxx5,100e7) throughput(maxx6,100e7) throughput(maxx7,100e7)];
    
    %Throughput of the algorithm selection
    thh=[th4temp;thwoa;thwa];
    barX = categorical({'Cluster 1','Cluster 2','Cluster3','Cluster 4','Cluster 5','Cluster 6','Cluster 7'});
    barX= reordercats(barX,{'Cluster 1','Cluster 2','Cluster3','Cluster 4','Cluster 5','Cluster 6','Cluster 7'});
    bar(barX,log10(transpose(thh)),0.4);
    ylabel('Sum rate in bits per second');
    legend({'LTE','Random Selection','Algorithm Selection'});
    title('Sum rate of the UE VBS network for 7 clusters');

    figure;
    bar([1 2 3 4 5 6 7],[q1 q2 q3 q4 q5 q6 q7],0.4);
    ylabel("Jain's fairness Index");
    xlabel("Cluster Number");
    title('Quality of Service for each cluster');
    zz=mean(th4temp,'all');
    zz1=mean(thwa,'all');
    disp(zz);
    disp(zz1);
    tp7=[zz zz1];
    %Returning the average throughput of the network
    
end