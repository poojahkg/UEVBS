%The algorithm implemented in all the cluster functions are the same with
%the only difference being the cluster size. Hence, the explanation will be
%given in this function block alone.
function [tp2, powcons2] = clustere2(X,disratee,num_ue)
    disrate=disratee;
    figure;
    alpha=0.5;
    beta=0.5;
    gnbc = [80,60]; %gNB - Base Station
    d=pdist2(gnbc,X,'euclidean')*20;
    clusterX = kmeans(X,2);
    plot(X(clusterX==1,1),X(clusterX==1,2),'r.','MarkerSize',15)
    hold on
    plot(X(clusterX==2,1),X(clusterX==2,2),'g.','MarkerSize',15)
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
    %The arrays defined above contain the UE details for each cluster
    
    pr1w =[];
    pr2w = [];
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
        else  
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
        end
    end
    %Power consumed in each cluster = (No.of UE in each cluster)*(BS power in watts) 
    % - (Pr of each UE from BS) + (No UEs-1)*(Pr of UE VBS) - (Pr of each UE from UE VBS)
    pc1h = numel(cone)*PTx - sum(pr1w);
    pc2h = numel(ctwo)*PTx - sum(pr2w);
   
    hold on;
    e1=[];
    e2=[];

    %Arrays ei are the ones holding the indices (in X) of the eligible UE VBSs. i
    %ranges from 1 to cluster-size
    ss1=[];
    ss2=[];
    
    sp1=[];
    sp2=[];
    
    %Arrays ssi are the ones holding the SINR of all the eligible UE VBSs in each cluster
    ra=randi(numel(cone(:,1)));
    rb=randi(numel(ctwo(:,1)));
    
    %Selecting UE VBS at random for comparison
    sr1=[];
    sr2=[];

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
    ss1=reshape(ss1,numel(e1),numel(a1));
    ss2=reshape(ss2,numel(e2),numel(a2));
    %ssi is converted into a 2D matrix with each row indicating eUE and
    %columns with SINR values from a particular eUE to a UE in the same
    %cluster.

    [maxx1,tow1]=sumsinr(ss1);
    [maxx2,tow2]=sumsinr(ss2);
    
    %The maximum mean SINR the eUEs is selected. maxx1 contains SINR value
    %and tow1 contains the index (wrt e1) of eUE which is selected
    
    sp1=reshape(sp1,numel(e1),numel(a1));
    sp2=reshape(sp2,numel(e2),numel(a2));

    [m1,t1]=sumse(sp1);
    [m2,t2]=sumse(sp2);

    q1=jfi(sp1,t1);
    q2=jfi(sp2,t2);
    %These are the QoS values of UE VBS in every cluster.
    
    [rt1,tt]=sumsinr(sr1);
    [rt2,tt]=sumsinr(sr2);
    
    %Same procedure for random UE
    
    thwoa=[];
    thwoa=[throughput(rt1,100e7) throughput(rt2,100e7)];
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
    
    %inde contains the index (in X) of the selected UE VBS.
    ss4g=[];
    
    for i=1:num_ue
        s14=sinr4g(X(i,1),X(i,2),X,0.25,Pr4g,d(i),num_ue);
        ss4g=[ss4g s14];
    end
    %ss4g is used to contain the SINR value for 4G network
    disp(ss4g);
    maxx4g=mean(ss4g);
    th4temp=[throughput(maxx4g,20e7) throughput(maxx4g,20e7)]; 
    %throughput of 4G network
    plot(cone(I1,1),cone(I1,2),'r*','MarkerSize',15,'DisplayName','UEVBS1');
    plot(cone(inde1,1),cone(inde1,2),'kd','MarkerSize',15,'DisplayName','UEVBS1 with SINR');
    plot(ctwo(I2,1),ctwo(I2,2),'g*','MarkerSize',15,'DisplayName','UEVBS2');
    plot(ctwo(inde2,1),ctwo(inde2,2),'kd','MarkerSize',15,'DisplayName','UEVBS2 with SINR');
    plot(gnbc(1),gnbc(2),'hk','MarkerSize',15,'DisplayName','Base Station');
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    title('Selection of UE-VBS');
    hold off;
    lgd=legend;
    lgd.NumColumns=2;
    
    rp1 = recpow(inde1, cone(inde1,1), cone(inde1,2), cone);
    rp2 = recpow(inde2, ctwo(inde2,1), ctwo(inde2,2), ctwo);
    
    pc1 = abs(pc1h + 10^((-1*cone(inde1,4)-30)/10)*(numel(cone)-1) - rp1);
    pc2 = abs(pc2h + 10^((-1*ctwo(inde2,4)-30)/10)*(numel(ctwo)-1) - rp2);
    
    pctot = pc1+pc2;
    
    powcons2 = 10*log10(pctot*1000);
    
    figure;
    bar([1 2], [pc1 pc2], 0.4);
    ylabel('Power consumed in each cluster in dBm');
    xlabel('Cluster number');
    title('Power consumption');

    figure;
    thh=[];
    thwa=[throughput(maxx1,100e7) throughput(maxx2,100e7)];
    
    %Throughput of the algorithm selection
    thh=[th4temp;thwoa;thwa];
    barX = categorical({'Cluster 1','Cluster 2'});
    barX= reordercats(barX,{'Cluster 1','Cluster 2',});
    bar(barX,log10(transpose(thh)),0.4);
    ylabel('Sum rate in bits per second');
    legend({'LTE','Random Selection','Algorithm Selection'});
    title('Sum rate of the UE VBS network for 2 clusters');

    figure;
    bar([1 2],[q1 q2],0.4);
    ylabel("Jain's fairness Index");
    xlabel("Cluster Number");
    title('Quality of Service for each cluster');
    zz=mean(th4temp,'all');
    zz1=mean(thwa,'all');
    disp(zz);
    disp(zz1);
    tp2=[zz zz1];
    %Returning the average throughput of the network
end