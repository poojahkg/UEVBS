function s=sinr4g(xb,yb,X,p,Pr4g,dist,num_ue)
   hnm=randn(1)+j*randn(1); %Rayleigh fading channel gain
   lambda=4;
   PL=4;
   Freq = 1.9e9;
   LightSpeedC = 3e8;
   TXAntennaGain=1;%db
   RXAntennaGain=1;%d
   PathLossExponent=4;%Line Of sight
   Wavelength4g=LightSpeedC/Freq;
   Dref=10;
   rstate = randn('state');
   GaussRandom= (randn*0.1+0);
   dr=[];
   nvar= -174+log10(20e7);
   nvar = 10^((nvar-30)/10);
   for i=1:num_ue
       if((X(i,1)~=xb) && (X(i,2)~=yb))
           M = Wavelength4g / (4 * pi * dist);
           Pr0=43 + TXAntennaGain + RXAntennaGain- (20*log10(1/M)); %43dBm is the standard BS power
           Pr4g=[Pr4g Pr0+(10*PathLossExponent* log10(dist/Dref))+GaussRandom];
           rtm=sqrt(((X(i,1)-xb)^2)-((X(i,2)-yb)^2)); %In 4G, the transmitter is the BS
           dr=[dr nvar+(randn(1)+j*randn(1))*(10^(-1*Pr4g(i)/10)/1000)*lambda*(rtm^(-1*lambda))];
       end
   end
   % Using the same formula, SINR is computed with a fixed
   % transmitter,i.e, the Base Station.
   Inter = sum(dr);
   rnm=dist;
   sign=hnm*p*lambda*(rnm^(-1*lambda));
   s = abs(sign/Inter);
end