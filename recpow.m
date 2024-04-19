function rp = recpow(UEn,xb,yb,X)
    lambda=4; %Frequency dependent loss and Path loss exponent are same in this case
    PL=4; %Path loss exponent
    dr=[];
    LightSpeedC=3e8; 
    Freq=5e9;
    TXAntennaGain=1;%db
    RXAntennaGain=1;%d
    PathLossExponent=4;%Line Of sight
    Wavelength=LightSpeedC/Freq;
    Dref=10;
    rstate = randn('state');
    GaussRandom= (randn*0.1+0);
    nvar= -174+10*log10(100e7);
    prw = [];
    for i=1:numel(X(:,1))
        if(i~=UEn && (X(i,1)~=xb) && (X(i,2)~=yb))
            rtm=sqrt(((X(i,1)-xb)^2)-((X(i,2)-yb)^2)); %Distance b/w two UEs
            M = Wavelength / (4 * pi * rtm);
            Pr01=X(i,4) + TXAntennaGain + RXAntennaGain- (20*log10(1/M));
            PrX1=Pr01+(10*PathLossExponent* log10(X(i,6)/Dref))+GaussRandom;
            prw = [prw 10^((-1*PrX1-30)/10)];
        end
    end
    rp = sum(prw);
end
