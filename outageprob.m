function op = outageprob(A,I,tow)
    Nr= numel(A(:,1))-1;
    Dr= Nr +((1+(2/29)*(0.12/A(I,3))^(1/2))*((pi/2)*sqrt(tow)));
    op= 1-(Nr/Dr);
end