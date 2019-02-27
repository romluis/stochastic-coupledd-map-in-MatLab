close all
clear all
m=200; lt=500; %system size and numer of steps
nu=40.0; gamma=0.0; xs=0.7; rr=0.7; ab=0.001;  TT=270; at=ab*TT; %parameters

l=8; lgh=2*l+1; %parameters average 
H = fspecial('average', [lgh lgh]);
%I = imfilter(I, H);


csr=ones(m,m); sgn=zeros(lt,1); % initial csr at one in all entries. $sgn will kepp the average

for it=1:lt
 pcsr=1.0./(1.0+(xs./csr).^nu); %matris of probability function at this step
 rnb=rand(m,m); %random matrix of numbers
 sp=rnb<pcsr; %logical comaprison indicating where there is a sp. Matrix of 0 and 1. 1 means spark
 csra=csr.*(1.0-rr.*sp); %releae step
 
 %csri=conv2(csra,kernel,'same');
 csri=imfilter(csra,H,'replicate'); %average step where matlab makes the average with neighbors using matrix H
                                     %matrix H is the filter. The option
                                     %replicates gives the boundary
                                     %condition..I set it at the one by
                                     %default.
 %csri=csra;
 %csri=mean2(csra).*ones(m,m);
 uu=at*(1.0-csri); %uptake
 
 csr=csri+uu; %new csr after step
 sgn(it)=mean2(csr); %mean of the whole matrix is saved

end

 figure(1)
imagesc(csr,[0.4 0.85]);

figure(2)
plot(sgn)

