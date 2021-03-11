clear
clc
addpath ./Data
addpath ./Functions

load background301.mat
load spokeobject301half.mat

object1=object;
object2=background;
object3=background;

pixelNum=size(object,1);

fobject1=fftshift(fft2(ifftshift(object1)));
fobject2=fftshift(fft2(ifftshift(object2)));
fobject3=fftshift(fft2(ifftshift(object3)));

na=1.4;
wavelength=500;
pixelSize=wavelength/20;
n=1.518;
alpha=asin(na/n);

cutoff=na*pixelNum*pixelSize/wavelength;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctf1=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2<=1);
ctfSignificantPix1=numel(find(abs(ctf1)>eps(class(ctf1))));
ifftscale1=numel(ctf1)/ctfSignificantPix1;
apsf1=fftshift(ifft2(ifftshift(ctf1)));
ipsf1=ifftscale1.*abs(apsf1).^2;
OTF1=fftshift(fft2(ifftshift(ipsf1)));
OTF1=abs(OTF1);

z0=-0.8*wavelength;
u0=8*pi/wavelength*z0*n*sin(alpha/2)^2;
cutoff=na*pixelNum*pixelSize/wavelength;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctf2=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2<=1);
ctf2=ctf2.*exp(1i.*u0./2.*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2));
ctfSignificantPix2=numel(find(abs(ctf2)>eps(class(ctf2))));
ifftscale2=numel(ctf2)/ctfSignificantPix2;
apsf2=fftshift(ifft2(ifftshift(ctf2)));
ipsf2=ifftscale2.*abs(apsf2).^2;
OTF2=fftshift(fft2(ifftshift(ipsf2)));
OTF2=abs(OTF2);

z0=0.8*wavelength;
u0=8*pi/wavelength*z0*n*sin(alpha/2)^2;
cutoff=na*pixelNum*pixelSize/wavelength;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctf3=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2<=1);
ctf3=ctf3.*exp(1i.*u0./2.*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2));
ctfSignificantPix3=numel(find(abs(ctf3)>eps(class(ctf3))));
ifftscale3=numel(ctf3)/ctfSignificantPix3;
apsf3=fftshift(ifft2(ifftshift(ctf3)));
ipsf3=ifftscale3.*abs(apsf3).^2;
OTF3=fftshift(fft2(ifftshift(ipsf3)));
OTF3=abs(OTF3);

pinholeRadius=round(1.25*0.61*wavelength/na/pixelSize);
ii=0;
imgConfocal1p25AU=0;
for ydd=-pinholeRadius:1:pinholeRadius
    for xdd=-pinholeRadius:1:pinholeRadius
        ii=ii+1;
        psfeff1=ipsf1.*moveElement(ipsf1,-xdd,-ydd);
        psfeff2=ipsf2.*moveElement(ipsf2,-xdd,-ydd);
        psfeff3=ipsf3.*moveElement(ipsf3,-xdd,-ydd);       
        otfeff1=fftshift(fft2(ifftshift(psfeff1)));
        otfeff2=fftshift(fft2(ifftshift(psfeff2)));
        otfeff3=fftshift(fft2(ifftshift(psfeff3))); 
        if ydd^2+xdd^2 <= pinholeRadius^2
           imgConfocal1p25AU=imgConfocal1p25AU+fftshift(ifft2(ifftshift(otfeff1.*fobject1)))+fftshift(ifft2(ifftshift(otfeff2.*fobject2)))+fftshift(ifft2(ifftshift(otfeff3.*fobject3)));    
        end
    end
end

figure
imagesc(imgConfocal1p25AU(151-125:151+125,151-125:151+125));colormap(hot)
title('Confocal 1.25AU')
