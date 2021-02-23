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
OTF1=abs(fftshift(fft2(ifftshift(ipsf1))));

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
OTF2=abs(fftshift(fft2(ifftshift(ipsf2))));

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
OTF3=abs(fftshift(fft2(ifftshift(ipsf3))));

imgWF=fftshift(ifft2(ifftshift(OTF2.*fobject2)))+fftshift(ifft2(ifftshift(OTF1.*fobject1)))+fftshift(ifft2(ifftshift(OTF3.*fobject3)));
figure
imagesc(imgWF(151-125:151+125,151-125:151+125));colormap(hot)
title('Wide field')

ipsfConfocal1=ipsf1.^2;
OTFConfocal1=fftshift(fft2(ifftshift(ipsfConfocal1)));
ipsfConfocal2=ipsf2.^2;
OTFConfocal2=fftshift(fft2(ifftshift(ipsfConfocal2)));
ipsfConfocal3=ipsf2.^2;
OTFConfocal3=fftshift(fft2(ifftshift(ipsfConfocal3)));

imgConfocal=fftshift(ifft2(ifftshift(OTFConfocal2.*fobject2)))+fftshift(ifft2(ifftshift(OTFConfocal1.*fobject1)))+fftshift(ifft2(ifftshift(OTFConfocal3.*fobject3)));
figure
imagesc(imgConfocal(151-125:151+125,151-125:151+125));colormap(hot)
title('Confocal 0AU')

