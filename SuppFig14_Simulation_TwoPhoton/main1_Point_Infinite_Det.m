clear
clc
addpath ./Data
load 2PSimulatedBeads251.mat

fobject=fftshift(fft2(ifftshift(object)));
pixelNum=size(object,1);
na=1.4;
pixelSize=10;

wavelengthDet=500;
cutoffDet=na*pixelNum*pixelSize/wavelengthDet;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctfDet=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoffDet).^2+((kxx-(pixelNum+1)/2)/cutoffDet).^2<=1);
ctfSignificantPixDet=numel(find(abs(ctfDet)>eps(class(ctfDet))));
ifftscaleDet=numel(ctfDet)/ctfSignificantPixDet;
apsfDet=fftshift(ifft2(ifftshift(ctfDet)));
ipsfDet=ifftscaleDet*abs(apsfDet).^2;
OTFDet=abs(fftshift(fft2(ifftshift(ipsfDet))));

wavelengthExc=1000;
cutoffExc=na*pixelNum*pixelSize/wavelengthExc;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctfExc=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoffExc).^2+((kxx-(pixelNum+1)/2)/cutoffExc).^2<=1);
ctfSignificantPixExc=numel(find(abs(ctfExc)>eps(class(ctfExc))));
ifftscaleExc=numel(ctfExc)/ctfSignificantPixExc;
apsfExc=fftshift(ifft2(ifftshift(ctfExc)));
ipsfExc=ifftscaleExc*abs(apsfExc).^4;
OTFex=abs(fftshift(fft2(ifftshift(ipsfExc))));

ipsfPointDet=ipsfDet.*ipsfExc;
OTFPointDet=abs(fftshift(fft2(ifftshift(ipsfPointDet))));

imgInifiteDet=fftshift(ifft2(ifftshift(fobject.*OTFex)));
imgPointDet=fftshift(ifft2(ifftshift(fobject.*OTFPointDet)));

figure;imagesc(imgInifiteDet);colormap(hot);set(gcf,'color','white');title('Infinite Detector');axis off
figure;imagesc(imgPointDet);colormap(hot);set(gcf,'color','white');title('Point Detector');axis off
