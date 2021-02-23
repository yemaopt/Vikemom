clear
clc
addpath ./Data
addpath ./Functions
load STEDSimulatedBeads251.mat
load PSFHollow.mat
ipsfhollow=ipsfhollow./max(max(ipsfhollow));
ipsfhollow=imresize(ipsfhollow,[501,501]);
ipsfhollow=ipsfhollow((251-125):(251+125),(251-125):(251+125));
ipsfhollow=ipsfhollow./max(max(ipsfhollow));
fobject=fftshift(fft2(ifftshift(object)));

pixelNum=size(object,1);
na=1.4;
pixelSize=5;

wavelengthDet=500;
cutoffDet=na*pixelNum*pixelSize/wavelengthDet;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctfDet=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoffDet).^2+((kxx-(pixelNum+1)/2)/cutoffDet).^2<=1);
ctfSignificantPixDet=numel(find(abs(ctfDet)>eps(class(ctfDet))));
ifftscaleDet=numel(ctfDet)/ctfSignificantPixDet;
apsfDet=fftshift(ifft2(ifftshift(ctfDet)));
ipsfDet=ifftscaleDet*abs(apsfDet).^2;
OTFDet=abs(fftshift(fft2(ifftshift(ipsfDet))));

wavelengthExc=500;
cutoffExc=na*pixelNum*pixelSize/wavelengthExc;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctfExc=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoffExc).^2+((kxx-(pixelNum+1)/2)/cutoffExc).^2<=1);
ctfSignificantPixExc=numel(find(abs(ctfExc)>eps(class(ctfExc))));
ifftscaleExc=numel(ctfExc)/ctfSignificantPixExc;
apsfExc=fftshift(ifft2(ifftshift(ctfExc)));
ipsfExc=ifftscaleExc*abs(apsfExc).^2;
s=5;
ipsfExc=ipsfExc.*exp(-ipsfhollow.*s);
OTFExc=abs(fftshift(fft2(ifftshift(ipsfExc))));

ipsfPointDet=ipsfDet.*ipsfExc;
OTFPointDet=abs(fftshift(fft2(ifftshift(ipsfPointDet))));

imgInfiniteDet=fftshift(ifft2(ifftshift(fobject.*OTFExc)));
imgPointDet=fftshift(ifft2(ifftshift(fobject.*OTFPointDet)));

figure
imagesc(imgInfiniteDet);colormap(hot);
set(gcf,'color','white');title('Infinite Detector')
axis off

figure
imagesc(imgPointDet);colormap(hot)
set(gcf,'color','white');;title('Point Detector')
axis off
