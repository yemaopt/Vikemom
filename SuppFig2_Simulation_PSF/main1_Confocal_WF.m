clear
clc
load SimulatedBeads251.mat

fObject=fftshift(fft2(ifftshift(object)));
[pixelNum,~]=size(object);

na=1.4;
pixelSize=10;
wavelength=500;
cutoff=na*pixelNum*pixelSize/wavelength;
[kyy,kxx]=meshgrid(1:pixelNum,1:pixelNum);
ctf=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2<=1);
ctfSignificantPix=numel(find(abs(ctf)>eps(class(ctf))));
ifftscale=numel(ctf)/ctfSignificantPix;
apsfWF=fftshift(ifft2(ifftshift(ctf)));
ipsfWF=ifftscale*abs(apsfWF).^2;
otfWF=abs(fftshift(fft2(ifftshift(ipsfWF))));

ipsfConfocal=ipsfWF.^2;
otfConfocal=abs(fftshift(fft2(ifftshift(ipsfConfocal))));

imageWF=fftshift(ifft2(ifftshift(fObject.*otfWF)));
imageConfocal=fftshift(ifft2(ifftshift(fObject.*otfConfocal)));

figure
imagesc(imageWF);colormap(hot);
set(gcf,'color','white');
axis off
title('WF')

figure
imagesc(imageConfocal);colormap(hot)
set(gcf,'color','white');
axis off
title('Confocal')
