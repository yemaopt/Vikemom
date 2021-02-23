clear
clc
addpath ./Data
addpath ./Functions
load objectspoke251.mat

fobject=fftshift(fft2(ifftshift(object)));
pixelNum=size(object,1);

na=1.4;
pixelSize=25;
wavelength=500;
cutoff=na*pixelNum*pixelSize/wavelength;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctf=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2<=1);
ctfSignificantPix=numel(find(abs(ctf)>eps(class(ctf))));
ifftscale=numel(ctf)/ctfSignificantPix;
apsf=fftshift(ifft2(ifftshift(ctf)));
ipsf=ifftscale*abs(apsf).^2;
OTF=abs(fftshift(fft2(ifftshift(ipsf))));

ipsfConfocal=ipsf.*ipsf;
otfConfocal=abs(fftshift(fft2(ifftshift(ipsfConfocal))));

imgConfocal=fftshift(ifft2(ifftshift(otfConfocal.*fobject)));
maxNum=max(max(imgConfocal));

mask=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff/4).^2+((kxx-(pixelNum+1)/2)/cutoff/4).^2<=1);

for noisePercent=[0 0.1 0.2 0.3 0.4 0.5]
    imgConfocal_Noise=imnoise(imgConfocal./maxNum,'gaussian',0,noisePercent^2).*maxNum;
    fimgConfocal=fftshift(fft2(ifftshift(imgConfocal_Noise)));
    imgConfocal_Noise_filtered=fftshift(ifft2(ifftshift(fimgConfocal.*mask)));
    figure
    subplot(121);imagesc(imgConfocal_Noise);title(['Noise ',num2str(noisePercent*100),'%'])
    subplot(122);imagesc(imgConfocal_Noise_filtered);title(['Noise ',num2str(noisePercent*100),'% filtered'])
    colormap(hot)
end
