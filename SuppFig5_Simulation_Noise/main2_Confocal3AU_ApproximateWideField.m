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

pinholeRadius=round(3*0.61*wavelength/na/pixelSize);

for noisePercent=0:0.1:0.5
    imgWF=0;
    noisePercent
    for ydd=-pinholeRadius:1:pinholeRadius
        for xdd=-pinholeRadius:1:pinholeRadius
            psfeff=ipsf.*moveElement(ipsf,-xdd,-ydd);
            otfeff=fftshift(fft2(ifftshift(psfeff)));
            sample=fftshift(ifft2(ifftshift(otfeff.*fobject)));
            maxnum=max(max((sample)));
            temp=imnoise(sample./maxnum,'gaussian',0,noisePercent^2).*maxnum;
            imgWF=imgWF+temp;
        end
    end
    figure
    imagesc(imgWF);colormap(hot);axis off
    title(['Noise ',num2str(noisePercent*100),'%'])
end
