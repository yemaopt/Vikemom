clear
clc

load modulatedImgSet.mat

pixelNum=size(im,1);
delta=0.0000001;
na=1.4;
pixelSize=10;
wavelength=500;
cutoff=na*pixelNum*pixelSize/wavelength;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctf=ones(pixelNum,pixelNum).*(((kxx-(pixelNum+1)/2)/cutoff).^2+((kyy-(pixelNum+1)/2)/cutoff).^2<=1);
ctfSignificantPix=numel(find(abs(ctf)>eps(class(ctf))));
ifftscale=numel(ctf)/ctfSignificantPix;
apsf=fftshift(ifft2(ifftshift(ctf)));
ipsf=ifftscale*abs(apsf).^2;
OTF=abs(fftshift(fft2(ifftshift(ipsf))));

[xx,yy]=meshgrid(-(pixelNum-1)/2:(pixelNum-1)/2,-(pixelNum-1)/2:(pixelNum-1)/2);
num=0;

cutoff_kg=4*pi*na*pixelSize/wavelength;
for kg=(0.9/4:0.9/4:0.9).*cutoff_kg
    for theta=0:(5/180*pi):pi
        num=num+1;
        kx=kg*cos(theta);
        ky=kg*sin(theta);
        pattern=1+cos(kx.*(xx)+ky.*(yy));
        fpattern=fftshift(fft2(ifftshift(pattern)));
        patterneff(:,:,num)=fftshift(ifft2(ifftshift(fpattern.*OTF)));
    end
end

[kxx, kyy]=meshgrid(1:pixelNum,1:pixelNum);
mask=ones(pixelNum,pixelNum).*(((kxx-(pixelNum+1)/2)/cutoff/2).^2+((kyy-(pixelNum+1)/2)/cutoff/2).^2<=1);

Iobj=im(:,:,1);
num=0;
for jj=1:20
    jj
    for num=1:size(im,3)
        fobj=fftshift(fft2(ifftshift(Iobj)));
        Iop=Iobj.*patterneff(:,:,num);
        fop=fftshift(fft2(ifftshift(Iop)));
        fim=fftshift(fft2(ifftshift(im(:,:,num))));
        ftemp1=mask.*(fim-OTF.*fop);
        Itemp1=fftshift(ifft2(ifftshift(ftemp1)));
        Itemp2=Itemp1.*patterneff(:,:,num);
        ftemp2=fftshift(fft2(ifftshift(Itemp2)));
        fobj=fobj+ftemp2;
        OTF=OTF+abs(fop).*conj(fop).*(fim-OTF.*fop)./max(max(abs(fop)))./(abs(fop).^2+delta);
        OTF=OTF.*mask;
        Iobj=fftshift(ifft2(ifftshift(fobj)));
        if mod(jj,2)==0
           Iobj(Iobj<0)=0;
        end
        figure(1)
        subplot(121);imagesc(abs(Iobj));
        subplot(122);imagesc(abs(OTF));
        colormap(hot)
        axis off
        pause(0.0001)
    end
end

