function [Iobj, OTFexc] = func1_VikemomRecon(im)

pixelNum=size(im,1);
delta=0.0000001;
na=1.4;
pixelSize=25;
wavelength=500;

cutoff=na*pixelNum*pixelSize/wavelength;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctfdet=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2<=1);
ctfdetSignificantPix=numel(find(abs(ctfdet)>eps(class(ctfdet))));
ifftscaledet=numel(ctfdet)/ctfdetSignificantPix;
apsfdet=fftshift(ifft2(ifftshift(ctfdet)));
ipsfdet=ifftscaledet*abs(apsfdet).^2;
OTFdet=fftshift(fft2(ifftshift(ipsfdet)));

[xx,yy]=meshgrid(-(pixelNum-1)/2:(pixelNum-1)/2,-(pixelNum-1)/2:(pixelNum-1)/2);
num=0;
cutoff_kg=2*2*pi/wavelength*na*pixelSize;
for kg=(0.9/4:0.9/4:0.9).*cutoff_kg
    for theta=0:(5/180*pi):pi
        num=num+1;
        kx=kg*cos(theta);
        ky=kg*sin(theta);
        pattern=1+cos(kx.*(xx)+ky.*(yy));
        fpattern=fftshift(fft2(ifftshift(pattern)));
        fpotf=fpattern.*OTFdet;
        patterneff(:,:,num)=fftshift(ifft2(ifftshift(fpotf)));
    end
end

[kxx, kyy]=meshgrid(1:pixelNum,1:pixelNum);
mask=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff/2).^2+((kxx-(pixelNum+1)/2)/cutoff/2).^2<=1);

Iobj=im(:,:,1);OTFexc=OTFdet;num=0;
for jj=1:60
    jj
    for num=1:size(im,3)
        fobj=fftshift(fft2(ifftshift(Iobj)));
        Iop=Iobj.*patterneff(:,:,num);
        fop=fftshift(fft2(ifftshift(Iop)));
        fim=fftshift(fft2(ifftshift(im(:,:,num))));
        ftemp1=mask.*(fim-OTFexc.*fop);
        Itemp1=fftshift(ifft2(ifftshift(ftemp1)));
        Itemp2=Itemp1.*patterneff(:,:,num);
        ftemp2=fftshift(fft2(ifftshift(Itemp2)));
        fobj=fobj+ftemp2;
        OTFexc=OTFexc+abs(fop).*conj(fop).*(fim-OTFexc.*fop)./max(max(abs(fop)))./(abs(fop).^2+delta);
        OTFexc=OTFexc.*mask;
        Iobj=fftshift(ifft2(ifftshift(fobj)));
    end
end

figure(102)
subplot(121);imagesc(abs(Iobj));
subplot(122);imagesc(abs(OTFexc));
colormap(hot)
axis off
pause(0.0001)

end


