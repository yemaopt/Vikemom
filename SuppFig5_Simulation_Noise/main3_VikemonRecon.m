clear
clc
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

pinholeRadius=round(1.25*0.61*wavelength/na/pixelSize);

for noisePercent=0:0.1:0.5
    noisePercent
    ii=0;
    for ydd=-pinholeRadius:1:pinholeRadius
        for xdd=-pinholeRadius:1:pinholeRadius
            ii=ii+1;
            psfeff(:,:,ii)=ipsf.*moveElement(ipsf,-xdd,-ydd);
            otfeff(:,:,ii)=fftshift(fft2(ifftshift(psfeff(:,:,ii))));
            sample(:,:,ii)=fftshift(ifft2(ifftshift(otfeff(:,:,ii).*fobject)));
            maxNum=max(max((sample(:,:,ii))));
            sample(:,:,ii)=imnoise(sample(:,:,ii)./maxNum,'gaussian',0,noisePercent^2).*maxNum;
            xd(ii)=xdd;
            yd(ii)=ydd;
        end
    end
    num=0;
    cutoff_kg=4*pi*na*pixelSize/wavelength;
    for kg=(0.9/4:0.9/4:0.9).*cutoff_kg
        for theta=0:(5/180*pi):pi
            num=num+1;
            kx=kg*cos(theta);
            ky=kg*sin(theta);
            for yscan=-(pixelNum-1)/2:(pixelNum-1)/2
                for xscan=-(pixelNum-1)/2:(pixelNum-1)/2
                    xx=xscan+(pixelNum-1)/2+1;
                    yy=yscan+(pixelNum-1)/2+1;
                    pattern=1+cos(kx.*(xd+xscan)+ky.*(yd+yscan));
                    im(yy,xx,num)=pattern*squeeze(double(sample(yy,xx,:)));
                end
            end
        end
    end
    [Iobj,OTFex]=func1_VikemomRecon(im);
end
