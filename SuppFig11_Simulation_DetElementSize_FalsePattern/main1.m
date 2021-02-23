clear
clc
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

pinholeRadius=round(1.25*0.61*wavelength/na/pixelSize);
for ydd=-pinholeRadius:1:pinholeRadius+1
    for xdd=-pinholeRadius:1:pinholeRadius+1
        psfeff{ydd+pinholeRadius+1,xdd+pinholeRadius+1}(:,:)=ipsf.*moveElement(ipsf,-xdd,-ydd);
        otfeff{ydd+pinholeRadius+1,xdd+pinholeRadius+1}(:,:)=fftshift(fft2(ifftshift(psfeff{ydd+pinholeRadius+1,xdd+pinholeRadius+1}(:,:))));
        sample{ydd+pinholeRadius+1,xdd+pinholeRadius+1}(:,:)=fftshift(ifft2(ifftshift(otfeff{ydd+pinholeRadius+1,xdd+pinholeRadius+1}(:,:).*fobject)));
        yd(ydd+pinholeRadius+1,xdd+pinholeRadius+1)=ydd;
        xd(ydd+pinholeRadius+1,xdd+pinholeRadius+1)=xdd;
    end
end

for yii=1:pixelNum
    for xjj=1:pixelNum
        for kk=1:24
            for mm=1:24
                sample1{yii,xjj}(kk,mm)=sample{kk,mm}(yii,xjj);
            end
        end
    end
end

cutoff_kg=4*pi*na*pixelSize/wavelength;
for binningRatio=[1 1/2 1/3 1/4 1/6 1/8]
    binningRatio
    num=0;
    xdbinning=0;
    ydbinning=0;
    pattern=0;
    im=0;
    clear sample2
    for kg=(0.9/4:0.9/4:0.9).*cutoff_kg
        for theta=10*(5/180*pi)
            num=num+1;
            kx=kg*cos(theta);
            ky=kg*sin(theta);
            xdbinning=imresize(xd(:,:),binningRatio,'box');
            ydbinning=imresize(yd(:,:),binningRatio,'box');
            for yscan=-(pixelNum-1)/2:(pixelNum-1)/2
                for xscan=-(pixelNum-1)/2:(pixelNum-1)/2          
                    pattern=1+cos(kx.*(xdbinning+xscan)+ky.*(ydbinning+yscan));
                    xx=xscan+(pixelNum-1)/2+1;
                    yy=yscan+(pixelNum-1)/2+1;                
                    sample2{yy,xx}(:,:)=imresize(sample1{yy,xx}(:,:),binningRatio,'box')./binningRatio./binningRatio;
                    im(yy,xx,num)=sum(sum(sample2{yy,xx}(:,:).*pattern(:,:)));
                end
            end
        end
        figure;imagesc(im(:,:,num));colormap(hot)
        title(['Detector element size ',num2str(1.25*0.61/na/11/binningRatio),'\lambda'])
        pause(0.01)
    end
end



