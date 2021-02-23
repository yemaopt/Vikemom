clear
clc
addpath ./Data
addpath ./Functions
load spokeobject251full.mat

fObject=fftshift(fft2(ifftshift(object)));

na=1.4;
wavelength=500;
pixelSize=25;
pixelNum=size(object,1);

cutoff=na*pixelNum*pixelSize/wavelength;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctf=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2<=1);
ctfSignificantPix=numel(find(abs(ctf)>eps(class(ctf))));
ifftscale=numel(ctf)/ctfSignificantPix;
apsf=fftshift(ifft2(ifftshift(ctf)));
ipsf=ifftscale.*abs(apsf).^2;
OTF=abs(fftshift(fft2(ifftshift(ipsf))));

for factor = [0.5, 0.75, 1, 1.25, 1.5, 2, 3]
    pinholeRadius=round(factor*0.61*wavelength/na/pixelSize);   
    ii=0;
    for ydd=-pinholeRadius:1:pinholeRadius
        for xdd=-pinholeRadius:1:pinholeRadius
            ii=ii+1;
            psfeff1=ipsf.*moveElement(ipsf,-xdd,-ydd);
            otfeff1=fftshift(fft2(ifftshift(psfeff1)));
            sample(:,:,ii)=fftshift(ifft2(ifftshift(otfeff1.*fObject)));
            xd(ii)=xdd;
            yd(ii)=ydd;
        end
    end  
    num=0;
    cutoff_kg=2*2*pi/wavelength*na*pixelSize;
    for kg=(0.9/4:0.9/4:0.9).*cutoff_kg
        for theta=0:(5/180*pi):pi
            num=num+1
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
    [Iobj, OTFexc] = func1_VikemomRecon(im);    
end
