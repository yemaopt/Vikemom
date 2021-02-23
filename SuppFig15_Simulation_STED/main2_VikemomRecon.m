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
OTFde=abs(fftshift(fft2(ifftshift(ipsfDet))));

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

pinholeRadius=round(3*0.61*wavelengthDet/na/pixelSize);

ii=0;
for ydd=-pinholeRadius:9:pinholeRadius
    for xdd=-pinholeRadius:9:pinholeRadius
        ii=ii+1;
        psfeff(:,:,ii)=ipsfExc.*moveElement(ipsfDet,-xdd,-ydd);
        otfeff(:,:,ii)=fftshift(fft2(ifftshift(psfeff(:,:,ii))));
        sample(:,:,ii)=fftshift(ifft2(ifftshift(otfeff(:,:,ii).*fobject)));
        xd(ii)=xdd;
        yd(ii)=ydd;
    end
end

num=0;
cutoff_kg=4*pi*na*pixelSize/wavelengthExc;
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

func1_VikemomRecon(im);
