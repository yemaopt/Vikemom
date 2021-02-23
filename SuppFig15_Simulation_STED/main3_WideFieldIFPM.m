clear
clc
addpath ./Data
addpath ./Functions
load STEDSimulatedbeads251.mat
load PSFHollow.mat
ipsfhollow=ipsfhollow./max(max(ipsfhollow));
ipsfhollow=imresize(ipsfhollow,[501,501]);
ipsfhollow=ipsfhollow((251-125):(251+125),(251-125):(251+125));
ipsfhollow=ipsfhollow./max(max(ipsfhollow));
fobject=fftshift(fft2(ifftshift(object)));

pixelNum=size(object,1);
na=1.4;
pixelSize=5;

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

wavelengthDet=500;
cutoffDet=na*pixelNum*pixelSize/wavelengthDet;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctfDet=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoffDet).^2+((kxx-(pixelNum+1)/2)/cutoffDet).^2<=1);
ctfSignificantPixDet=numel(find(abs(ctfDet)>eps(class(ctfDet))));
ifftscaleDet=numel(ctfDet)/ctfSignificantPixDet;
apsfDet=fftshift(ifft2(ifftshift(ctfDet)));
ipsfDet=ifftscaleDet*abs(apsfDet).^2;
OTFDet=abs(fftshift(fft2(ifftshift(ipsfDet))));

x=-(pixelNum-1)/2:(pixelNum-1)/2;
y=-(pixelNum-1)/2:(pixelNum-1)/2;
[xx yy]=meshgrid(x,y);

cutoff_kg=4*pi*na*pixelSize/wavelengthExc;
num=0;
k=0;
for k=(0.9/4:0.9/4:0.9).*cutoff_kg
     for theta=0:(5/180*pi):pi
         num=num+1;
         kx=k*cos(theta);
         ky=k*sin(theta);
         pattern=1+cos(kx.*(xx)+ky.*(yy)); 
         fp=fftshift(fft2(ifftshift(pattern)));
         fpotf=fp.*OTFDet;
         patterneff=fftshift(ifft2(ifftshift(fpotf)));
         temp1=object.*patterneff;
         ftemp1=fftshift(fft2(ifftshift(temp1)));
         ftemp2=ftemp1.*OTFExc;
         temp3=fftshift(ifft2(ifftshift(ftemp2)));       
         im(:,:,num)=temp3;
     end
end

func1_VikemomRecon(im);
figure;imagesc(OTFExc);colormap(hot)

 