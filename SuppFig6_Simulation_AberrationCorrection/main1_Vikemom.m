clear
clc
load objectspoke251.mat
addpath ./Functions

fobject=fftshift(fft2(ifftshift(object)));
[pixelNum,~]=size(object);

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

k=2*pi/wavelength;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
rou=(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2).^0.5;
costheta=(kyy-(pixelNum+1)/2)/cutoff./rou;
costheta(126,126)=0;
ctfexc=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2<=1);

% Spherical aberration
% W=0.5*wavelength;
% ctfexc=ctfexc.*exp(1i.*k.*W.*rou.^4);
%
% Astigmatism
% W=0.5*wavelength;
% ctfexc=ctfexc.*exp(1i.*k.*W.*rou.^2.*costheta.^2);
%
% Defocus
% W=0.3*wavelength;
% ctfexc=ctfexc.*exp(1i.*k.*W.*rou.^2);
%
% Mixed
% W1=0.2*wavelength;
% W2=0.2*wavelength;
% W3=0.1*wavelength;
% ctfexc=ctfexc.*exp(1i.*k.*W1.*rou.^4).*exp(1i.*k.*W2.*rou.^2.*costheta.^2).*exp(1i.*k.*W3.*rou.^2);

% field curvature
W=0.5*wavelength;
ctfexc=ctfexc.*exp(1i.*k.*W.*rou.^2);

ctfexcSignificantPix=numel(find(abs(ctfexc)>eps(class(ctfexc))));
ifftscaleexc=numel(ctfexc)/ctfexcSignificantPix;
apsfexc=fftshift(ifft2(ifftshift(ctfexc)));
ipsfexc=ifftscaleexc*abs(apsfexc).^2;
OTFexc=fftshift(fft2(ifftshift(ipsfexc)));

pinholeRadius=round(3*0.61*wavelength/na/pixelSize);

ii=0;
for ydd=-pinholeRadius:2:pinholeRadius
    for xdd=-pinholeRadius:2:pinholeRadius
        ii=ii+1;
        psfeff(:,:,ii)=ipsfexc.*moveElement(ipsfdet,-xdd,-ydd);
        otfeff(:,:,ii)=fftshift(fft2(ifftshift(psfeff(:,:,ii))));
        sample(:,:,ii)=fftshift(ifft2(ifftshift(otfeff(:,:,ii).*fobject)));
        xd(ii)=xdd;
        yd(ii)=ydd;
    end
end

num=0;
cutoff_kg=4*pi*na*pixelSize/wavelength;
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

figure(101);imagesc(abs(OTFexc));colormap(hot)
[IobjRecon,OTFexcRecon]=func1_VikemomRecon(im);

