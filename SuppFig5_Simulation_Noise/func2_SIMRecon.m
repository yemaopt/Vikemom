function [IReconSIM] = func2_SIMRecon(im)

na=1.4;
pixelSize=25;
wavelength=500;
pixelNum=size(im,1);
cutoff=na*pixelNum*pixelSize/wavelength;
[kxx,kyy]=meshgrid(1:pixelNum,1:pixelNum);
ctf=ones(pixelNum,pixelNum).*(((kyy-(pixelNum+1)/2)/cutoff).^2+((kxx-(pixelNum+1)/2)/cutoff).^2<=1);
ctfSignificantPix=numel(find(abs(ctf)>eps(class(ctf))));
ifftscale=numel(ctf)/ctfSignificantPix;
apsf=fftshift(ifft2(ifftshift(ctf)));
ipsf=ifftscale*abs(apsf).^2;
OTF=fftshift(fft2(ifftshift(ipsf)));

N=3;
cutoff_kg=4*pi*na*pixelSize/wavelength;
kg=0.9*cutoff_kg;

for n=0:N-1
    kx(n+1)=kg*cos(pi/N*n);
    ky(n+1)=kg*sin(pi/N*n);
end

phi1=0;
phi2=2*pi/3;
phi3=4*pi/3;

weight=(interp2(1:pixelNum,1:pixelNum,OTF,kg*pixelNum/2/pi+(pixelNum+1)/2,(pixelNum+1)/2));
coeffMatrix=[1 0.5*weight*exp(-1i*phi1) 0.5*weight*exp(1i*phi1);1 0.5*weight*exp(-1i*phi2) 0.5*weight*exp(1i*phi2);1 0.5*weight*exp(-1i*phi3) 0.5*weight*exp(1i*phi3)];
coeffMatrixInv=inv(coeffMatrix);

for mm=1:N
    D1{mm}=fftshift(fft2(ifftshift(im(:,:,1+3*(mm-1)))));
    D2{mm}=fftshift(fft2(ifftshift(im(:,:,2+3*(mm-1)))));
    D3{mm}=fftshift(fft2(ifftshift(im(:,:,3+3*(mm-1)))));
    D0{mm}=coeffMatrixInv(1,1).*D1{mm}+coeffMatrixInv(1,2).*D2{mm}+coeffMatrixInv(1,3).*D3{mm};
    Dn1{mm}=coeffMatrixInv(2,1).*D1{mm}+coeffMatrixInv(2,2).*D2{mm}+coeffMatrixInv(2,3).*D3{mm};
    Dp1{mm}=coeffMatrixInv(3,1).*D1{mm}+coeffMatrixInv(3,2).*D2{mm}+coeffMatrixInv(3,3).*D3{mm};
end

%directly sum in Spatial domain
[XX,YYY]=meshgrid(-(pixelNum-1)/2:(pixelNum-1)/2,-(pixelNum-1)/2:(pixelNum-1)/2);
IReconSIM=0;
for mm=1:N
    In1=fftshift(ifft2(ifftshift(Dn1{mm})));
    Ip1=fftshift(ifft2(ifftshift(Dp1{mm})));
    IReconSIM=IReconSIM+In1.*exp(1i.*(kx(mm).*XX+ky(mm).*YYY))+Ip1.*exp(1i.*(-kx(mm).*XX-ky(mm).*YYY));
end
I0=fftshift(ifft2(ifftshift(D0{1})));
IReconSIM=IReconSIM+I0;
figure
imagesc(abs(IReconSIM));colormap(hot)
axis off



