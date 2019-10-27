import DVRPackage.*
%run('C:\Users\liny\OneDrive\Matlab functions me\UnitsConstants.m')

%Set trap parameters
mass = UnitsConstants.mRb87;
w0 = 707*UnitsConstants.nm;
asep = 900*UnitsConstants.nm/w0; % in the unit of the waist
V0 = 13.2*UnitsConstants.kHz; % from 21MHz/7.15mW*0.0045mW
Nbands=4;

%%%%%%%%Double-well Gaussian potential
%range of DVR space (units of w_0)
ax=2.26;
dx=0.05;
[xvals,Nx] = GetGridDx(ax,dx); %xvals are the grid valus, and Nx is the number of grids

BareScal=EnergyFromLengthScale(w0,mass); %Scaling of second derivative
Gaussf=@(x) -V0*(exp(-2.0*(x-0.5*asep).^2))-V0*(exp(-2.0*(x+0.5*asep).^2)); %Double-Gaussian potential with x in waist units
% set the location of the Gaussian to be +-0.5asep so that total sep is asep
% all lengths in the unit of the waist so that I~Exp(-2*r^2/w0^2)
% 9.80100(2) https://digital.library.unt.edu/ark:/67531/metadc13155/m2/1/high_res_d/NBS%20Monograph%20107.pdf

%Get states/energies
[evecs,G_evals] = DVR_1D(xvals,BareScal,'p',Nbands,Gaussf); %p is solving the symetric wavefunctions? since it's even
[ovecs,G_ovals] = DVR_1D(xvals,BareScal,'m',Nbands,Gaussf); %m is solving the assymetric wavefunctions? since it's odd

disp('Even Double-well energies')
G_evals/UnitsConstants.kHz
disp('Odd Double-well energies')
G_ovals/UnitsConstants.kHz

% eigen wavefuncitons. here G_gs / G_es are the states for the ground band,
% and G_gs2 / G_es2 are the states for the excited band
G_gs = Unpack1DState(evecs(1,:),'p')/sqrt(dx); 
G_es = Unpack1DState(ovecs(1,:),'m')/sqrt(dx);
G_gs2 = Unpack1DState(evecs(2,:),'p')/sqrt(dx);
G_es2 = Unpack1DState(ovecs(2,:),'m')/sqrt(dx);
long_G_xvals=dx*[-Nx:Nx];

figure(11)
subplot(2,1,1)    
plot(long_G_xvals,G_gs,'r')
hold
plot(long_G_xvals,G_es,'b')
subplot(2,1,2)    
plot(long_G_xvals,G_gs2,'r')
hold
plot(long_G_xvals,G_es2,'b')


