clear all
import DVRPackage.*
%run('C:\Users\liny\OneDrive\Matlab functions me\UnitsConstants.m')

%Set trap parameters
mass = UnitsConstants.mRb87;
w0 = 707*UnitsConstants.nm;
asep = 0.209*(80.405-76.095)*1000*UnitsConstants.nm/w0; % in the unit of the waist. conversion of 0.209(3) um/MHz
scale = 1.4; % scale 1.4 and scaleDelta 0.87 gives 4.11ms tunneling 2pi time
V0 = scale*13.2*UnitsConstants.kHz; % from 21MHz/7.15mW*0.0045mW
scaleDelta=0.87;
Delta= (0.53-0.47)/0.5*V0*scaleDelta;
Nbands=8;

%%%%%%%%Biased Double-well Gaussian potential
%range of DVR space (units of w_0)
ax=2.26;
dx=0.05;
[xvals,Nx] = GetASGridDx(ax,dx);

BareScal=EnergyFromLengthScale(w0,mass); %Scaling of second derivative
Gaussf=@(x) -(V0+0.5*Delta)*exp(-2.0*(x-0.5*asep).^2)-(V0-0.5*Delta)*exp(-2.0*(x+0.5*asep).^2)+mRb87*9.8*x*w0/hbar/2/pi; %Gaussian potential with x in waist units
%Get states/energies
[vecs,G_vals] = DVR_AS1D(xvals,BareScal,Nbands,Gaussf);

disp('Double-well energies')
G_vals/UnitsConstants.kHz

G_gs = vecs(1,:)/sqrt(dx);
G_es = vecs(2,:)/sqrt(dx);
G_gs2 = vecs(3,:)/sqrt(dx);
G_es2 = vecs(4,:)/sqrt(dx);
long_G_xvals=dx*[-Nx:Nx];

figure(1)
hold off
subplot(2,1,1)    
plot(long_G_xvals,G_gs,'r',long_G_xvals,G_es,'b')
subplot(2,1,2)    
plot(long_G_xvals,G_gs2,'r',long_G_xvals,G_es2,'b')

figure(2)
hold off
plot(long_G_xvals,Gaussf(long_G_xvals))

1/(G_vals(1)-G_vals(2))*1000


