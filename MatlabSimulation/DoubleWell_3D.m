import DVRPackage.*
clc
clear all
close all

plotWhileCalculating = true;
NumPts = 20;
UMatrix=zeros(NumPts);
JMatrix=zeros(NumPts );
OIMatrix=zeros(NumPts );
separation=zeros(NumPts ,1);
potential=zeros(NumPts ,1);
KHZ=UnitsConstants.kHz;
% number of energy eigestates to calculate
Nbands=4;
mass = UnitsConstants.mRb87;
w0 = 707*UnitsConstants.nm;
% Shoudl calculate this based on w0 and lambda...
zR = 2170*UnitsConstants.nm;
%Scaling of second derivative in waist units
Bare_wScal = EnergyFromLengthScale(w0, mass); 
%Scaling of second derivative in Rayleigh-range units
Bare_zScal = EnergyFromLengthScale(zR, mass); 

disp('Starting Calculation')
for separationInc=1:1:NumPts 
    for depthInc=1:1:NumPts 
        disp("Run: " + string(separationInc) + ", " + string(depthInc))
        %% Set trap parameters
        asep = 600*UnitsConstants.nm/w0*((separationInc-1)*0.08+1);
        V0 = 1000*UnitsConstants.kHz*depthInc;
        separation(separationInc)=asep;
        potential(depthInc)=V0;
        %% Not used currently
        % single well trap frequency
        %effOmega=sqrt(UnitsConstants.h*V0*4.0/(mass*w0^2));
        %holen=sqrt(UnitsConstants.hbar/(mass*effOmega));
        %% More setup
        %Double-well Gaussian potential
        %range of DVR spaces (units of w_0 for x,y, units of zR for z)
        ax=2; dx=0.04; [xvals,Nx] = GetGridDx(ax,dx);
        ay=1.5; dy=0.04; [yvals,Ny] = GetGridDx(ay,dy);
        az=1; dz=0.022; [zvals,Nz] = GetGridDx(az,dz);
        inv_z_f=@(z) 1.0/(1.0+z.^2);
        Gaussf=@(x,y,z) -V0*inv_z_f(z)*exp(-2.0*inv_z_f(z)*y.^2)* ...
                        (exp(-2.0*inv_z_f(z)*(x-0.5*asep).^2)+...
                         exp(-2.0*inv_z_f(z)*(x+0.5*asep).^2)); 
        %Double-Gaussian potential with x,y in waist units, z in 
        % Rayleigh-range units
        %% Get even and odd parity harmonic oscillator states/energies
        % 'p' is for even symmetry solutions, so here get the solutions
        % which are even in all three axes. I think this is doing the hard
        % calculations. 
        [ppp_vecs,ppp_vals] = DVR_3D(xvals, yvals, zvals, ...
                                     Bare_wScal, Bare_wScal, Bare_zScal, ...
                                     'p','p','p',Nbands,Gaussf);
        % 'm' is for even symmetry solutions, so get the solutions which
        % are odd in the x axis and even in y and z.
        [mppvecs, mppvals] = DVR_3D(xvals, yvals, zvals, ...
                                    Bare_wScal, Bare_wScal, Bare_zScal, ...
                                    'm','p','p',Nbands,Gaussf);
        % disp('Even Double-well energies')
        % pppvals/UnitsConstants.kHz
        % disp('Odd Double-well energies')
        % mppvals/UnitsConstants.kHz

        %Unpack states from the parity-adapted representation
        % naming here is suggestive of ground vs. excited states, but that
        % seems a bit misleading, G_gs2 I think is higher energy than 
        % G_es, really g = even, e = odd is more straightforward I think.
        G_gs = Unpack3DState(ppp_vecs(1,:),Nx,Ny,Nz,'p','p','p')/sqrt(dx*dy*dz);
        G_es = Unpack3DState(mppvecs(1,:),Nx,Ny,Nz,'m','p','p')/sqrt(dx*dy*dz);
        G_gs2 = Unpack3DState(ppp_vecs(2,:),Nx,Ny,Nz,'p','p','p')/sqrt(dx*dy*dz);
        G_es2 = Unpack3DState(mppvecs(2,:),Nx,Ny,Nz,'m','p','p')/sqrt(dx*dy*dz);
        
        if(false)
            long_G_xvals=dx*[-Nx:Nx]; 
            long_G_yvals=dy*[-Ny:Ny]; 
            long_G_zvals=dz*[-Nz:Nz];

            [Xg,Yg,Zg] = NearestSeparableState(G_gs);
            [Xe,Ye,Ze] = NearestSeparableState(G_es);

            figure
            subplot(3,1,1)
            semilogy(long_G_xvals,abs(Xg),'r')
            hold
            semilogy(long_G_xvals,abs(Xe),'b')
            subplot(3,1,2)
            semilogy(long_G_yvals,abs(Yg),'r')
            hold
            semilogy(long_G_yvals,abs(Ye),'b')
            subplot(3,1,3)
            semilogy(long_G_zvals,abs(Zg),'r')
            hold
            semilogy(long_G_zvals,abs(Ze),'b')
        end

        %Define states to localize (along x)
        psis=cell(1,2);
        psis{1}=G_gs;
        psis{2}=G_es;
        t = cputime;
        [outvecs,S] = LocalizeState(psis,Nx,Ny,Nz,dx,dy,dz);
        wannier=zeros(size(psis{1}));
        for i=1:size(outvecs,1)
            wannier=wannier+outvecs(2,i)*psis{i};
        end

        %Hubbard model
        %Define Hamiltonian in basis of states 1,2
        Hami=diag([ppp_vals(1), mppvals(1)])/KHZ;
        %Rotate Hamiltonian to this basis to obtain (non-interacting) Hubbard model
        HHubb=(outvecs')*Hami*outvecs;
        disp('Hubbard parameters E1, E2, J:')
        disp([ HHubb(1,1) ,HHubb(2,2), HHubb(1,2)])
        JMatrix(separationInc,depthInc)=HHubb(1,2);
        OIMatrix(separationInc,depthInc)=HHubb(2,2)-HHubb(1,1);

        %Interaction energy
        UMatrix(separationInc,depthInc) = InteractionFromLengthScale(w0,w0,zR,mass)...
                       *swaveintegral(wannier,wannier,dx,dy,dz)/KHZ;
        % disp('Hubbard interaction U' )
        % disp(U)

        %Nearest separable
        if (false)
            [X,Y,Z] = NearestSeparableState(wannier);
            figs{1} = figure(1);
            hold on
            subplot(3,1,1);
            plot(long_G_xvals,X,'r')
            hold on
            subplot(3,1,2)    
            plot(long_G_yvals,Y,'r')
            hold on
            subplot(3,1,3)    
            plot(long_G_zvals,Z,'r')
            drawnow;
        end
        %V0
        %asep
        %UMatrix(l,k)
        %JMatrix(l,k)
        if (plotWhileCalculating)
           figure(2);
            imagesc(abs(UMatrix))
            xlabel('depth')
            ylabel('Separation')
            figure(3);
            imagesc(abs(JMatrix))
            xlabel('depth')
            ylabel('Separation')
            %figure(4);
            %imagesc(abs(OIMatrix)) 
            %xlabel('depth')
            %ylabel('Separation')
        end
    end
end

figure(2);
imagesc(abs(UMatrix))
xlabel('depth')
ylabel('Separation')
figure(3);
imagesc(abs(JMatrix))
xlabel('depth')
ylabel('Separation')
%figure(4);
%imagesc(abs(OIMatrix))
%xlabel('depth')
%ylabel('Separation')

%%
save('C:\Users\Mark-Brown\Code\JandUofTweezers\MarksData.mat','JMatrix','UMatrix','OIMatrix','separation','potential')
