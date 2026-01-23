N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);
Lz = 4000;

im = InternalModesWKBSpectral(N2=N2,zIn=[-Lz 0],latitude=31,rho0=1025);

% profile on; [F,G,h] = im.modesAtQuadraturePoints(nPoints=64,k=2*pi/1e3); profile viewer;