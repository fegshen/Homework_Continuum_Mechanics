clear;clc;close all;
%--- Simulate the propagation of wave
%-- Image size
NX = 512;
NY = 512;

%-- numble of speckles(S)
S = 1200;

%-- speckle size(a)
a = 2;

%-- peak intensity of each speckle
I0 = 1;

%-- wave number and angular frequency
KX = pi/50;
KY = pi/50;
w = 2*pi;

%-- deformation time
d_ts = 0:1.2:6;

%-- deformation amplititude
d_amp = 10;

%main program
rng(1); % random seed
%- generate the speckle randomly
X=NY*rand(S,1);
Y=NY*rand(S,1);
%- wave propgation
%- calculate the deformated speckle image in a fast way
num_img = 0;
N_vector = [KX,KY]/norm([KX,KY]); N_vectorx = N_vector(1);N_vectory = N_vector(2);
for d_t = d_ts
    phase = w*d_t - [X,Y] * [KX;KY]; phase(phase<0) = 0;
    XY_deform = [X,Y]+d_amp*sin(phase).* N_vector;
    X_deform = XY_deform(:,1);
    Y_deform = XY_deform(:,2);
    %- generation of speckles image after deformation
    i = 1:NX;j=1:NY;
    X_cover = erf((i-X_deform(1:end))./a)-erf((i+1-X_deform(1:end))./a);
    Y_cover = erf((j-Y_deform(1:end))./a)-erf((j+1-Y_deform(1:end))./a);
    I_defrom= I0.*pi./4.*a.^2.*(X_cover'*Y_cover);
    A=double(I_defrom);
    
    if ~exist('max_int','var')
        max_int = max(A(:));
    end
    G = A/max_int;G(G>1)=1;
    imwrite(G,[num2str(num_img),'.jpg'])
    
    num_img = num_img +1;
end

%-- Calculate F,C and 2st P-K stree
num_img = 0;
c1 = 0.10445;c2=6.86123;c3=0.001;c4=0.00491;d=2000;rho0 = 1.12*10^3;
[gridX,gridY] = meshgrid(1:NX,1:NY);
for d_t = d_ts
    phase = w*d_t - KX*gridX-KY*gridY; phase(phase<0) = 0;
    %- gradient of displacement
    tmp = cos(phase).*(phase>0);
    Uxx = -KX*d_amp*N_vectorx*tmp;
    Uxy = -KY*d_amp*N_vectorx*tmp;
    Uyx = -KX*d_amp*N_vectory*tmp;
    Uyy = -KY*d_amp*N_vectory*tmp;
    
    %- F
    Fxx = Uxx+1;
    Fxy = Uxy;
    Fyy = Uyy+1;
    Fyx = Uyx;
    
    %- C
    Cxx = Fxx.*Fxx+Fyx.*Fyx;
    Cxy = Fxx.*Fxy+Fyx.*Fyy;
    Cyx = Cxy;
    Cyy = Fxy.*Fxy+Fyy.*Fyy;
    
    %- priciple direction
    Main_strain = N_vectorx^2*Cxx+2*N_vectorx*N_vectory*Cxy+N_vectory^2*Cyy;
    %-- Draw Strain
    figure;
    clims = [0 3.3];
    im_tmp = imagesc(Main_strain,clims);
    colorbar;
    saveas(im_tmp,['strain',num2str(num_img),'.jpg']);
    
    
    J = Fxx.*Fyy-Fxy.*Fyx;
    I_1 = Main_strain + 1 + 1; % three princple stress
    I_4 = Main_strain;
    S = rho0*2/d*(J-1).*J./Main_strain + 2*c1*rho0*(I_1-3).*exp(c2/2*(I_1-3)) + 2*c3*rho0*(I_4-1).*exp(c4*(I_4-1).^2);
    %-- Draw Stress
    figure;
    im_tmp2 = imagesc(S);
    colorbar;
    saveas(im_tmp2,['stress',num2str(num_img),'.jpg']);
    num_img = num_img +1;
end


