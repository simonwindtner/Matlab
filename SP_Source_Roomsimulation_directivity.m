%%--------------------------------------------------------------------------
% short script to simulate sound pressure level trend 
% for given amount of microphones in a room 
% sources are simulated for each x,y and z coordinate in the room 
% Test
%--------------------------------------------------------------------------

clear all;
clc;

%% define number of sources and mic positions
% mic position between 0..1 for x,y and z coordinates 
% all vectors musst be the same size
mic_x = [.2, .5, .8];
mic_y = [.5, .2, .5];
mic_z = [.5, .5, .5];
roomscale = 1;
Z_cs = 1 * 10;

a = 0.63;
b = 0.37;
n_r = 0.5;


%% meshgrid to simulate source positions
%  radius will be referenced to the mid radius of the room r(0.5,0.5,0.5)
%  SPL is simulated for each microphone seperatley 
[X,Y,Z] = meshgrid(0:.1:1 * roomscale);
X1 = squeeze(X(:,:,Z_cs));
Y1 = squeeze(Y(:,:,Z_cs));
SPL = zeros(11,11,11);
SPL(:,:,:) = -100;
n_db = zeros(11,11,11);
n_db(:,:,:) = -100;
niere = zeros(11,11,11);
niere(:,:,:) =  -100;
for Z_cs = 1:10 
    for n = 2:2%length(mic_x)
        dx = X - mic_x(n);  
        dy = Y - mic_y(n);
        dz = Z - mic_z(n);
        r = sqrt(dx.^2 + dy.^2 + dz.^2) + 1;
        r_niere = sqrt(dx.^2 + dy.^2) + 1;
        alpha = acosd(dy ./ sqrt(dx.^2 + dy.^2) );
        for col = 1:11
            for row = 1:11 
                if (alpha(row,col, Z_cs) > 90) || (r_niere(row,col,Z_cs) > ((n_r * (a + b*cosd(alpha(row,col,Z_cs))) + 1)))
                    n_db(row,col, Z_cs) = -3;
                else
                    n_db(row,col, Z_cs) = 0;
                end
            end
        end
        niere = max(n_db,niere);
        surf(X1,Y1,niere(:,:,1));
        S = -20*log10(r / (sqrt(3)/2));
        SPL = max(SPL,S) + niere;
    end
    SPL_print = squeeze(SPL(:,:,Z_cs));
    temp = figure(Z_cs)
    %surf(X1,Y1,SPL_print);
    xlabel('X')
    ylabel('Y')
    caxis([-10 0])
    view(2)
    colorbar
    for n = 2:2%length(mic_x)
        hold on 
        plot(mic_x(n),mic_y(n) ,'r*')
    end   
% for exporting the plots storage path has to be changed
%     title(['SPL trend for Z - level ' num2str(Z_cs/10)])
%     print(temp, '-djpeg' , ['/Users/simonwindtner/Documents/ET-TI/SMA/Ensemble/'...
%                      'Z' num2str(Z_cs)]);
end
%%
a = 0.63;
b = 0.37;
phi = 0:0.15:pi/2;
%phi = 0:0.15:2*pi/4;
theta = 0:0.1:2*pi
r = 0.4 * (a + b*cos(phi));
% [x,y] = pol2cart(theta,rho)
figure(10)
polarplot(phi,r)
hold on
% h = plot(mic_x(1),mic_y(1) ,'r*')
% uistack(h,down);
% surf(X1,Y1,SPL_print);
%%

%%Richtcharakteristik für Niere 

[X,Y,Z]=sphere(10);
[tau,phi,r] = cart2sph(X,Y,Z);                 
gamma = acos(Z);                              
richtniere=abs(0.63+0.37*cos(gamma));         
[x_niere,y_niere,z_niere]=sph2cart(tau,phi,richtniere);
figure;
surf(x_niere,y_niere,z_niere); 
axis equal;
view(-30,10); 
title('Richtcharakteristik Niere');

r_n = sqrt(x_niere.^2 + y_niere.^2 + z_niere.^2);    
