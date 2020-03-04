%%-------------------------------------------------------------------------
% Windtner Simon 2020-02
% release 1.0
% short script to simulate sound pressure level trend 
% for given amount of microphones in a room 
% sources are simulated for each x,y and z coordinate in the room 
% consideration of cardioid characteristics of microphone 
% consideration of direction of microphones
%--------------------------------------------------------------------------

clear all;
clc;


%% define number of sources and mic positions
% mic position between 0..1 for x,y and z coordinates 
% all vectors musst be the same size!

mic_x = [.2, .4, .6, .8];
mic_y = [.5, .2, .2, .5];
mic_z = [.5, .5, .5, .5];

% mic_direction is seen like this
% Y
% 1   --------------------------------
%    |     + <---- 0° ----> -         |
%    |                                |
%    |         \   |   /              |
%    |          \  |  /               |
%    |             x                  |
%    |                                |
% 0   --------------------------------
% X  0                                1  
mic_direction = [30, 0, 0, -30];

% folowing parameters can be changed:
% roomscale -> size of room (standard: 1x1x1)
roomscale = 1;
% gridscale -> points of meshgrid resolutions
gridscale = 0.01;

Z_level = 1 / gridscale;


%% Simulation algorithm
%  meshgrid to simulate source positions
%  SPL is simulated for each microphone seperatley 

% preallocation of meshgrid and vectors
[X,Y,Z] = meshgrid(0:gridscale:1 * roomscale);
SPL = X;
SPL_cardioid_char = X;
SPL_cardioid_char_total = X;
alpha = X;
SPL(:,:,:) = -100;
SPL_cardioid_char(:,:,:) = -100;
SPL_cardioid_char_total(:,:,:) =  -100;


for Z_level = 1:1:Z_level 
    for n = 1:length(mic_x)
        % calculation of distance between microphone and simulated source
        % position -> radius
        dx = X - mic_x(n);  
        dy = Y - mic_y(n);
        dz = Z - mic_z(n);
        r_source = sqrt(dx.^2 + dy.^2 + dz.^2) + 1;
        % calculation of angle between microphone and source position
        for col = 1:length(X)
            for row = 1:length(X)
                if dx(row,col,Z_level) < 0 
                    alpha(row,col,Z_level) = 360 - acosd(dy(row,col,Z_level) / sqrt(dx(row,col,Z_level)^2 + dy(row,col,Z_level)^2));
                else 
                    alpha(row,col,Z_level) =  acosd(dy(row,col,Z_level) / sqrt(dx(row,col,Z_level)^2 + dy(row,col,Z_level)^2));
                end    
            end
        end
        
        % consideration of microphone direction
        alpha = wrapTo360(alpha + mic_direction(n));
        
        % calculation of cardioid characterisation
        for col = 1:length(X)
            for row = 1:length(X)
                if alpha(row,col,Z_level) == 360.0
                    alpha(row,col,Z_level) = 0;
                end    
                SPL_cardioid_char(row,col,Z_level) = 20*log10((1 + cosd(alpha(row,col,Z_level)))/2);
            end
        end
        
        SPL_cardioid_char_total = max(SPL_cardioid_char,SPL_cardioid_char_total);
        SPL_radius = -20*log10(r_source);
        SPL = max(SPL,SPL_radius);
        SPL_total = SPL + SPL_cardioid_char_total;
    end
    % surfaceplot of SPL over X,Y and for current Z - level
    temp = figure(Z_level)
    surf(X(:,:,Z_level),Y(:,:,Z_level),SPL_total(:,:,Z_level));
    xlabel('X')
    ylabel('Y')
    shading interp
    caxis([-10 0])
    view(2)
    colorbar
    
    % plot of microphone position and arrows for direction
    for n = 1:length(mic_x)
        hold on 
        plot(mic_x(n),mic_y(n) ,'r*')
        p0 = [mic_x(n) mic_y(n) 0.0];
        p1 = [(mic_x(n) - tand(mic_direction(n) * 0.25)) (mic_y(n) + .25) 0.0];
        vectarrow(p0,p1);
    end
%% Exporting of Plots    
% for exporting the plots storage path has to be adapted
     title(['SPL trend for Z - level ' num2str(Z_level * gridscale)])
%     print(temp, '-djpeg' , ['/Users/simonwindtner/Documents/ET-TI/SMA/Ensemble/'...
%                      'Z' num2str(Z_cs)]);
end
