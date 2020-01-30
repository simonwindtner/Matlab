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
% mic_direction is seen like this
% - ----> 0°C <---- +
mic_direction = [45, 0, -45];


roomscale = 1;
gridscale = 0.01;
Z_cs = 1 * 100;

arrow_x = [.2 .2; .45 .45; .8 .8];
arrow_y = [.5 .8; .2 .5; .5 .8];


%% meshgrid to simulate source positions
%  radius will be referenced to the mid radius of the room r(0.5,0.5,0.5)
%  SPL is simulated for each microphone seperatley 
[X,Y,Z] = meshgrid(0:gridscale:1 * roomscale);
SPL = X;
n_db = X;
niere = X;
SPL(:,:,:) = -100;
n_db(:,:,:) = -100;
niere(:,:,:) =  -100;
for Z_cs = 1:10:Z_cs 
    for n = 1:length(mic_x)
        dx = X - mic_x(n);  
        dy = Y - mic_y(n);
        dz = Z - mic_z(n);
        %dx_rot = dx .* cosd(mic_direction(n)) - dy .* sind(mic_direction(n));
        %dy_rot = dx .* sind(mic_direction(n)) - dy .* cosd(mic_direction(n));
        r_source = sqrt(dx.^2 + dy.^2 + dz.^2) + 1;
        alpha = acosd(dy ./ sqrt(dx.^2 + dy.^2) );%+ mic_direction(n);
        for col = 1:length(X)
            for row = 1:length(X)
                n_db(row,col,Z_cs) = 20*log10((1 + cosd(alpha(row,col,Z_cs)))/2);
            end
        end
        niere = max(n_db,niere);
        %surf(X(:,:,Z_cs),Y(:,:,Z_cs),niere(:,:,Z_cs));
        S = -20*log10(r_source);
        % / (sqrt(3)/2)
        SPL = max(SPL,S);
        SPL_neu = SPL + niere;
    end    
    temp = figure(Z_cs)
    surf(X(:,:,Z_cs),Y(:,:,Z_cs),SPL_neu(:,:,Z_cs));
    xlabel('X')
    ylabel('Y')
    caxis([-10 0])
    view(2)
    colorbar
    for n = 1:length(mic_x)
        hold on 
        plot(mic_x(n),mic_y(n) ,'r*')
        arrow = annotation('arrow',arrow_x(n,1:end), arrow_y(n,1:end));
        arrow.Color = 'red';
        arrow.LineWidth = 0.85; 
    end
% for exporting the plots storage path has to be changed
     title(['SPL trend for Z - level ' num2str(Z_cs)])
%     print(temp, '-djpeg' , ['/Users/simonwindtner/Documents/ET-TI/SMA/Ensemble/'...
%                      'Z' num2str(Z_cs)]);
end
