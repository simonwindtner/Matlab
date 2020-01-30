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

%% meshgrid to simulate source positions
%  radius will be referenced to the mid radius of the room r(0.5,0.5,0.5)
%  SPL is simulated for each microphone seperatley 
[X,Y,Z] = meshgrid(0:.01:1 * roomscale);
SPL = X;
SPL = -100;
for Z_cs = 1:10
    for n = 1:length(mic_x)
        dx = X - mic_x(n);  
        dy = Y - mic_y(n);
        dz = Z - mic_z(n);
        r = sqrt(dx.^2 + dy.^2 + dz.^2) + 1;
        S = -20*log10(r / (sqrt(3)/2));
        SPL = max(SPL,S);
    end
    temp = figure(Z_cs)
    surf(X(:,:,Z_cs),Y(:,:,Z_cs),SPL(:,:,Z_cs));
    xlabel('X')
    ylabel('Y')
    caxis([-10 0])
    view(2)
    colorbar
    for n = 1:length(mic_x)
        hold on 
        plot(mic_x(n),mic_y(n) ,'r*')
    end   
% for exporting the plots storage path has to be changed
    title(['SPL trend for Z - level ' num2str(Z_cs/10)])
    %print(temp, '-djpeg' , ['/Users/simonwindtner/Documents/ET-TI/SMA/Ensemble/'...
                    % 'Z' num2str(Z_cs)]);
end
%%









