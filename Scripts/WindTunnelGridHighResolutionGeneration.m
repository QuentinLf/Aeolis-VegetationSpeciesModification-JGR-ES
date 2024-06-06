% Create the wind tunnel grid with a 1mm resolution
Resolution = 0.001; % in meter
TunelLong  = 9.55;  % in meter
TunelWidth = 1;     % in meter
[XGridHR, YGridHR] = meshgrid(0:Resolution:TunelLong,0:Resolution:TunelWidth);
save('D:\Documents\Article\AeolisPaper\DataInfo\Scripts\HighResGridZarnetske1mm.mat','XGridHR','YGridHR');