function VegGridGeneration(Species, Density, Method, WindSpeed, Resolution)

    % Input varaibles
    %----------------
    % Species = 'AMAM' (for Ammophila arenaria) or 'AMBR' (for Ammphila breviligulata) or 'LEMO' (for Leymus mollis)
    % Density = number of tiller/m²
    % Method = 'Frontal'
    % WindSpeed = value of wind speed in m/s
    % Resolution = 0.05 % Resolution of the grid in meter (here 5cm)

    % Load Data
    %----------
    load(['D:\Documents\Article\AeolisPaper\DataInfo\Scripts\Zarnetske_',Species,'_Metrics.mat']);
    load(['D:\Documents\Article\AeolisPaper\DataInfo\Scripts\FlexureEquationWindTiller',Species,'.mat']);
    load(['D:\Documents\Article\AeolisPaper\DataInfo\VegetationGirds\RandomVegGrid_',Method,'_',Species,num2str(Density),'_Position.mat']);
    load('D:\Documents\Article\AeolisPaper\DataInfo\Scripts\HighResGridZarnetske1mm.mat');
    
    % Creation of the DownScaled resolution Grid
    %-------------------------------------------
    TunelLong = 9.55; % in meter
    TunelWidth = 1;   % in meter
    xLR = (0:Resolution:TunelLong);
    yLR = (0:Resolution:TunelWidth);
    [XGridLR, YGridLR] = meshgrid(0:Resolution:xLR(end),0:Resolution:yLR(end));
    
    % Position of the Veg box in the Grid
    %------------------------------------
    VegStart = 3.38;
    VegEnd = 4.38;
    
    posVegStartLR = find(XGridLR(1,:) == Resolution*round(VegStart/Resolution));
    posVegEndLR   = find(XGridLR(1,:) == Resolution*round(VegEnd/Resolution));
    
    % Computation of the DownScaling factor
    %--------------------------------------
    DownScale = Resolution/(XGridHR(1,2)-XGridHR(1,1));


    % Adjust leaf height dependig on wind speed
    %------------------------------------------
    Flexure = coefa * WindSpeed + coefb * Density + coefc;
    LeafHeightAjusted = LeafHeight * (Flexure/100);

    % Check if the ajusted height is higher to initial height or inf to 0 
    %--------------------------------------------------------------------
    if LeafHeightAjusted > LeafHeight
        LeafHeightAjusted = LeafHeight;
    elseif LeafHeightAjusted < 0
        LeafHeightAjusted = 0;
    end

    % Computation of the surface according to the Frontal method
    %-----------------------------------------------------------
    Surface      = (( StemWidth * StemHeight  ) + ( ( LeafWidth * LeafHeightAjusted) * LeafNumber )) * 1e+6; % mm²
    SurfaceRayon = sqrt(Surface / pi); % mm

    % Computation of the MAX surface according to the Frontal method
    %---------------------------------------------------------------
    SurfaceMax      = (( StemWidth * StemHeight  ) + ( ( LeafWidth * LeafHeight) * LeafNumber )) * 1e+6; % mm²
    SurfaceRayonMax = sqrt(SurfaceMax / pi) ; % mm


    % Filling of the high definition grid
    %------------------------------------
    veg2DHR = zeros(size(XGridHR));

    for v = 1:length(posX)
        i = posX(v,1);
        j = posY(v,1);
        SubX = i-round(SurfaceRayonMax)-1 : i+round(SurfaceRayonMax)+1; SubX(SubX <= 0) = 1; 
        SubY = j-round(SurfaceRayonMax)-1 : j+round(SurfaceRayonMax)+1; SubY(SubY <= 0) = 1;  SubY(SubY > size(XGridHR,1)) = size(XGridHR,1);
        Dist = sqrt((XGridHR(SubY, SubX) - XGridHR(j, i)).^2 + (YGridHR(SubY, SubX) - YGridHR(j, i)).^2);

        % On parcourt tous les indices de la sous-matrice inférieurs à SurfaceRayon
        for ii = 1:length(SubX)
            for jj = 1:length(SubY)
                if Dist(jj,ii) <= SurfaceRayon/1000
                    % On calcule l'indice correspondant dans la matrice veg2DHR
                    idx_i = i - round(SurfaceRayonMax) - 1 + ii;
                    idx_j = j - round(SurfaceRayonMax) - 1 + jj;
                    % On remplit la matrice avec des 1 pour cet indice
                    veg2DHR(idx_j, idx_i) = 1;
                end
            end
        end
    end

        % Filling the DownScaled resolution Grid 
        %---------------------------------------
        nbCell = DownScale * DownScale;
        
        veg2D = zeros(size(XGridLR));
        for j = 1:length(yLR)-1
            for i = posVegStartLR:posVegEndLR-1%posVegStartLR:posVegEndLR-1
                posVeg = find(veg2DHR((j-1)*DownScale+1:j*DownScale,(i-1)*DownScale+1:i*DownScale) == 1);
                veg2D(j,i) = length(posVeg)/nbCell;%((length(posVeg) * 100)/nbCell)/100; % /100 for [0 1]
                clear posVeg
            end
        end
        
    % Save the vegetation grid
    %-------------------------
    if Resolution < 1
        save(['D:\Documents\Article\AeolisPaper\DataInfo\VegetationGirds\RandomVegGrid_Flexure_',Method,'_',Species,num2str(Density),'_',num2str(round(Resolution*100)),'cm.mat'],'veg2D','WindSpeed');
    elseif Resolution >= 1
        save(['D:\Documents\Article\AeolisPaper\DataInfo\VegetationGirds\RandomVegGrid_Flexure_',Method,'_',Species,num2str(Density),'_',num2str(Resolution),'m.mat'],'veg2D','WindSpeed');
    else
        error('.....!!!!! The Saving Resolution Is Wrong !!!!!.....');
    end
end