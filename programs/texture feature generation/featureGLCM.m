function [sdx, sdy, CON, DIS, HOM, ASM, COR, ENT] = featureGLCM(imgNoisy)
% [sigma_x sigma_y CON DIS HOM ASM COR ENT] = featureGLCM(image, offsets, gi):
%   
% @Parameters@
%   image = input array with double, the image which is a difference image
%   of the denoised image and original image
%   offsets =  A p-by-2 array of offsets specifying the distance between the 
%   pixel-of-interest and its neighbor. Each row in the array is a two-element
%   vector, [ROW_OFFSET COL_OFFSET], that specifies the elationship, or 'Offset',
%   between a pair of pixels. 
%
% @Returns@
%   sigma_x = Standard deviation
%   sigma_y = Standard deviation
%   CON = Contrast, one of the GLCM features
%   DIS = Dissimilarity, one of the GLCM features
%   HOM = Homogeneity, one of the GLCM features
%   ASM = Angular Second Moment, one of the GLCM features
%   COR = Correlation, one of the GLCM features
%   ENT = Entropy, one of the GLCM features
%
% Authors: Yu-Ning Chang, Yu-Ju Lin
% Copyright (c) 2014 Herbert H. Chang, Ph.D.
% Computational Biomedical Engineering Laboratory (CBEL)
% Department of Engineering Science and Ocean Engineering
% National Taiwan University, Taipei, Taiwan
% ------------------------------------------

imgNoisy = double(imgNoisy);

% compute GLCMs
offsets = [0 1;-1 1;-1 0;-1 -1];
[GLCMS,SI] = graycomatrix(imgNoisy,'Of',offsets,'NumLevels', 8,'Symmetric', true , 'G', [min(imgNoisy(:)) max(imgNoisy(:))]); 

for gi = 1:4
   % GLCMs of 4 different degree
    glcm = GLCMS(:,:,gi);
    S=size(glcm,1);
    % Normalization
    M = glcm/sum(glcm(:));
    py = sum(M,1);
    px = sum(M,2);
    % GLCM Mean
    for i = 1:S
        ux(i,1) = i*px(i,1);
        uy(1,i) = i*py(1,i);
    end
    ux = sum(ux(:));
    uy = sum(uy(:));   
    % GLCM Variance
    for i = 1:S
        for j= 1:S
            sigma_x(i,j)= ((i-ux)^2)*M(i,j);
            sigma_y(i,j)= ((i-uy)^2)*M(i,j);
        end
    end
    %Standard Deviation
    sdx(gi,1) = sqrt(sum(sum(sigma_x)));
    sdy(gi,1) = sqrt(sum(sum(sigma_y)));

    for i=1:S
        for j=1:S
            %Contrast
            Con(i,j)=((i-j)^2)*M(i,j);
            %Dissimilarity (DIS)
            D(i,j)=abs(i-j)*M(i,j);
            %Homogeneity (HOM)
            H(i,j)=(1/(1+(i-j)^2))*M(i,j);
            % Energy(ASM),Homogeneity (Angular Second Moment) 
            A(i,j)=M(i,j)^2;
            %Entropy
            E(i,j)=M(i,j).*log(M(i,j));  
            %Correlation
            Cor(i,j)=((i*j)*M(i,j)-ux*uy)/(sdx(gi,1)*sdy(gi,1));        
        end
    end

    CON(gi,1) = sum(sum(Con));
    DIS(gi,1) = sum(sum(D));
    HOM(gi,1) = sum(sum(H));
    ASM(gi,1) = sum(sum(A));
    COR(gi,1) = sum(sum(Cor));

    E=M.*log(M);
    E1=~isnan(E);
    E2=E(E1);
    ENT(gi,1)=-sum(sum(E2)); 
end


end