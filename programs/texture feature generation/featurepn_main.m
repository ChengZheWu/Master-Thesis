clear all;
%featurepn_main: t-test by noise level
%

% Some strings for the directories and filenames.
fileType1 = '.dcm';
fileType2 = '.tif';

model = 'msles';        %'msles'
sequence = 't1';         % t1, t2, pd
thickness = ['1mm';'3mm';'5mm';];%'7mm';'9mm';
slice = [150; 60; 36;];% 26; 20
   
[tx,ty] = size(thickness);

for count = 1:3     

        
    thick = thickness(count, :);
    imageFolder = ['/home/greent911/cbel/SBD/' model '/' sequence '/']; %Linux
%     imageFolder = ['/Users/greent/cbel/SBD/' model '/' sequence '/']; %OSX
    
%     imageFolder = [imageFolder thick '\'];
%     dataresultFolder = [model '\' thick '_noise features(pn)\'];    
    imageFolder = [imageFolder thick '/']; %Linux
    dataresultFolder = [model '/' thick '_noise features(pn)/']; %Linux
        
%%
for pnnumber = 0 
    pnnumberstr = num2str(pnnumber);
    for filenumber = 0
        
        filenumberstr = num2str(filenumber);   
        noiseType = ['pn' pnnumberstr '_rf' filenumberstr];   % pn0_rf0
%         folder = [imageFolder noiseType '\'];
        folder = [imageFolder noiseType '/']; %Linux
 
        for number = 1:slice(count)
            numberstr = num2str(number,'%03d');    
            sliceNumber = ['_' numberstr]; % _012, _013, _014, _015, _016, _017, _018, _019, _020, _021 etc
            filenameNoisy = [folder noiseType sliceNumber fileType1];

            % Read the input images.
            noisyInfo = dicominfo((filenameNoisy));
            imgNoisy = dicomread(noisyInfo);

            % Convert the images into double for numerical computation.
            imgNoisy = double(imgNoisy);           
            % The maximum possible intensity of the input image, assuming 8-bit resolution.
            maxPossible = max(max(imgNoisy)); %4095
            % Obtain the size of the input image.
            [height, width] = size(imgNoisy);   

            %features
            %glcm
            [sigma_x, sigma_y, CON, DIS, HOM, ASM, COR, ENT] = featureGLCM(imgNoisy,maxPossible);
            for gi = 1:4
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SDX' num2str(gi) '.data'], sigma_x(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SDY' num2str(gi) '.data'], sigma_y(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_CON' num2str(gi) '.data'], CON(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_DIS' num2str(gi) '.data'], DIS(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_HOM' num2str(gi) '.data'], HOM(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_ASM' num2str(gi) '.data'], ASM(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_COR' num2str(gi) '.data'], COR(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_ENT' num2str(gi) '.data'], ENT(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
            end
            
            %glrlm
            [SRE, LRE, GLN, RLN, RP, LGRE, HGRE, SRLGE, SRHGE, LRLGE, LRHGE] = featureGLRLM(imgNoisy,maxPossible);
            for gi =1:4
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SRE' num2str(gi) '.data'], SRE(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LRE' num2str(gi) '.data'], LRE(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_GLN' num2str(gi) '.data'], GLN(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_RLN' num2str(gi) '.data'], RLN(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_RP' num2str(gi) '.data'], RP(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LGRE' num2str(gi) '.data'], LGRE(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_HGRE' num2str(gi) '.data'], HGRE(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SRLGE' num2str(gi) '.data'],  SRLGE(gi,1),'-append','delimiter','\t','precision','%.9f') ;
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SRHGE' num2str(gi) '.data'], SRHGE(gi,1),'-append','delimiter','\t','precision','%.9f') ;
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LRLGE' num2str(gi) '.data'], LRLGE(gi,1),'-append','delimiter','\t','precision','%.9f') ;
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LRHGE' num2str(gi) '.data'], LRHGE(gi,1),'-append','delimiter','\t','precision','%.9f') ;
            end

            %tamura
            [coa, con, dir] =featureTamura(imgNoisy, maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_coa.data'], coa,'-append','delimiter','\t','precision','%.9f') ; 
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_con.data'], con,'-append','delimiter','\t','precision','%.9f') ; 
%             dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_dir.data'], dir,'-append','delimiter','\t','precision','%.9f') ; 

            %image
            [mea, ent, stdv, var] = feature_image(imgNoisy,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_mea.data'], mea,'-append','delimiter','\t','precision','%.9f') ;
%             dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_ent.data'], ent,'-append','delimiter','\t','precision','%.9f') ;      
%             dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_std.data'], std,'-append','delimiter','\t','precision','%.9f') ;      
%             dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_var.data'], var,'-append','delimiter','\t','precision','%.9f') ;      
            % (1) Gradient Convariance Matrix Eigenvalue 0.23
            % [Gx, Gy] = gradient(imgNoisy);
            % G = [Gx Gy];
            % C = (G*G');
            % EV = eig(C);
            % GCM=sum(EV)
            
            %(2)Laplaican
            M_LAP=[1 -2 1; -2 4 -2; 1 -2 1];
%             M=[13 -2 -7 -2 13;-2 -17 -22 -17 -2;-7 -22 148 -22 -7;-2 -17 -22 -17 -2;13 -2 -7 -2 13];
            M_LOG=[0 0 -1 0 0;0 -1 -2 -1 0;-1 -2 16 -2 -1;0 -1 -2 -1 0;0 0 -1 0 0];
            LAP=sum(sum(abs(conv2(imgNoisy, M_LAP))));
            LAP=LAP*sqrt(0.5*pi)./(6*(width-2)*(height-2));
            LOG=sum(sum(abs(conv2(imgNoisy, M_LOG))));
            LOG=LOG*sqrt(0.5*pi)./(6*(width-2)*(height-2));
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LAP.data'], LAP,'-append','delimiter','\t','precision','%.9f');      
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LOG.data'], LOG,'-append','delimiter','\t','precision','%.9f');      
            
            % (3)diff operator
            [Gx, Gy] = gradient(imgNoisy);
            [Gxx, Gxy] = gradient(Gx);
            [Gyx, Gyy] = gradient(Gy);
            % 1 0.078
            du=(Gx.^2).*Gyy - (2*Gx.*Gy.*Gxy) + (Gy.^2).*Gxx;
            dd=(Gx.^2)+(Gy.^2);
            d=du./dd;
            d(isnan(d))=0;
            RDF1=sum(sum(abs(d)))/(width*height);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_RDF1.data'], RDF1,'-append','delimiter','\t','precision','%.9f');      
            % 3 0.089
            H = fspecial('average', [3 3]);
            aGx2 = imfilter((Gx.^2), H);
            aGy2 = imfilter((Gy.^2), H);
            a2Gxy = imfilter((Gx.*Gy), H);
            a2Gxy = a2Gxy.^2;
            du = aGx2.*aGy2 - a2Gxy;
            dd = aGx2 + aGy2;
            d=du./dd;
            d(isnan(d))=0;
            RDF2=sum(sum(abs(d)))/(width*height);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_RDF2.data'], RDF2,'-append','delimiter','\t','precision','%.9f');                 

            % (4)
            AJANE3=AjaNE(imgNoisy,3,maxPossible,3);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_AJANE3.data'], AJANE3,'-append','delimiter','\t','precision','%.9f');                             
            % (5)
            BRUMM=brummer(imgNoisy,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_BRUMM.data'], BRUMM,'-append','delimiter','\t','precision','%.9f');                             
            % (6)
            BRUM1=noise_M1(imgNoisy,3,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_BRUM1.data'], BRUM1,'-append','delimiter','\t','precision','%.9f');                             
            BRUM2=noise_M2(imgNoisy,3,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_BRUM2.data'], BRUM2,'-append','delimiter','\t','precision','%.9f');                             
            % (7) 
            BSQM2=noise_SQM2(imgNoisy,3,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_BSQM2.data'], BSQM2,'-append','delimiter','\t','precision','%.9f');                             
            % (8) sij 
            SIJBER=sijbers(imgNoisy,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SIJBER.data'], SIJBER,'-append','delimiter','\t','precision','%.9f');                             
            % (9) sijM2 
            SIJM2=sijbersM2(imgNoisy,maxPossible,3);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SIJM2.data'], SIJM2,'-append','delimiter','\t','precision','%.9f');                             
            SSQM2=sijbersSQM2(imgNoisy,maxPossible,3);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SSQM2.data'], SSQM2,'-append','delimiter','\t','precision','%.9f');                             
            % (10) 
%             GCM=Nsigmaest1(imgNoisy)
            % (11) eneric Transfer Function based Technique forEstimating
            M=[1 -2 1; -2 4 -2; 1 -2 1];
            Nlap=abs(conv2(imgNoisy, M));
            gx=[1 2 1;0 0 0;-1 -2 -1];
            gy=gx';
            Ngx=abs(conv2(imgNoisy, gx));
            Ngy=abs(conv2(imgNoisy, gy));
            Ngx=0.1*Ngx;
            Ngy=0.1*Ngy;
            Thx=max(max(Ngx));
            Thy=max(max(Ngy));
            Nn=Nlap;
            Nn(Nn>Thx)=Thx;
            Nn(Nn>Thy)=Thy;
            GTFLAP=sum(sum(Nn))/(width*height);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_GTFLAP.data'], GTFLAP,'-append','delimiter','\t','precision','%.9f');                             
            % (12) 
            imgNoisy2 = imgNoisy.*imgNoisy;
            r1100 = sum(imgNoisy2(:));
            r1100 = r1100/(height*width);
            temp = imgNoisy(:,2:end);
            temp(:,end+1)=0;
            tempm = imgNoisy.*temp;
            r1110 = sum(tempm(:));
            r1110 = r1110/(height*width);
            temp = temp(:,2:end);
            temp(:,end+1)=0;
            tempm = imgNoisy.*temp;
            r1120 = sum(tempm(:));
            r1120 = r1120/(height*width);
            r11nf = 0.5*(3*r1110-r1120);
            r11nf = r1110;
            mea = mean2(imgNoisy);
            mea2 = mea*mea;
            ACF = (r11nf - mea2)/(r1100-r11nf);
            
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_ACF.data'], ACF,'-append','delimiter','\t','precision','%.9f') ;      
            AJANE1=AjaNE(imgNoisy,3,maxPossible,1);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_AJANE1.data'], AJANE1,'-append','delimiter','\t','precision','%.9f');                             
            AJANE2=AjaNE(imgNoisy,3,maxPossible,2);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_AJANE2.data'], AJANE2,'-append','delimiter','\t','precision','%.9f');                             
            AJANE4=AjaNE(imgNoisy,3,maxPossible,4);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_AJANE4.data'], AJANE4,'-append','delimiter','\t','precision','%.9f');                             
        end
    end
end


%%
for pnnumber = 1:2:9
    pnnumberstr = num2str(pnnumber);
    for filenumber = 0:20:40
        filenumberstr = num2str(filenumber);  
        noiseType = ['pn' pnnumberstr '_rf' filenumberstr];   % pn1_rf0, pn1_rf20, pn1_rf40
%         folder = [imageFolder noiseType '\'];
        folder = [imageFolder noiseType '/']; %Linux

        for number = 1:slice(count)
            numberstr = num2str(number,'%03d');    
            sliceNumber = ['_' numberstr]; % _012, _013, _014, _015, _016, _017, _018, _019, _020, _021 etc
            filenameNoisy = [folder noiseType sliceNumber fileType1];

            % Read the input images.
            noisyInfo = dicominfo((filenameNoisy));
            imgNoisy = dicomread(noisyInfo);

            % Convert the images into double for numerical computation.
            imgNoisy = double(imgNoisy);           
            % The maximum possible intensity of the input image, assuming 8-bit resolution.
            maxPossible = max(max(imgNoisy)); %4095
            % Obtain the size of the input image.
            [height, width] = size(imgNoisy);   

            %features
            %glcm
            [sigma_x, sigma_y, CON, DIS, HOM, ASM, COR, ENT] = featureGLCM(imgNoisy,maxPossible);
            for gi = 1:4
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SDX' num2str(gi) '.data'], sigma_x(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SDY' num2str(gi) '.data'], sigma_y(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_CON' num2str(gi) '.data'], CON(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_DIS' num2str(gi) '.data'], DIS(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_HOM' num2str(gi) '.data'], HOM(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_ASM' num2str(gi) '.data'], ASM(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_COR' num2str(gi) '.data'], COR(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_ENT' num2str(gi) '.data'], ENT(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
            end

            %glrlm
            [SRE, LRE, GLN, RLN, RP, LGRE, HGRE, SRLGE, SRHGE, LRLGE, LRHGE] = featureGLRLM(imgNoisy,maxPossible);
            for gi =1:4
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SRE' num2str(gi) '.data'], SRE(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LRE' num2str(gi) '.data'], LRE(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_GLN' num2str(gi) '.data'], GLN(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_RLN' num2str(gi) '.data'], RLN(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_RP' num2str(gi) '.data'], RP(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LGRE' num2str(gi) '.data'], LGRE(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_HGRE' num2str(gi) '.data'], HGRE(gi,1),'-append','delimiter','\t','precision','%.9f') ; 
                dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SRLGE' num2str(gi) '.data'],  SRLGE(gi,1),'-append','delimiter','\t','precision','%.9f') ;
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SRHGE' num2str(gi) '.data'], SRHGE(gi,1),'-append','delimiter','\t','precision','%.9f') ;
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LRLGE' num2str(gi) '.data'], LRLGE(gi,1),'-append','delimiter','\t','precision','%.9f') ;
%                 dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LRHGE' num2str(gi) '.data'], LRHGE(gi,1),'-append','delimiter','\t','precision','%.9f') ;
            end
            
            %tamura
            [coa, con, dir] =featureTamura(imgNoisy, maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_coa.data'], coa,'-append','delimiter','\t','precision','%.9f') ; 
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_con.data'], con,'-append','delimiter','\t','precision','%.9f') ; 
%             dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_dir.data'], dir,'-append','delimiter','\t','precision','%.9f') ; 
             
            %image
            [mea, ent, stdv, var] = feature_image(imgNoisy,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_mea.data'], mea,'-append','delimiter','\t','precision','%.9f') ;
%             dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_ent.data'], ent,'-append','delimiter','\t','precision','%.9f') ;      
%             dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_std.data'], std,'-append','delimiter','\t','precision','%.9f') ;      
%             dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_var.data'], var,'-append','delimiter','\t','precision','%.9f') ;      
            % (1) Gradient Convariance Matrix Eigenvalue
            % [Gx, Gy] = gradient(imgNoisy);
            % G = [Gx Gy];
            % C = (G*G');
            % EV = eig(C);
            % GCM=sum(EV)
            
            % (2)Laplaican
            M_LAP=[1 -2 1; -2 4 -2; 1 -2 1];
%             M=[13 -2 -7 -2 13;-2 -17 -22 -17 -2;-7 -22 148 -22 -7;-2 -17 -22 -17 -2;13 -2 -7 -2 13];
            M_LOG=[0 0 -1 0 0;0 -1 -2 -1 0;-1 -2 16 -2 -1;0 -1 -2 -1 0;0 0 -1 0 0];
            LAP=sum(sum(abs(conv2(imgNoisy, M_LAP))));
            LAP=LAP*sqrt(0.5*pi)./(6*(width-2)*(height-2));
            LOG=sum(sum(abs(conv2(imgNoisy, M_LOG))));
            LOG=LOG*sqrt(0.5*pi)./(6*(width-2)*(height-2));
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LAP.data'], LAP,'-append','delimiter','\t','precision','%.9f');      
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_LOG.data'], LOG,'-append','delimiter','\t','precision','%.9f');      
            
            % (3)diff operator
            [Gx, Gy] = gradient(imgNoisy);
            [Gxx, Gxy] = gradient(Gx);
            [Gyx, Gyy] = gradient(Gy);
            % 1 
            du=(Gx.^2).*Gyy - (2*Gx.*Gy.*Gxy) + (Gy.^2).*Gxx;
            dd=(Gx.^2)+(Gy.^2);
            d=du./dd;
            d(isnan(d))=0;
            RDF1=sum(sum(abs(d)))/(width*height);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_RDF1.data'], RDF1,'-append','delimiter','\t','precision','%.9f');      
       
            % 3
            H = fspecial('average', [3 3]);
            aGx2 = imfilter((Gx.^2), H);
            aGy2 = imfilter((Gy.^2), H);
            a2Gxy = imfilter((Gx.*Gy), H);
            a2Gxy = a2Gxy.^2;
            du = aGx2.*aGy2 - a2Gxy;
            dd = aGx2 + aGy2;
            d=du./dd;
            d(isnan(d))=0;
            RDF2=sum(sum(abs(d)))/(width*height);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_RDF2.data'], RDF2,'-append','delimiter','\t','precision','%.9f');                 

            % (4) 
            AJANE3=AjaNE(imgNoisy,3,maxPossible,3);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_AJANE3.data'], AJANE3,'-append','delimiter','\t','precision','%.9f');                             
            % (5)
            BRUMM=brummer(imgNoisy,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_BRUMM.data'], BRUMM,'-append','delimiter','\t','precision','%.9f');                             
            % (6)
            BRUM1=noise_M1(imgNoisy,3,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_BRUM1.data'], BRUM1,'-append','delimiter','\t','precision','%.9f');                             
            BRUM2=noise_M2(imgNoisy,3,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_BRUM2.data'], BRUM2,'-append','delimiter','\t','precision','%.9f');                             
            % (7) 
            BSQM2=noise_SQM2(imgNoisy,3,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_BSQM2.data'], BSQM2,'-append','delimiter','\t','precision','%.9f');                             
            % (8) sij
            SIJBER=sijbers(imgNoisy,maxPossible);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SIJBER.data'], SIJBER,'-append','delimiter','\t','precision','%.9f');                             
            % (9) sijM2 
            SIJM2=sijbersM2(imgNoisy,maxPossible,3);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SIJM2.data'], SIJM2,'-append','delimiter','\t','precision','%.9f');                             
            SSQM2=sijbersSQM2(imgNoisy,maxPossible,3);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_SSQM2.data'], SSQM2,'-append','delimiter','\t','precision','%.9f');                             
            % (10)
%             GCM=Nsigmaest1(imgNoisy)
            % (11) eneric Transfer Function based Technique forEstimating
            % Noise from Image 
            M=[1 -2 1; -2 4 -2; 1 -2 1];
            Nlap=abs(conv2(imgNoisy, M));
            gx=[1 2 1;0 0 0;-1 -2 -1];
            gy=gx';
            Ngx=abs(conv2(imgNoisy, gx));
            Ngy=abs(conv2(imgNoisy, gy));
            Ngx=0.1*Ngx;
            Ngy=0.1*Ngy;
            Thx=max(max(Ngx));
            Thy=max(max(Ngy));
            Nn=Nlap;
            Nn(Nn>Thx)=Thx;
            Nn(Nn>Thy)=Thy;
            %GTFLAP = median(std(im2col(Nn,[11 11])));
            GTFLAP=sum(sum(Nn))/(width*height);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_GTFLAP.data'], GTFLAP,'-append','delimiter','\t','precision','%.9f'); 
            % (12)
            imgNoisy2 = imgNoisy.*imgNoisy;
            r1100 = sum(imgNoisy2(:));
            r1100 = r1100/(height*width);
            temp = imgNoisy(:,2:end);
            temp(:,end+1)=0;
            tempm = imgNoisy.*temp;
            r1110 = sum(tempm(:));
            r1110 = r1110/(height*width);
            temp = temp(:,2:end);
            temp(:,end+1)=0;
            tempm = imgNoisy.*temp;
            r1120 = sum(tempm(:));
            r1120 = r1120/(height*width);
            r11nf = 0.5*(3*r1110-r1120);
            r11nf = r1110;
            mea = mean2(imgNoisy);
            mea2 = mea*mea;
            ACF = (r11nf - mea2)/(r1100-r11nf);
            
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_ACF.data'], ACF,'-append','delimiter','\t','precision','%.9f') ;      
            AJANE1=AjaNE(imgNoisy,3,maxPossible,1);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_AJANE1.data'], AJANE1,'-append','delimiter','\t','precision','%.9f');                             
            AJANE2=AjaNE(imgNoisy,3,maxPossible,2);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_AJANE2.data'], AJANE2,'-append','delimiter','\t','precision','%.9f');                             
            AJANE4=AjaNE(imgNoisy,3,maxPossible,4);
            dlmwrite ([  dataresultFolder 'pn' pnnumberstr '_rf' filenumberstr '_AJANE4.data'], AJANE4,'-append','delimiter','\t','precision','%.9f');                             
        end
    end
end

ttest_fileName = ['ttest_' thick '_' model '(pn).txt'];
ttest_pn(ttest_fileName, model, thick);

end
