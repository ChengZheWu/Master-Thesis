function ttest_pn(filename, modeltype, thickness)
fileType = '.data';

%glcm
% model_1 = ['SDX1';'SDX2';'SDX3';'SDX4';'SDY1';'SDY2';'SDY3';'SDY4';'CON1';'CON2';'CON3';'CON4';'DIS1';'DIS2';'DIS3';'DIS4';'HOM1';'HOM2';'HOM3';'HOM4';
%          'ASM1';'ASM2';'ASM3';'ASM4';'COR1';'COR2';'COR3';'COR4';'ENT1';'ENT2';'ENT3';'ENT4';];
model_1 = ['COR1';'COR2';'COR3';'COR4';];
%basic
model_2 = ['mea';];
%tamura
model_3 = ['coa';'con';];
%glrlm
model_4 = ['GLN1';'GLN2';'GLN3';'GLN4';'RLN1';'RLN2';'RLN3';'RLN4';];
model_5 = ['LGRE1';'LGRE2';'LGRE3';'LGRE4';];
model_6 = ['SRLGE1';'SRLGE2';'SRLGE3';'SRLGE4';];
model_7 = ['RP1';'RP2';'RP3';'RP4';];

model_8 = ['LAP';'LOG';];
model_9 = ['RDF1';'RDF2';];
model_10 = ['BRUMM';'BRUM1';'BRUM2';'BSQM2';'SIJM2';'SSQM2';];
model_11 = ['ACF';];
model_12 = ['GTFLAP';'AJANE3';'AJANE1';'AJANE2';'AJANE4';'SIJBER';];

% dataresultFolder = ['t-test\' modeltype '\'];
dataresultFolder = ['t-test/' modeltype '/']; %Linux

for count = 1:12
    if(count == 1)
        model_name = model_1 ;
    elseif(count == 2)
        model_name = model_2 ;
    elseif(count == 3)
        model_name = model_3 ;
    elseif(count == 4)
        model_name = model_4 ;
    elseif(count == 5)
        model_name = model_5 ;
    elseif(count == 6)
        model_name = model_6 ;
    elseif(count == 7)
        model_name = model_7 ;
    elseif(count == 8)
        model_name = model_8 ;
    elseif(count == 9)
        model_name = model_9 ;
    elseif(count == 10)
        model_name = model_10 ;
    elseif(count == 11)
        model_name = model_11 ;
    else
        model_name = model_12;
    end
    fpn = 9;
    [tx,ty] = size(model_name);
    for type=1:tx
        model= model_name(type,:);
        cc=0;
        sig_all=0;
        for pnnumber = 1:2:fpn
            pnnumberstr = num2str(pnnumber);
            for filenumber = 0:20:40
                filenumberstr = num2str(filenumber);   

%                  imageFolder = [modeltype '\' thickness '_noise features(pn)\']; %noise features_glcm, noise features1 
                 imageFolder = [modeltype '/' thickness '_noise features(pn)/']; %Linux

                 noiseType = ['pn' pnnumberstr '_rf' filenumberstr '_' model];   % pn1_rf0CON, pn1_rf20CON, pn1_rf40CON
                 a = load ([imageFolder noiseType fileType]);

                 noiseType_0= ['pn0_rf0_' model];   % pn10_rf0_CON1, pn1_rf20CON, pn1_rf40CON
                 d = load ([imageFolder noiseType_0 fileType]);

                 [h,sig,ci] = ttest2(a,d);
                 cc = cc+1;
                 sig_all = sig_all + sig;

                 trf = filenumber + 20;
                 if(filenumber <= 40)
                     for pfilenumber = trf:20:40
                         pfilenumberstr = num2str(pfilenumber);   
                         pnoiseType = ['pn' pnnumberstr '_rf' pfilenumberstr '_' model];   % pn1_rf0CON, pn1_rf20CON, pn1_rf40CON
                         c = load ([imageFolder pnoiseType fileType]);
                         [h,sig,ci] = ttest2(a,c);
                         cc = cc+1;
                         sig_all = sig_all + sig;
                     end
                 end

                 tpn = pnnumber + 2;
                 if(tpn <= fpn)        
                     for tpnnumber = tpn:2:fpn;
                         tpnnumberstr = num2str(tpnnumber);
                        for tfilenumber = 0:20:40 
                            tfilenumberstr = num2str(tfilenumber);   
                            tnoiseType = ['pn' tpnnumberstr '_rf' tfilenumberstr '_' model];   % pn1_rf0CON, pn1_rf20CON, pn1_rf40CON
                            b = load ([imageFolder tnoiseType fileType]);               
                            [h,sig,ci] = ttest2(a,b);
                            cc = cc+1;
                            sig_all = sig_all + sig;

                        end
                     end
                 end
            end
        end
        
        cc
        A =model;
        sig_all = sig_all /cc
        file = [dataresultFolder filename];
        fid = fopen(file,'a+'); 
        fprintf(fid,'%s\t%.6f\n',A,sig_all);
        fclose(fid);
        
    end
end
   
