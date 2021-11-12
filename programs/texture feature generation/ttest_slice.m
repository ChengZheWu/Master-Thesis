function ttest_slice(filename, modeltype, thick, fs)
fileType = '.data';


%glcm
model_1 = ['COR1';'COR2';'COR3';'COR4';];
%basic
model_2 = ['mea';'UNI';];
%tamura
model_3 = ['coa';'con';];
%glrlm
model_4 = ['GLN1';'GLN2';'GLN3';'GLN4';'RLN1';'RLN2';'RLN3';'RLN4';];
model_5 = ['LGRE1';'LGRE2';'LGRE3';'LGRE4';];
model_6 = ['SRLGE1';'SRLGE2';'SRLGE3';'SRLGE4';];
model_7 = ['RP1';'RP2';'RP3';'RP4';];

model_8 = ['LAP';'LOG';];
model_9 = ['RDF1';'RDF2';'Mean';];
model_10 = ['BRUMM';'BRUM1';'BRUM2';'BSQM2';'SIJM2';'SSQM2';'HistP';];
model_11 = ['ACF';];
model_12 = ['GTFLAP';'AJANE3';'AJANE1';'AJANE2';'AJANE4';'SIJBER';'Energy';];
model_13 = ['Entropy';];
model_14 = ['Variance';'Skewness';'Kurtosis';];

% dataresultFolder = ['t-test\' modeltype '\'];
dataresultFolder = ['t-test/' modeltype '/']; %linux 

for count = 1:14
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
    elseif(count == 12)
        model_name = model_12 ;
    elseif(count == 13)
        model_name = model_13 ;    
    else
        model_name = model_14;
    end

    [tx,ty] = size(model_name);
    for type=1:tx
        model= model_name(type,:);
        cc=0;
        sig_all=0;
        for number = 1:fs
             numberstr = num2str(number,'%03d');    
             sliceNumber = ['_' numberstr]; % _012, _013, _014, _015, _016, _017, _018, _019, _020, _021
%              imageFolder = [modeltype '\' thick '_noise features(slice)\']; %noise features_glcm, noise features1
             imageFolder = [modeltype '/' thick '_noise features(slice)/']; %Linux   

             noiseType = ['number' sliceNumber '_' model];   % number_001img_mean
             a = load ([imageFolder noiseType fileType]);

             tn = number + 1;
             if(tn <= fs)        
                 for tnumber = tn:fs
                     tnumberstr = num2str(tnumber,'%03d');    
                     tsliceNumber = ['_' tnumberstr]; % _012, _013, _014, _015, _016, _017, _018, _019, _020, _021
                     tnoiseType = ['number' tsliceNumber '_' model];   % number_001img_mean   b = load ([imageFolder tnoiseType fileType]);               
                      b = load ([imageFolder tnoiseType fileType]);
                     [h,sig,ci] = ttest2(a,b);
                     cc = cc+1;
                     sig_all = sig_all + sig;

                  end
             end
        end

        cc
        A =model;
        sig_all = sig_all /cc;
        file = [dataresultFolder filename];
        fid = fopen(file,'a+'); 
        fprintf(fid,'%s\t%.6f\n',A,sig_all);
        fclose(fid);

    end
end
