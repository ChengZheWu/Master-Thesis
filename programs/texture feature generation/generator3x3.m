clear ;
CoreNum = 4; % 設定CPU core數量
if isempty(gcp('nocreate'))
    parpool(CoreNum);
end

% get the texture images
imgsPath = dir('./images/*.png');
len = length(imgsPath);
height = 32;
width = 32;
ENTimgs = zeros(len, height, width);
HOMimgs = zeros(len, height, width);
GLNimgs = zeros(len, height, width);
RLNImgs = zeros(len, height, width);
RPimgs = zeros(len, height, width);
SREimgs = zeros(len, height, width);
COAimgs = zeros(len, height, width);
f = waitbar(0, 'please wait');
for  i = 1:len
    str=['run...',num2str((i/len)*100),'%'];
    waitbar(i/len, f, str);
    imgpath = strcat(imgsPath(i).folder, '\', imgsPath(i).name);
    img = imread(imgpath);
%     [height, width] = size(img);
    N = 1;
    % zero padding
    imgpad = padarray(img, [N N], 0);
    % Span of the filter.
    span = 2 * N;
    parfor x = 1:height
        for y = 1:width
            % Extract pixel with its neighbors for computation.
            kernel = imgpad(x:x+span, y:y+span);
            [sigma_x, sigma_y, CON, DIS, HOM, ASM, COR, ENT] = featureGLCM(kernel);
            [SRE, LRE, GLN, RLN, RP, LGRE, HGRE, SRLGE, SRHGE, LRLGE, LRHGE] = featureGLRLM(kernel);
            [coa, con, dir] =featureTamura(kernel);
            if img(x, y) == 0
                ENTimgs(i, x, y) = 0;
                HOMimgs(i, x, y) = 0;
                GLNimgs(i, x, y) = 0;
                RLNImgs(i, x, y) = 0;
                RPimgs(i, x, y) = 0;
                SREimgs(i, x, y) = 0;
                COAimgs(i, x, y) = 0;
            else
                ENTimgs(i, x, y) = ENT(1);
                HOMimgs(i, x, y) = HOM(1);
                GLNimgs(i, x, y) = GLN(1);
                RLNImgs(i, x, y) = RLN(1);
                RPimgs(i, x, y) = RP(1);
                SREimgs(i, x, y) = SRE(1);
                COAimgs(i, x, y) = coa;
            end
        end
    end
end
delete(f);

% normalize 0~1
ENTimgs = ENTimgs / max(ENTimgs, [], 'all');
HOMimgs = HOMimgs / max(HOMimgs, [], 'all');
GLNimgs = GLNimgs / max(GLNimgs, [], 'all');
RLNImgs = RLNImgs / max(RLNImgs, [], 'all');
RPimgs = RPimgs / max(RPimgs, [], 'all');
SREimgs = SREimgs / max(SREimgs, [], 'all');
COAimgs = COAimgs / max(COAimgs, [], 'all');

% save images
f = waitbar(0, 'please wait');
for i = 1:len
    str=['run...',num2str((i/len)*100),'%'];
    waitbar(i/len, f, str);
    imgName = imgsPath(i).name;
    ENTimg = reshape(ENTimgs(i, :, :), [height, width]);
    imwrite(ENTimg, strcat('C:\Users\Zhe\Desktop\TextureComputation\Matlab\ENT\', imgName));
    HOMimg = reshape(HOMimgs(i, :, :), [height, width]);
    imwrite(HOMimg, strcat('C:\Users\Zhe\Desktop\TextureComputation\Matlab\HOM\', imgName));
    GLNimg = reshape(GLNimgs(i, :, :), [height, width]);
    imwrite(GLNimg, strcat('C:\Users\Zhe\Desktop\TextureComputation\Matlab\GLN\', imgName));
    RLNImg = reshape(RLNImgs(i, :, :), [height, width]);
    imwrite(RLNImg, strcat('C:\Users\Zhe\Desktop\TextureComputation\Matlab\RLN\', imgName));
    RPimg = reshape(RPimgs(i, :, :), [height, width]);
    imwrite(RPimg, strcat('C:\Users\Zhe\Desktop\TextureComputation\Matlab\RP\', imgName));
    SREimg = reshape(SREimgs(i, :, :), [height, width]);
    imwrite(SREimg, strcat('C:\Users\Zhe\Desktop\TextureComputation\Matlab\SRE\', imgName));
    COAimg = reshape(COAimgs(i, :, :), [height, width]);
    imwrite(COAimg, strcat('C:\Users\Zhe\Desktop\TextureComputation\Matlab\COA\', imgName));
end
delete(f);