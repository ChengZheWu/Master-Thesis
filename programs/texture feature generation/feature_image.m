function [mea ent std var] = feature_image(imgNoisy,maxPossible)
%��v���Ƕ���S�ʰ��B��
%mean
m_imgNoisy = imgNoisy/maxPossible;
mea = mean2(m_imgNoisy);
% entropy(I)
ent = entropy(m_imgNoisy);
std = std2(m_imgNoisy);
var = (std).^2;           
          
end