# Multiview residual selective kernel networks (MRSKNet) for lung nodule classification associated with texture features
# Abstract
Lung cancer is one of the leading causes of death worldwide. This thesis is dedicated to improving the computer-aided diagnosis (CAD) to predict the likelihood of malignant nodules in computed tomography (CT) images. Because the lung nodules, which have various sizes and shapes, are hard to recognize, we present a multiview residual selective kernel network (MRSKNet), which integrates the advantages of ResNet for feature reuse and SKNet for adaptive receptive field (RF) size selection. The MRSKNet is evaluated on Lung Image Database Consortium and Image Database Resource Initiative (LIDC-IDRI), which is a public database of lung CT images with 877 (447 benign and 430 malignant) nodules. Based on ten-fold cross validation, the experiment with original images achieve area under the receiver operating characteristic curve (AUC) of 0.9696, accuracy of 0.9349, sensitivity of 0.9346. Additionally, seven kinds of texture features calculated from the gray-level co-occurrence matrix (GLCM), gray-level run length matrix (GLRLM), and Tamura methods are exploited. We concatenated the original images with the texture features to diversify the input data. Among them, the homogeneity (HOM) improved the performance most. The AUC, accuracy, and sensitivity which were more considerable were increased to 0.9711, 0.9366, and 0.9556, respectively. Finally, we compared with other baseline models and previous works to validate our proposed method. Experimental results indicated that the classification ability of our MRSKNet with HOM outperformed the most state-of-the-art methods.  
<img src="https://github.com/ChengZheWu/Master-Thesis/blob/main/figures/Lung_Nodule_Classification.png">  






[Multiview residual selective kernel networks for lung nodule classification associated with texture features](https://www.airitilibrary.com/Publication/alDetailedMesh1?DocID=U0001-2510202115112000#Summary)
