# IACV Project Politecnico di Milano 2022-2023
This project aims to analyse the geometric characteristics of clips inside the KT3DMoSeg dataset [[1]](#1), built upon the KITTI benchmark, that can be found in the [Data](./Data/) folder.
KT3DMoSeg dataset was created by manually selecting 22 sequences and labelling each individual foreground object. In the dataset are present sequences with more significant camera translation so camera mounted on moving cars are preferred.
Clips with more than 3 motions are also chosen, as long as these moving objects contain enough features for forming motion hypotheses. 22 short clips, each with 10-20 frames, are chosen for evaluation. We extract dense trajectories from each sequence using [[2]](#2) and prune out trajectories shorter than 5 frames. 

## Project structure
The original repository was created by Xun Xu, Loong-Fah Cheong and Zhuwen Li. Note that our work for this project is almost entirely in the [Visualize](./Visualize/) folder, and other minor changes for the correct introduction of the additional geometric models Fundamental Translational, Fundamental Affine and SubsetOnlyHF and for the visualization of the kmeans results with the correct alpha parameters.

* **Data**: the KT3DMoSeg dataset. One `.mat` file per clip describing all its characteristics. 
* **OriginalSequence**: the KT3DMoSeg dataset in images format.
* **RepresentativeResults**: here you can find some of the results to grasp the utility of our work. All the results can be found in the relative folder as described below.
* **Results**: this folder is used to save the results obtained from the following operations:
  - Sample hypotheses (Affine, Homography and Fundamental Matrix)
  - Compute the Ordered Residual Kernel (ORK) 
  - Motion Segmentation results using single-view and multi-view spectral clustering (CoReg, Subset, SubsetOnlyHF). 
* **Tools**: libraries and support code, including:
  - MATLAB Functions for Multiple View Geometry
  - Accelerated Hypothesis Generation for Multi-Structure Robust Fitting
  - Sparse Subspace Clustering (SSC)
* **script**: the main folder containing the following subfolders:
  - Hypo: Scripts to sample hypotheses from each model (Affine, Homography Fundamental, FundamentalT and FundamentalA).
  - Kernel: Compute Ordered Residual Kernel (ORK).
  - MoSeg: Run the single-view (script_X_MoSeg_RandSamp.m) or multi-view motion segmentation (script_AHF_X_RandSamp.m, where X is KerAdd/CoReg/Subset/SubsetOnlyHF). 
  - CheckPerf: check the motion segmentation results. Due to randomness in hyphoteses sampling, the performance may vary in each run.
  - Visualize: visualize and save all the elements for the analysis of the clips.
## Contributions
In this part we highlight specifically all the code contributions we developed in the repository:
* [**Hypo**](./script/Hypo),[**Kernel**](./script/Kernel),[**MoSeg**](./script/MoSeg): this is the part of the algorithm of the paper [[1]](#1), we added the MATLAB files for the computation of the Hypothesis, Kernel and Motion Segmentation for the models SubsetOnlyHF, FundamentalT and FundamentalA.
  - SubsetOnlyHF: `script_AHF_SubsetOnlyFH_RandSamp.m`
  - FundamentalT:  `script_FundamentalT_Hypo_RandSamp.m`, `script_FundamentalT_ORK_RandSamp.m`, `script_FundamentalT_MoSeg_RandSamp.m`.
  - FundamentalA:  `script_FundamentalA_Hypo_RandSamp.m`, `script_FundamentalA_ORK_RandSamp.m`, `script_FundamentalA_MoSeg_RandSamp.m`.
* [**Results**](./Results): In [MoSeg](./Results/MoSeg) we stored an archive with 8 full experiments computed following the indication of the paper. This data will be used for all the anlysis.
* [**Tools**](./Tools): we adapted the model specific functions for fitting, computing the residuals and handling the degeneracies of the particular models in [Model Specific](./Tools/multigs/model_specific)
  - FundamentalT:  `fundamentalT_fit.m`, `fundamentalT_res.m`, `fundamentalT_degen.m`.
  - FundamentalA:  `fundamentalT_fit.m`, `fundamentalT_res.m`, `fundamentalT_degen.m`.
  - We also adapted the functions for Homography and Fundamental to our needs.
* [**CheckPerf**](./script/CheckPerf): our contribution in this folder is related to the tabular representations of the results and the search for alpha parameters:
  - `script_CheckCommonAlpha.m` : Find the optimal alpha for each model considering all the experiments in separately and combine the  segmentation results together in a table
  - `script_CheckCommonAlphaAverage.m` : Find the optimal alpha for each model considering all the experiments averaging the results and combine the segmentation results together in a table.
  - `script_CheckPerf_RandSamp.m` : Find the optimal alpha for each experiment (alphas are not constant among models) and combine the segmentation results together in a table.
  - You can also find here the results table formatted in pretty way (`datiBestAlpha.xls`, `datiSameAlpha.xls` and `datiSameAlphaAverage.xls`)
* [**Visualize**](./script/Visualize): in this folder is present the majority of our work conducted for the project. We created different `.m` scripts for the visualization of metrics such as epipolar lines, silhouette scores and evaluation clustering results. in particular:
  - [Epipolar Lines](./script/Visualize/EpipolarLines): Here you can find:
    - `VisualizeEpipolarLines.m` : Visualize and save in [EpipolarLinesVideosOfBestF](./script/Visualize/EpipolarLines/EpipolarLinesVideosOfBestF) the epipolar lines computed with the minimum number of points.
    - `VisualizeEpipolarLinesWithAllPoints.m` : Visualize and save in [EpipolarLinesVideosWithAllPoints](./script/Visualize/EpipolarLines/EpipolarLinesVideosWithAllPoints) the epipolar lines computed with the maximum number of points.
    - `VisualizeEpipolarLinesWithMotionPoints.m` : Visualize and save in [EpipolarLinesVideosWithMotionPoints](./script/Visualize/EpipolarLines/EpipolarLinesVideosWithMotionPoints) the epipolar lines computed with the maximum number of points divided per motion.
  - [**HyperparameterTuning**](./script/Visualize/HyperparameterTuning): Here you can find:
    - `ModelSelectionHypertuning.m` : Perform the process of parameter tuning for lambda1, lambda2 and lambda3 of the GRIC. Results are saved in  `BestLambdasNormalized.m` and `BestLambdas.m`
  - [**Kmeans and Embeddings**](./script/Visualize/Kmeans%20and%20Embeddings) : Here you can find:
    - `VisualizeKmeansInteractive.m` : Perform kmeans with the choices of the model and experiment, outputs the video in [SegmentationVideos](./script/Visualize/Kmeans%20and%20Embeddings/SegmentationVideos).
    - `VisualizeSilhouettePlots.m` : Plot all the silhouette plots with respect of the embeddings for all the models.
  - [**Model Silhouette**](./script/Visualize/Model%20Silhouettes) : Here you can find:
    - `VisualizeModelsSilhouettePlots.m` /  `VisualizeModelsSilhouettePlotsNormalized.m` :  Save for each geometric model, foe each clip, for each couple of frames the silhouette plots derived from the geometric models. Results saved in [ModelSilouettes](./script/Visualize/Model%20Silhouettes/ModelSilouettes).
    - `VisualizeVideoModelSilhouettes.m` : Combine the results of the previous files in a video in [ModelSilhouettesNormalized](./script/Visualize/Model%20Silhouettes/ModelSilouettesNormalized).
  - [**Model Selection**](./script/Visualize/Model%20Selection) : Here you can find:
    - `ModelSelection.m` / `ModelSelectionNormalized.m` : Perform the process of model selection wrt various selection criteria and output a file `seq_best_model.mat` with the models for each motion. Results are stored in [ModelSelectionResults](./script/Visualize/Model%20Selection/ModelSelectionResults).
    - `ModelSelectionSilhouettePlotsNormalized.m` : Starting from the output of the selection process compute the silhouette plots and save them in [ModelSelectionSilhouettesNormalized](./script/Visualize/Model%20Selection/ModelSelectionSilhouettesNormalized).
    - `VideoModelSelction.m` : Combine the results of the previous file in a video in [ModelSelectionSilhouettesNormalized](./script/Visualize/Model%20Selection/ModelSelectionSilhouettesNormalized).

# Usage
1. Refer to the Original README section to produce the experiment results.
2. Store the results in Archivio.
3. Visualize the results from [CheckPerf](./script/CheckPerf) folder.
4. There is not a specific order in the analysis we offer in [Visualize](./script/Visualize) folder. In general you can follow the order presented above.
# References
<a id="1">[1]</a> 
X. Xu, L.F. Cheong, and Z. Li. Motion segmentation by exploiting complementary geometric models, In CVPR 2018.

<a id="2">[2]</a> 
N. Sundaram, T. Brox, and K. Keutzer. Dense point trajectories by GPU-accelerated large displacement optical flow. In ECCV, 2010.


# Original README
Here we report the readme of the original repository.
### Dataset and code released for
# Xun Xu, Loong-Fah Cheong and Zhuwen Li, Motion segmentation by exploiting complementary geometric models, CVPR 2018
#
# By: Xun Xu, Sep 2018
# Contact: alex.xun.xu@gmail.com

Dataset:

The KT3DMoSeg dataset [1] was created upon the KITTI benchmark [2] by manually selecting 22 sequences and labelling each individual foreground object. We select sequence with more significant camera translation so camera mounted on moving cars are preferred. We are interested in the interplay of multiple motions, so clips with more than 3 motions are also chosen, as long as these moving objects contain enough features for forming motion hypotheses. 22 short clips, each with 10-20 frames, are chosen for evaluation. We extract dense trajectories from each sequence using [3] and prune out trajectories shorter than 5 frames.

Raw Sequence:

We provide the raw sequence under ./OriginalSequence/

Data Format:

All meta data are saved as Matlab mat files under ./Data/ . All information are fields of struct varaible "Data". All dense trajectories are indexed by "Data.yAll" and all sparse trajectories are indexed by "Data.ySparse". The ground-truth label for all sparse trajectories is indexed by "Data.GtLabel". Motion segmentation is only evaluated on sparse trajectories. The visible mask is indexed by "Data.visibleSparse/visibleAll". The visibleSparse/visibleAll is a boolean matrix with NxL dimension where N is the number of trajectories and L is the total number of frames. Additional original video sequences can be accessed at https://www.dropbox.com/sh/4u1p3xwe6v48kww/AACLNzNSWg_YYZ6dsIXtLRoTa?dl=0

Demo Code:

A demo code to sample hypotheses (Affine, Homography and Fundamental Matrix), Compute Ordered Residual Kernel (ORK) [4] and single-view or multi-view spectral clustering is provided. The code is based on sampling hypotheses from all sparse trajectories and test/evaluate on sparsely labelled trajectories.
To run this demo, follow the steps:
(1) Run hypotheses generation code from each model under ./script/Hypo/
(2) Run ORK kernel computing code under ./script/Kernel/
(3) Run single-view (script_X_MoSeg_RandSamp.m, where X is Affine/Homography/Fundamental) or multi-view (script_AHF_X_RandSamp.m, where X is KerAdd/CoReg/Subset) motion segmentation under ./script/MoSeg/
(4) Check motion segmentation results by running the code under ./script/CheckPerf/
Due to randomness in hypotheses sampling, the performance may vary in each run. An optimal alpha is around 10.


Please cite [1] if you wish to use this dataset and demo code.

The libraries and code to support this demo include:

MATLAB Functions for Multiple View Geometry
David Capel, Andrew Fitzgibbon, Peter Kovesi, Tomas Werner, Yoni Wexler, and Andrew Zisserman
http://www.robots.ox.ac.uk/~vgg/hzbook/code/

Accelerated Hypothesis Generation for Multi-Structure Robust Fitting
T.-J. Chin
https://cs.adelaide.edu.au/~tjchin/doku.php?id=publications_code

Sparse Subspace Clustering (SSC)
Ehsan Elhamifar
http://www.ccs.neu.edu/home/eelhami/codes.htm

[1] X. Xu, L.F. Cheong, and Z. Li. Motion segmentation by exploiting 
complementary geometric models, In CVPR 2018.
[2] A. Geiger, P. Lenz, C. Stiller, and R. Urtasun. Vision meets robotics: The kitti dataset. International Journal of Robotics Research, 2013.
[3] N. Sundaram, T. Brox, and K. Keutzer. Dense point trajectories by GPU-accelerated large displacement optical flow. In ECCV, 2010.
[4] T. Chin, H. Wang and D. Suter. The ordered residual kernel for robust motion subspace clustering. In NIPS, 2009.
