%%
% In this file we compute firstly all the selection criteria and then we find the 
% Silhouette score of the sequences wrt the model selected by each criterium
% In this version we don't normalize the process of model fitting in the
% second phase
%%
%% Obtain silhouette graphs and tables for all clips, for all models
clear all;
close all;

max_NumHypoPerFrame = 500;
motions = ["motion1", "motion2","motion3","motion4", "motion5"];
n1 = 1;
n2 = 7;

%% Load Seq Information
temp = load('../../../Data/SeqList.mat');
SeqList = temp.SeqList;
FrameGap = 1;
addpath '..\..\..\Tools\multigs\model_specific';
%% models silhouette
%Input theg desired geometric model to analyze
models = ["fundamental", "fundamentala", "fundamentalt", "homography","affine"];
scoresMotion = []; 
silhouetteMotion = [];
silhouetteMotionAvg = [];
silhouetteAvg = [];

ps = [7, 4, 2, 4, 3]; % also called p: minimum number of correspondences needed to estimate the constraint
ds = [3, 3, 3, 2, 2]; % dimension of the constraint:3 for fundamentals 2 for homography
ks = [7, 4, 2, 8, 6]; % number of parameters of the model:7 4 2 for fundamentals 8 for homography

% model selection algorithms
names = ["GRIC","GRIC Tuned", "AIC", "BIC", "MDL", "GBIC","Residuals","RobustResiduals","AllFundamental","AllHomography"];

gric = [];
aic = [];
bic = [];
mdl = [];
gbic = [];
lambdas = [];

for seq = 1:length(SeqList)
    seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
    Data = seqData.Data;
    gt = seqData.Data.GtLabel;
    nMotions = max(gt);

    for f_i = 1:Data.nFrames-FrameGap
        for model=1:size(models,2)
            model_type = models(model);
            [fitModel, resModel, ~, ~, ~] = getModelParam(model_type);
            %% Select points for every motion visible on both frames
            visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);
            Models = [];
            matchedPoints1 = [];
            matchedPoints2 = [];
            % obtain a binary array containing indices of point divided per
            % motions and for every motion estimate the F
            sizes = [];
            for m = 1:nMotions
                pointsInMotion = visible_pts_ind & (gt==m);
                % matched points of motion m visible in frames f_i and f_i+1
                matchedPoints1{m} = Data.ySparse(:,pointsInMotion,f_i);
                matchedPoints2{m} = Data.ySparse(:,pointsInMotion,f_i+FrameGap);
                if size(matchedPoints1{m}, 2) > 1
                    matchedPoints1{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i));
                    matchedPoints2{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i+FrameGap));
                end
                sizes = [sizes size(matchedPoints1{m},2)];
                % in case points of a certain motion are not visible in the
                % pair of frames or correspondences are less than 8, 
                % model is [0 0 0; 0 0 0; 0 0 1]
                if size(matchedPoints1{m}, 2) < 8
                    fprintf("Warning: %s in frames pair (%d-%d) for motion %d has %d correspondences\n", ...
                        SeqList{seq}, f_i, f_i+FrameGap, m, size(matchedPoints1{m},2));
                end
                try
                    Models{m} = feval(fitModel, [matchedPoints1{m}; matchedPoints2{m}]);
                catch 
                    Models{m} = NaN;
                end
            end  
            

            ReswrtMotion = [];
            variances = [];
            ReswrtMotionNormalized = [];
            robustRes = [];
            r = 4;
            d = ds(model); 
            k = ks(model);
            p = ps(model);

            for m=1:nMotions
                if isnan(Models{m})
                    grictunedtmp{f_i,model,m} = Inf;
                    grictmp{f_i,model,m} = Inf;
                    aictmp{f_i,model,m} = Inf;
                    bictmp{f_i,model,m} = Inf;
                    mdltmp{f_i,model,m} = Inf;
                    gbictmp{f_i,model,m} = Inf;
                    residualstmp{f_i,model,m} = Inf;
                    robustresidualstmp{f_i,model,m} = Inf;
                else
                    ReswrtMotion{m} = feval(resModel, Models{m},[matchedPoints1{m}; matchedPoints2{m}]);
                    variances{m} = var(ReswrtMotion{m});
                    ReswrtMotionNormalized{m} = ReswrtMotion{m}.^2./variances{m};
                    n = size(ReswrtMotionNormalized{m},1);
                    % hyperparameter tuning for GRIC
                    l = load("BestLambdas.mat").lambdas;
                    lambda1T = l(1);
                    lambda2T = l(2);
                    lambda3T = l(3);
                    
                    lambda1 = 2;
                    lambda2 = 2;
                    lambda3 =  2;
                    lambdas = [lambda1, lambda2, lambda3];
                    tmp = lambda2*(r-d);
                    tmpT = lambda2T*(r-d);
                    robustResT{m} = min(ReswrtMotionNormalized{m}, tmpT);
                    robustRes{m} = min(ReswrtMotionNormalized{m}, tmp);
                    grictmp{f_i,model,m} = sum(robustRes{m}) + lambda1*d*n + lambda3*k;
                    grictunedtmp{f_i,model,m} = sum(robustResT{m}) + lambda1T*d*n + lambda3T*k;
                    aictmp{f_i,model,m} = sum(ReswrtMotionNormalized{m}) + 2*p;
                    bictmp{f_i,model,m} = sum(ReswrtMotionNormalized{m}) + p*log(r*n);
                    mdltmp{f_i,model,m} = sum(ReswrtMotionNormalized{m}) + p*0.5*log2(r*n);
                    gbictmp{f_i,model,m} = sum(ReswrtMotionNormalized{m}) + log(r)*d*n + log(r*n)*k;
                    residualstmp{f_i,model,m} = sum(ReswrtMotionNormalized{m});
                    robustresidualstmp{f_i,model,m} = sum(robustRes{m});
                end
            end
        end
    end
    grictuned{seq} = grictunedtmp;
    gric{seq} = grictmp;
    aic{seq} = aictmp;
    bic{seq} = bictmp;
    mdl{seq} = mdltmp;
    gbic{seq} = gbictmp;
    residuals{seq} = residualstmp;
    robustresiduals{seq} = robustresidualstmp;
end 
% --------------------------------------------------------------------------------
mod_sel = containers.Map([1,2,3,4,5,6,7,8], {gric, grictuned aic, bic, mdl, gbic, residuals,robustresiduals});
% find the best model clip - frame - motion for each criterium 
for i=1:length(mod_sel)
    currModel = mod_sel(i);
    fprintf("Testing "+names(i)+"\n");
    for seq=1:length(SeqList)
        seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Data = seqData.Data;
        nMotions = max(Data.GtLabel);
        results = currModel{seq};
        best_models = [];
        fr_list = [];
        for fr=1:Data.nFrames-FrameGap
            fr_list = [fr_list; "frame "+int2str(fr)+"-"+int2str(fr+1)];
            m = cell2mat(squeeze(results(fr,:,:)));
            [val,idx] = min(m);
            curr_best=[];
            for motion=1:nMotions
                best_model = models(idx(motion));
                curr_best=[curr_best best_model];
            end
            best_models=[best_models;curr_best];
        end
        T = array2table(best_models);
        T.Properties.VariableNames = motions(1:nMotions);
        T.Properties.RowNames = [fr_list];
        T.Properties.DimensionNames(1) = SeqList(seq);
        if ~exist("ModelSelectionResults", 'dir')
            mkdir("ModelSelectionResults");
        end
        %writetable(T,"ModelSelectionResults/"+names(i)+".xls",'WriteRowNames',true,'Sheet',seq)
        seq_best_models{i,seq} = best_models;
    end
end
%% add all fundamental and all homography to the criteria for benchmark
constantModels = ["Fundamental", "Homography"];
for m = 1:size(constantModels,2)
    for i = 1:22 
        nrows = size(seq_best_models{1,i},1);
        ncols = size(seq_best_models{1,i},2);
        M = strings([nrows,ncols]);
        M(1:nrows,1:ncols) = constantModels(m);
        seq_best_models{8+m,i} = M;
    end
end
totalmodels = length(mod_sel) + length(constantModels);

% save the best models obtained
save("seq_best_models","seq_best_models");

%% Silhouette
for i=1:totalmodels
    for seq = 1:length(SeqList)
        seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Data = seqData.Data;
        gt = seqData.Data.GtLabel;
        nMotions = max(gt);
        silhouettes = [];
        curr_best_models = seq_best_models{i, seq};
        for f_i = 1:Data.nFrames-FrameGap
            %% Select points for every motion visible on both frames
            visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);
            Models = [];
            matchedPoints1 = [];
            matchedPoints2 = [];
            sizes = [];
            % obtain a binary array containing indices of point divided per
            % motions and for every motion estimate the F
            for m = 1:nMotions
                pointsInMotion = visible_pts_ind & (gt==m);
                % matched points of motion m visible in frames f_i and f_i+1
                matchedPoints1{m} = Data.ySparse(:,pointsInMotion,f_i);
                matchedPoints2{m} = Data.ySparse(:,pointsInMotion,f_i+FrameGap);
                if size(matchedPoints1{m}, 2) > 1
                    matchedPoints1{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i));
                    matchedPoints2{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i+FrameGap));
                end
                sizes = [sizes size(matchedPoints1{m},2)];

              
                model_type = curr_best_models(f_i,m);
                [fitModel, resModel, ~, ~, ~] = getModelParam(model_type);
                try

                    Models{m} = feval(fitModel, [matchedPoints1{m}; matchedPoints2{m}]);
                catch 
                    Models{m} = NaN;
                end
            end
           
            % now compute residuals of every point wrt the model for every motion.
            % This is needed to compute a and b for every point
            % ReswrtMotion: matrix m x n: every row corresponds to the residuals of
            % a motion wrt all the m motions. e.g. 2 motions => 2x2:
            % [res(points_motion_1,model_1) res(points_motion_1,model_2);
            % res(points_motion_2,model_1) res(points_motion_2,model_2) ]
            ReswrtMotion = [];
            for m=1:nMotions
                for n=1:nMotions
                    model_type = curr_best_models(f_i,n);
                    [~, resModel, ~, ~, ~] = getModelParam(model_type);
                    if isnan(Models{n})
                        ReswrtMotion{m}{n} = Inf(size(matchedPoints1{m},2),1);
                    else
                        ReswrtMotion{m}{n} = feval(resModel, Models{n},[matchedPoints1{m}; matchedPoints2{m}]);
                    end
                end
            end


            %now compute a and b for every point and its silhouette score     
            for m=1:nMotions
                b = [];
                for n=1:nMotions
                    if m==n % vector of a's for motion m
                        a = ReswrtMotion{m}{n};                    
                    else %residuals of points of motion m wrt motion n
                        if isempty(b)
                            b = ReswrtMotion{m}{n};
                        else
                            b = min(b, ReswrtMotion{m}{n});
                        end
                    end
                end
                % we have a and b for every point of motion m. Now silhouette
                scoresMotion{seq}{f_i}{m} = (b - a) ./ max(a, b);
                silhouetteMotion{seq}{f_i}{m} = mean(scoresMotion{seq}{f_i}{m});
                silhouettes(m, f_i) = mean(scoresMotion{seq}{f_i}{m});
                frames(f_i) = "frames "+f_i+" "+(f_i+1);
            end
        end
        silhouetteMotionAvg{seq,i} = nanmean(silhouettes,2); 
        silhouetteMethodAvg{seq,i} = mean(silhouetteMotionAvg{seq, i});
    end
end

meanSilhouettes = mean(cell2mat(silhouetteMethodAvg));


TTotal = array2table([cell2mat(silhouetteMethodAvg); mean(cell2mat(silhouetteMethodAvg))]);
TTotal.Properties.VariableNames = names;
TTotal.Properties.RowNames = [SeqList "avg"];
writetable(TTotal,"silhouetteScoresModelSelection.xls",'WriteRowNames',true,'Sheet',1);
TTotal
%% Residui
for i=1:totalmodels
    for seq = 1:length(SeqList)
        seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Data = seqData.Data;
        gt = seqData.Data.GtLabel;
        nMotions = max(gt);
        silhouettes = [];
        curr_best_models = seq_best_models{i, seq};
        for f_i = 1:Data.nFrames-FrameGap
            %% Select points for every motion visible on both frames
            visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);
            Models = [];
            matchedPoints1 = [];
            matchedPoints2 = [];
            ReswrtMotion = [];
            % obtain a binary array containing indices of point divided per
            % motions and for every motion estimate the F
            for m = 1:nMotions
                pointsInMotion = visible_pts_ind & (gt==m);
                % matched points of motion m visible in frames f_i and f_i+1
                matchedPoints1{m} = Data.ySparse(:,pointsInMotion,f_i);
                matchedPoints2{m} = Data.ySparse(:,pointsInMotion,f_i+FrameGap);
                if size(matchedPoints1{m}, 2) > 1
                    matchedPoints1{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i));
                    matchedPoints2{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i+FrameGap));
                end
                
                model_type = curr_best_models(f_i,m);
                [fitModel, resModel, ~, ~, ~] = getModelParam(model_type);
                try
                    Models{m} = feval(fitModel, [matchedPoints1{m}; matchedPoints2{m}]);
                    ReswrtMotion{m} = feval(resModel, Models{m},[matchedPoints1{m}; matchedPoints2{m}]);                    
                catch 
                    Models{m} = NaN;
                    ReswrtMotion{m} = NaN(size(matchedPoints1{m},2),1);
                end
                avgRes{m} = mean(ReswrtMotion{m});
            end     
 
            avgResFrames = nanmean(cell2mat(avgRes));
            
            % we have a and b for every point of motion m. Now silhouette
            scoresMotion{seq}{f_i} = avgResFrames;
           
        end
        scoresMotionAvg{seq,i} = mean(cell2mat(scoresMotion{seq}));
    end   
end


meanFinal = cell2mat(scoresMotionAvg);


TTotalRes = array2table([meanFinal; nanmean(meanFinal)]);
TTotalRes.Properties.VariableNames = names;
TTotalRes.Properties.RowNames = [SeqList "avg"];
writetable(TTotalRes,"residualsScoresModelSelection.xls",'WriteRowNames',true,'Sheet',1);
TTotalRes







