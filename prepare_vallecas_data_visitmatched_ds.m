clear
datadir = 'C:\Users\justoh\Data\vallecas_study';
datafn = 'Summarytable_20200429.xlsx';

% Sheet 1: Orig converters
% Sheet 2: New converters
% Sheet 3: Converters_2020
% Sheet 4: V4+Converters at V1
% Sheet 5: Late_matched_controls 
% Sheet 6: Controls

for i = 1:6
   [num{i},txt{i},raw{i}] = xlsread(fullfile(datadir,datafn),i);
end
% check if the data appears to be matching between the different sheets
% THIS IS OK, CHECKED 9TH MARCH 2020
if 0
for i = 1:2
    matches{i} = zeros(size(raw{1},2));
    for j = 1:size(raw{1},2)
        for k = 1:size(raw{1},2)
            matches{i}(j,k) = strcmp(raw{i}{1,j},raw{i + 1}{1,k});
        end
    end
end
end
% ******************************
% age is the first column of num
% ******************************

sz = [23 26 16 28 42 720];
for i = 1:6
    for j = 1:sz(i)
         subjid{i}(j) = str2num(raw{i}{j + 1,1});
    end
end

% ****************************
% subjid's do not intersect, checked 9th march 2020 
% *****************************

for i = 1:6
    age{i} = num{i}(:,1);
    gender{i} = num{i}(:,2);
end
% !!!!!!!!!!!!!!!!!!
% Now controls are in sheet 6 !!!!!!!!
% !!!!!!!!!!!!!!!!!!
% find age and gender matched control sample
% idx = vallecas_matched_sample([age{2}; age{3}],[gender{2}; gender{3}],age{4},gender{4},0);
% idx_subjid = subjid{4}(idx);
% save(fullfile(datadir,'matched_controlidx_newconverters'),'idx','idx_subjid');
% total 1264 ? variables, see below
% idx's refer to num
% demographics 
varidx{1} = [1 2 4 5];
% demographics without APOE
varidx{2} = [1 2 5];
% neuropsychology; this is just 9 variables send by Linda (see below)
% corresponding to 1264 variables mentioned in the paper
varidx{3} = [6 7 11 12 13 14 15 16 18];

% MRI hippocampus total and GM-density entorhinal cortex
varidx{4} = [23 25:1272];
% % handle NaNs

% ***************************************************************
% variable names
% ***************************************************************
for k = 1:length(varidx)
    varnames{k} = raw{1}(1,varidx{k} + 3);
    
    for j = 2:6
        rrr = raw{j}(1,varidx{k} + 3 + (j == 4));
        for i = 1:length(varnames{k})
            ss{k,j}(i) = strcmp(varnames{k}{i},rrr{i});
        end
    end
end
        

disp('Varnames demographics');
disp(varnames{1})
disp('Varnames np');
disp(varnames{2})

% save(fullfile(datadir,'vdata_visitmatched_v4_varnames'),'varnames');

    
% LABELS:
% 1: Train converters
% 2: Train controls
% 3: Test converters
% 4: Test controls
% 5: v4+ converters at v1 (Not used)

% 1264 variables from Linda's email::
% 1248 EC GM-density values
% 1 hippocampal volume
% 4 demographics (age, sex, education, APOE)
% 9 neuropsych (MMSE, FAQ, 3 Rey complex figure scores, FCSRTrecdiftot only, semantic fluency, phonetic fluency, STAI trait)
% 1 subject ID
% 1 diagnosis (control or converter)

vdata.label_desc = '1: Train converters, 2: Train control 3: Test converters 4: Test controls 5: V4+ converters at V1 (not used)';
vdata.data_desc = '1 demographics + apoe, 2 demographics, 3:neuropsychology,4: MRI';
vdata.label = [ones(size(num{1},1),1); 3*ones(size(num{2},1),1); 3*ones(size(num{3},1),1); 5*ones(size(num{4},1),1); 4*ones(size(num{5},1),1); 2*ones(size(num{6},1),1)];
vdata.subjno =[subjid{1},subjid{2},subjid{3},subjid{4},subjid{5},subjid{6}];
% vdata.label(idx + size(num{1},1) + size(num{2},1) + size(num{3},1)) = 4;
% test controls now num{5}
for i = 1:length(varidx)
       vdata.data{i} = [num{1}(:,varidx{i}); num{2}(:,varidx{i}); num{3}(:,varidx{i}); num{4}(:,varidx{i}); num{5}(:,varidx{i}); num{6}(:,varidx{i})];
end
idx12 = find(vdata.label < 3);
for i = 1:length(varidx)
    vdata.data{i} = knn_imputation(vdata.data{i},idx12,1:length(varidx{i}),3);
end

    
save(fullfile(datadir,'vdata_visitmatched_v4_subjno'),'vdata');
idx12 = find(vdata.label < 3);
idx52 = find(vdata.label == 2 | vdata.label == 5); 
for iii = 1:100
    foldid{iii} = balanced_crossval(vdata.label(idx12),10,[],0,0);
end
vdata.foldid12 = foldid;
for iii = 1:100
    foldid{iii} = balanced_crossval(vdata.label(idx52),10,[],0,0);
end
vdata.foldid52 = foldid;
save(fullfile(datadir,'vdata_visitmatched_v4_subjno'),'vdata');

