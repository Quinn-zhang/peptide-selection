function Clustering_v2(varargin)

%Usage:
%   clustering('inname','filename.txt','indir','path')
%   clustering(...,'number_dif',n)
%   clustering(...,'min_clustersize',m)
%   clustering(...,'cter','NNNN')
%   clutsering(...,'stringency',s)
%   clustering(...,'min_abun',n)
%   clustering(...,'logos','off')
%   clustering(...,'gappen',n)
%
%Input (optional):
%   clustering('inname','filename.txt','indir','path') indicates the file 
%       name and path to the folder where it is located. If not specified, 
%       a dialog box to choose the file will open.
%       clustering can read the output files from Step2, LoopLengths and
%       FindSeq. Requisites: data is on the format peptide seq - abundance
%       - nucleotide seq.
%   clustering(...,'number_dif',n) indicates how many different sequences
%       will be clustered. If not specified, n = 200, i.e. top 200 most
%       abundant sequences will be clustered.
%   clustering(...,'min_clustersize',m) indicates the minimum number of
%       sequences within a cluster to be considered. If a cluster has less than
%       m sequences, it will be transferred to the "mixed" cluster.
%   clustering(...'cter','XXXX') indicates a constant C-terminal region of the
%       peptide.Peptides without this constant region will be not considered.
%   clustering(...,'min_abun',n) indicates the minimum abundance for clones
%       to be considered.
%   clutsering(...,'stringency',s) allows to fine-tune the clustering of
%       the script to different datasets. In general, higher values of
%       stringency will lead to more similar peptides within each cluster and
%       more sequences in the mixed cluster. Lower values of stringency allow
%       more differences within each cluster and as a result less sequences go
%       to the mixed cluster.
%   clustering(...,'logos','off') disables the generation of logo picture
%       files
%   clustering(...,'gappen',n) changes the value of gap opening and gap
%       extension penalties. Default value is 8.
%
%Output:
%   o	Clusters_filename.txt file within the same folder as the input file.
%   o	A series of .jpg files corresponding to the sequence logos of each group within the same folder as the input file.
   




%% INPUT SECTION
inname = '';  
outdir = ''; % default save directory is the same as input directory
outname= 'Clusters';   % default save name
indir = '';
MAXCLUST = '';
cter = 0;
min_clustersize = 3;
stringency = 0.5;
CutOff = '';
min_abun = 0;
logos = 'on';
GapOpen = 8;
RemNter = '';
RemCter = '';

% check for input variable
if exist('varargin','var')
    L = length(varargin);
    %if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); 
    %end

    % read input variables
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'inname', inname = varargin{ni+1};
            case 'outdir', outdir = varargin{ni+1};
            case 'indir', indir=varargin{ni+1};
            case 'number_dif', CutOff=varargin{ni+1}; %number of different sequences you want to end up having
            case 'maxclust', MAXCLUST=varargin{ni+1};
            case 'min_clustersize', min_clustersize=varargin{ni+1};
            case 'cter', cter=varargin{ni+1};
            case 'stringency', stringency=varargin{ni+1};
            case 'min_abun', min_abun=varargin{ni+1};
            case 'logos', logos=varargin{ni+1};
            case 'gappen', GapOpen=varargin{ni+1};
            case 'remconstant', RemNter=varargin{ni+1}; RemCter=varargin{ni+2};
        end
    end
end

% check whether inname was defined
if strcmp(inname,'')
    [inname,indir,~] = uigetfile('*.txt','Select file');
else
    [~,message] = fopen(fullfile(indir, inname));
    if strcmp(message,'') == 0
        display('File not found, a dialog box will open...');
        [inname,indir,~] = uigetfile('*.txt','Select file');
    end;
end;

if strcmp(outdir,'') == 1
    outdir = indir;
end

% check how many different sequences will be compared
if min_abun == 0
    if strcmp(CutOff,'') == 1
        CutOff = 200;
        display('Comparing top 200 abundant sequences. To change, specify number_dif');
    else
        display(['Comparing top ' num2str(CutOff) ' abundant sequences.']);
    end;
else
    display(['Comparing sequences whose abundance is higher or equal to ' num2str(min_abun)]);
end;

% check clustering parameters:
display(['Minimum cluster size = ' num2str(min_clustersize) ', (smaller clusters will be merged in "cluster mixed")']);
display(['Stringency = ' num2str(stringency)]);


% check whether C-ter has been specified:
if cter == 0
    display('No c-ter specified.')
else
    display(['Considering only sequences containing: ' cter '.']);
end;

% check additional parameters
if GapOpen == 8
    display('Using standard gap penalty = 8');
else
    display(['Using gap penalty = ' num2str(GapOpen)]);
end;


%% DATA READING
%open file and read data

file = fopen(fullfile(indir, inname));
AllVar = textscan(file, '%s %d %s %*[^\n]');
fclose('all');

AllSeq = AllVar{1}; %Sequences are stored as a cell array of strings
AllOccur = AllVar{2};
AllNtd = AllVar{3};

% Keep only sequences with correct C-ter
if cter == 0
else
    KEEP = ~cellfun('isempty',(strfind(AllSeq,cter)));
    AllSeq = AllSeq(KEEP);
    AllOccur = AllOccur(KEEP);
    AllNtd = AllNtd(KEEP);
end;
if min_abun == 0
        if numel(AllOccur) > CutOff
            CUTOFF = AllOccur(CutOff);
%             if CUTOFF < 3
%                 CUTOFF = 3;
%             else
%             end;
        else CUTOFF = 3;
        end;
        display(['Minimum abundance = ' num2str(CUTOFF)]);
else
    CUTOFF = min_abun;
end;

KEEP = find(AllOccur>=CUTOFF); %I discard sequences that appeared less than 'CUTOFF' times (for naive libraries set CUTOFF to 1!!)
AllSeq = AllSeq(KEEP);
AllOccur = AllOccur(KEEP);
AllNtd = AllNtd(KEEP);

total_seq_considered = sum(AllOccur);
display(['Number of total sequences considered = ' num2str(total_seq_considered)]);
total_dif_sequences = numel(AllSeq);
display(['Number of different sequences = ' num2str(total_dif_sequences)]);

%Removing constant N-ter and C-ter
if strcmp(RemNter,'')
else
    AllSeq=strrep(AllSeq,RemNter,'');
end;
if strcmp(RemCter,'')
    AllSeq=strrep(AllSeq,RemCter,'');
end;


%Setting the maximum number of clusters accordin to the number of total different sequences
if strcmp(MAXCLUST,'') == 1
    if total_dif_sequences < 500
        MAXCLUST = round(total_dif_sequences * stringency);
    else
        MAXCLUST = 250;
    end;
else
end;


%creating the names for the phylotree
AllNames = cell(numel(AllSeq),1);
for i=1:numel(AllSeq)
    AllNames{i} = [AllSeq{i} '           ' num2str(AllOccur(i)) '       ' AllNtd{i}];
end;

clear('AllVar');


%% DATA ANALYSIS

t_total=tic;

%Calculate distances
tic;
%Dist = seqpdist(ModSeq,'Method','p-distance'); %vector with all the distances 1-2 1-3 1-4 ... 1-N 2-1 2-3 ... 
Dist2 = seqpdist(AllSeq,'Method','alignment-score','GapOpen',GapOpen);
display(['Distance calculated in ' num2str(toc) ' sec']);

%Make tree
tic;
tree = seqlinkage(Dist2,'average',AllNames);
% tree = seqneighjoin(Dist2,'equivar',AllNames);
%Possible criteria = single, complete, average, weighted, centroid, median
%view(tree)
display(['Tree done in ' num2str(toc) ' sec']);


% Validate the clusters in the tree and find the best partition
tic;
[i,~,~] = cluster(tree, [],'criterion','average','maxclust',MAXCLUST);
%Possible criteria = maximum, median, average, ratio, gain, silhouette
[LeafNames] = get(tree,'LeafNames');
[i,IX] = sort(i);
LeafNames = LeafNames(IX);
display(['Clustering done in ' num2str(toc) ' sec']);



%% PRINTING THE FILES

%Printing stats file
%fh = fopen(fullfile(outdir,[outname '_stats.txt']),'a');
A=zeros(1,20);
B=zeros(1,21);
B(1)=0;
for iv=1:MAXCLUST+1
    A(iv)=numel(find(i==iv));%number of elements in each cluster
    %fprintf(fh, '%d\t', A(iv));
    B(iv+1)=B(iv)+A(iv);
end;
%C=numel(AllSeq);
%fprintf(fh, '%d \r\n', C);
%fclose('all');

%Printing the consensus file
fh = fopen(fullfile(outdir,[outname '_' inname ]),'w');

v=0; % counter to 0
for iv=2:MAXCLUST+2
    if B(iv)-B(iv-1) >= min_clustersize
        v=v+1;
        fprintf(fh, '\r\n Group%d\r\n', (v));
        if B(iv)==0
            continue
        else
            for iii=B(iv-1)+1:B(iv)
                fprintf(fh, '%s\r\n', [LeafNames{iii}]);
            end;
        end;
    else
    end;
end;
fclose('all');

fh = fopen(fullfile(outdir,[outname '_MIXED.txt']),'w');
for iv=2:MAXCLUST+2
    if B(iv)-B(iv-1) < min_clustersize
        if B(iv)==0
            continue
        else
            for iii=B(iv-1)+1:B(iv)
                fprintf(fh, '%s\r\n', [LeafNames{iii}]);
            end;
        end;
    else
    end;
end;

fclose('all');


%% 2nd clustering with MIXED

%Opening the consensus file
flogo = fopen(fullfile(outdir,[outname '_MIXED.txt']));

AllVar = textscan(flogo, '%s %d %s %*[^\n]');
fclose('all');

AllSeq = AllVar{1}; %Sequences are stored as a cell array of strings
AllOccur = AllVar{2};
AllNtd = AllVar{3};

total_dif_sequences = numel(AllSeq);
display(['Sequences remaining for the second clustering = ' num2str(total_dif_sequences)]);


%creating the names for the phylotree
AllNames = cell(numel(AllSeq),1);
for i=1:numel(AllSeq)
    AllNames{i} = [AllSeq{i} '           ' num2str(AllOccur(i)) '       ' AllNtd{i}];
end;

clear('AllVar');

%Setting the maximum number of clusters accordin to the number of total different sequences
MAXCLUST = round(total_dif_sequences * stringency);

if total_dif_sequences == 1
    fh = fopen(fullfile(outdir,[outname '_' inname]),'a');
    fprintf(fh, '%s\r\n', [AllNames{1}]);
    fclose('all');
    return
elseif total_dif_sequences == 0
    return
end;


%% DATA ANALYSIS


%Calculate distances
tic;
%Dist = seqpdist(ModSeq,'Method','p-distance'); %vector with all the distances 1-2 1-3 1-4 ... 1-N 2-1 2-3 ... 
Dist2 = seqpdist(AllSeq,'Method','alignment-score','GapOpen',GapOpen);
display(['Distance calculated in ' num2str(toc) ' sec']);

%Make tree
tic;
tree = seqlinkage(Dist2,'average',AllNames);
%Possible criteria = single, complete, average, weighted, centroid, median
display(['Tree done in ' num2str(toc) ' sec']);


% Validate the clusters in the tree and find the best partition
tic;
[i,~,~] = cluster(tree, [],'criterion','average','maxclust',MAXCLUST);
%Possible criteria = maximum, median, average, ratio, gain, silhouette
[LeafNames] = get(tree,'LeafNames');
[i,IX] = sort(i);
LeafNames = LeafNames(IX);
display(['Clustering done in ' num2str(toc) ' sec']);


display(['total time = ' num2str(toc(t_total)) ' sec']);

%% PRINTING THE FILES

%Printing stats file
%fh = fopen(fullfile(outdir,[outname '_stats.txt']),'a');
A=zeros(1,20);
B=zeros(1,21);
B(1)=0;
for iv=1:MAXCLUST+1
    A(iv)=numel(find(i==iv));%number of elements in each cluster
    %fprintf(fh, '%d\t', A(iv));
    B(iv+1)=B(iv)+A(iv);
end;
%C=numel(AllSeq);
%fprintf(fh, '%d \r\n', C);
%fclose('all');

%Printing the consensus file
fh = fopen(fullfile(outdir,[outname '_' inname]),'a');

for iv=2:MAXCLUST+2
    if B(iv)-B(iv-1) >= min_clustersize
        v=v+1;
        fprintf(fh, '\r\n Group%d\r\n', (v));
        if B(iv)==0
            continue
        else
            for iii=B(iv-1)+1:B(iv)
                fprintf(fh, '%s\r\n', [LeafNames{iii}]);
            end;
        end;
    else
    end;
end;

fprintf(fh, '\r\n GroupMIXED\r\n');
for iv=2:MAXCLUST+2
    if B(iv)-B(iv-1) < min_clustersize
        if B(iv)==0
            continue
        else
            for iii=B(iv-1)+1:B(iv)
                fprintf(fh, '%s\r\n', [LeafNames{iii}]);
            end;
        end;
    else
    end;
end;

fclose('all');

%Eliminating temporary MIXED file
delete(fullfile(outdir,[outname '_MIXED.txt']));


%% Create logos

if strcmp(logos,'off') == 1

else
    % Open the file
    fh = fopen(fullfile(outdir,[outname '_' inname]),'r');

    AllVar = textscan(fh, '%s %d %s %*[^\n]');
    fclose('all');

    AllSeq = AllVar{1}; %Sequences are stored as a cell array of strings
    AllSeq{numel(AllSeq)+1}='Gr'; %Append a last line to facilitate the logo loop

    %The first line of each group is 'GroupN'
    Groups = strfind(AllSeq,'Gr');
    Index = find(~cellfun('isempty',Groups));

    for i=1:numel(Index)-1
        LogoSeq = AllSeq(Index(i)+1:Index(i+1)-1);
        [~,LogoHandle]=seqlogo(LogoSeq,'Alphabet','AA','DisplayLogo',true,'SScorrection',false);
        Name = ['Group ' num2str(i)];
        set(LogoHandle,'Name',Name);
        print(LogoHandle,fullfile(outdir,['Logo_Group ' num2str(i)]),'-djpeg');
        close all hidden;
    end;
end;










