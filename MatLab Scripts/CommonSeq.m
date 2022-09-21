function CommonSeq(varargin)
%Usage:

% CommonSeq('inname1','filename1.txt','inname2','filename2.txt','inname3','filename3.txt','indir1,'path1','indir2','path2','indir3','path3')
% CommonSeq(...,'cutoff',n)
% CommonSeq(...,'top',m)
% CommonSeq('cter','XXX')

%input(optional)
%
% CommonSeq('inname1','filename1.txt','inname2','filename2.txt','inname3','filename3.txt','indir1,'path1','indir2','path2','indir3','path3')
%       specifies three files and three paths corresponding to them. If not
%       specified, dialog boxes will open for each.
% CommonSeq(...,'cutoff',n), where n is the minimum abundance to be considered
% CommonSeq(...,'top',m), alternative to the previous one, it indicates the top m abundant sequences of each file will be considered
% CommonSeq('cter','XXX'), specifies constant C-terminal residues (allows the removal of frame-shifted clones that do not have them)

%output
%
% A new folder named "comparison" within the folder containing the FIRST FILE. The following files are generated:
%       Comparison_seq1.txt, Comparison_seq2.txt, Comparison_seq3.txt = contain sequences that appeared only in the first, second and third file respectively
%       Comparison_seq12.txt, Comparison_seq13.txt, Comparison_seq23.txt = contain sequences that appeared in two of the files
%       Comparison_seq123.txt = contains sequences that appeared in the three files
%       Comparison_stats.txt = contains the number of total/different sequences considered in each case and the number of different sequences assigned to each file


%% INPUT SECTION
inname1 = '';
inname2 = '';
inname3 = '';
outdir = ''; % default save directory is the same as input directory
outname= 'Comparison';   % default save name
top = '';
CUTOFF = 0;
Cter = '';
M = 0; %sequences for dataset 3
N = 0; %sequences for dataset 3

% check for input variable
if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); 
    end

    % read input variables
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'inname1', inname1 = varargin{ni+1};
            case 'inname2', inname2 = varargin{ni+1};
            case 'inname3', inname3 = varargin{ni+1};
            case 'outdir', outdir = varargin{ni+1};
            case 'indir1', indir1=varargin{ni+1};
            case 'indir2', indir2=varargin{ni+1};
            case 'indir3', indir3=varargin{ni+1};
            case 'cutoff', CUTOFF=varargin{ni+1};
            case 'top', top=varargin{ni+1};
            case 'cter', Cter=varargin{ni+1};
        end
    end
end

% check whether inname was defined
if strcmp(inname1,'')
    [inname1,indir1,~] = uigetfile('*.txt','Select first file');
else
    [~,message] = fopen(fullfile(indir1, inname1));
    if strcmp(message,'') == 0
        display('File not found, a dialog box will open...');
        [inname1,indir1,~] = uigetfile('*.txt','Select first file');
    end;
end;

if strcmp(inname2,'')
    [inname2,indir2,~] = uigetfile('*.txt','Select second file');
else
    [~,message] = fopen(fullfile(indir2, inname2));
    if strcmp(message,'') == 0
        display('File not found, a dialog box will open...');
        [inname2,indir2,~] = uigetfile('*.txt','Select second file');
    end;
end;

if strcmp(inname3,'')
    [inname3,indir3,filterindex] = uigetfile('*.txt','Select third file or press cancel to compare only two files');
else
    [~,message] = fopen(fullfile(indir3, inname3));
    if strcmp(message,'') == 0
        display('File not found, a dialog box will open...');
        [inname3,indir3,filterindex] = uigetfile('*.txt','Select third file or press cancel to compare only two files');
    end;
end;


if strcmp(outdir,'')
    mkdir(indir1, 'Comparison');
    outdir = 'Comparison';
end


%% DATA READING
%open file and read data

%Reading 1st file
file = fopen(fullfile(indir1, inname1));
AllVar1 = textscan(file, '%s %d %s %*[^\n]');
fclose('all');

AllSeq1 = AllVar1{1}; %Sequences are stored as a cell array of strings
AllOccur1 = AllVar1{2};
AllNtd1 = AllVar1{3};

if strcmp(Cter,'')
else
    KEEP = ~cellfun('isempty',(strfind(AllSeq1,Cter))); %I also discard sequences withtout GGSG
    AllSeq1 = AllSeq1(KEEP);
    AllOccur1 = AllOccur1(KEEP);
    AllNtd1 = AllNtd1(KEEP);
end;

if strcmp(top,'') % no top number sequences identified    
    KEEP = find(AllOccur1>=CUTOFF); %I discard sequences that appeared less than 'CUTOFF' times
    AllSeq1 = AllSeq1(KEEP);
    AllOccur1 = AllOccur1(KEEP);
    AllNtd1 = AllNtd1(KEEP);
else
    if numel(AllSeq1) < top
        CUTOFF = 0;
    else
    CUTOFF = AllOccur1(top);
    display(['Minimum abundance for dataset 1 = ' num2str(CUTOFF)]);
    end;
    KEEP = find(AllOccur1>=CUTOFF); %I discard sequences that appeared less than 'CUTOFF'
    AllSeq1 = AllSeq1(KEEP);
    AllOccur1 = AllOccur1(KEEP);
    AllNtd1 = AllNtd1(KEEP);
end;


K = numel(AllSeq1); %total different sequences
display(['File 1: considering ' num2str(K) ' different sequences']);
I = sum(AllOccur1); %total sequences considered
display(['File 1: considering ' num2str(I) ' total sequences']);

%Reading 2nd file
file = fopen(fullfile(indir2, inname2));
AllVar2 = textscan(file, '%s %d %s %*[^\n]');
fclose('all');

AllSeq2 = AllVar2{1}; %Sequences are stored as a cell array of strings
AllOccur2 = AllVar2{2};
AllNtd2 = AllVar2{3};

if strcmp(Cter,'')
else
KEEP = ~cellfun('isempty',(strfind(AllSeq2,Cter))); %I also discard sequences withtout GGSG
AllSeq2 = AllSeq2(KEEP);
AllOccur2 = AllOccur2(KEEP);
AllNtd2 = AllNtd2(KEEP);
end;

if strcmp(top,'') % no top number sequences identified    
    KEEP = find(AllOccur2>=CUTOFF); %I discard sequences that appeared less than 'CUTOFF' times
    AllSeq2 = AllSeq2(KEEP);
    AllOccur2 = AllOccur2(KEEP);
    AllNtd2 = AllNtd2(KEEP);
else
    if numel(AllSeq2) < top
        CUTOFF = 0;
    else
    CUTOFF = AllOccur2(top);
    display(['Minimum abundance for dataset 2 = ' num2str(CUTOFF)]);
    end;
    KEEP = find(AllOccur2>=CUTOFF); %I discard sequences that appeared less than 'CUTOFF' times
    AllSeq2 = AllSeq2(KEEP);
    AllOccur2 = AllOccur2(KEEP);
    AllNtd2 = AllNtd2(KEEP);
end;

L = numel(AllSeq2); %total different sequences
display(['File 2: considering ' num2str(L) ' different sequences']);
J = sum(AllOccur2); %total sequences considered
display(['File 2: considering ' num2str(J) ' total sequences']);


%Reading 3rd file
if filterindex == 0
    AllSeq3 = {''};
    AllOccur3 = 0;
    AllNtd3 = {''};
else
    file = fopen(fullfile(indir3, inname3));
    AllVar3 = textscan(file, '%s %d %s %*[^\n]');
    fclose('all');

    AllSeq3 = AllVar3{1}; %Sequences are stored as a cell array of strings
    AllOccur3 = AllVar3{2};
    AllNtd3 = AllVar3{3};
    
    if strcmp(Cter,'')
    else
        KEEP = ~cellfun('isempty',(strfind(AllSeq1,Cter))); %I also discard sequences withtout GGSG
        AllSeq1 = AllSeq1(KEEP);
        AllOccur1 = AllOccur1(KEEP);
        AllNtd1 = AllNtd1(KEEP);
    end;
    
    if strcmp(top,'') % no top number sequences identified    
        KEEP = find(AllOccur3>=CUTOFF); %I discard sequences that appeared less than 'CUTOFF' times
        AllSeq3 = AllSeq3(KEEP);
        AllOccur3 = AllOccur3(KEEP);
        AllNtd3 = AllNtd3(KEEP);
    else
        if numel(AllSeq3) < top
        CUTOFF = 0;
        else
        CUTOFF = AllOccur3(top);
        display(['Minimum abundance for dataset 3 = ' num2str(CUTOFF)]);
        end;
        
        KEEP = find(AllOccur3>=CUTOFF); %I discard sequences that appeared less than 'CUTOFF' times
        AllSeq3 = AllSeq3(KEEP);
        AllOccur3 = AllOccur3(KEEP);
        AllNtd3 = AllNtd3(KEEP);
    end;
    
    M = numel(AllSeq3); %total different sequences
    display(['File 3: considering ' num2str(M) ' different sequences']);
    N = sum(AllOccur3); %total sequences considered
    display(['File 3: considering ' num2str(N) ' total sequences']);
end;




%% DATA ANALYSIS

ii=0; iii=0; iv=0; v=0; vi=0; vii=0;
Seq1={}; Seq2={}; Seq3={}; Seq12={}; Seq13={}; Seq23={}; Seq123={};

%Found how many times each element of dataset1 is found in the other datasets
for i=1:numel(AllNtd1)
    Repe12 = strcmp(AllNtd1{i},AllNtd2);
    Repe13 = strcmp(AllNtd1{i},AllNtd3);
    if sum(Repe12) == 0; %sequence not present in 2
        if sum(Repe13) == 0; %sequence present only in 1
            ii=ii+1;
            Seq1{ii}= [AllSeq1{i} '      ' num2str(AllOccur1(i)) '      ' '0' '      ' '0' '      ' AllNtd1{i}];
        elseif sum(Repe13) == 1; %sequence present in 1 and 3
            iii=iii+1;
            Position = find(Repe13);
            Seq13{iii}= [AllSeq1{i} '      ' num2str(AllOccur1(i)) '      ' '0' '      ' num2str(AllOccur3(Position)) '      ' AllNtd1{i}];
            %eliminate this sequence from dataset3
            Repe13 = Repe13-1; %I transform 1s to 0s
            KEEP = find(Repe13);
            AllSeq3 = AllSeq3(KEEP);
            AllNtd3 = AllNtd3(KEEP);
            AllOccur3 = AllOccur3(KEEP);            
        end;
        
    elseif sum(Repe12) == 1; %sequence present in 2
        if sum(Repe13) == 0; %sequence present in 1 and 2
            iv=iv+1;
            Position = find(Repe12);
            Seq12{iv}= [AllSeq1{i} '      ' num2str(AllOccur1(i)) '      ' num2str(AllOccur2(Position)) '      ' '0' '      ' AllNtd1{i}];
        elseif sum(Repe13) == 1; %sequence present in 1,2 and 3
            v=v+1;
            Position2 = find(Repe12);
            Position3 = find(Repe13);
            Seq123{v}= [AllSeq1{i} '      ' num2str(AllOccur1(i)) '      ' num2str(AllOccur2(Position2)) '      ' num2str(AllOccur3(Position3)) '      ' AllNtd1{i}];
            %eliminate this sequence from dataset3
            Repe13 = Repe13-1; %I transform 1s to 0s
            KEEP = find(Repe13);
            AllSeq3 = AllSeq3(KEEP);
            AllNtd3 = AllNtd3(KEEP);
            AllOccur3 = AllOccur3(KEEP); 
        end;
        %eliminate this sequence from dataset2
        Repe12 = Repe12-1; %I transform 1s to 0s
        KEEP = find(Repe12);
        AllSeq2 = AllSeq2(KEEP);
        AllNtd2 = AllNtd2(KEEP);
        AllOccur2 = AllOccur2(KEEP);         
    end;
end;
display('Finished dataset 1');

%IDEM with the remaining sequences of dataset2
for i=1:numel(AllNtd2)
    Repe23 = strcmp(AllNtd2{i},AllNtd3);
    if sum(Repe23) == 0; %sequence present only in 2
        vi=vi+1;
        Seq2{vi}= [AllSeq2{i} '      ' '0' '      ' num2str(AllOccur2(i)) '      ' '0' '      ' num2str(AllNtd2{i})];
    elseif sum(Repe23) == 1; %sequence presen in 2 and 3
        vii=vii+1;
        Position = find(Repe23);
        Seq23{vii} = [AllSeq2{i} '      ' '0' '      ' num2str(AllOccur2(i)) '      ' num2str(AllOccur3(Position)) '      ' AllNtd2{i}];
        %eliminate this sequence from dataset3
            Repe23 = Repe23-1; %I transform 1s to 0s
            KEEP = find(Repe23);
            AllSeq3 = AllSeq3(KEEP);
            AllNtd3 = AllNtd3(KEEP);
            AllOccur3 = AllOccur3(KEEP); 
    end;
end;
display('Finished dataset 2');

%IDEM remaining sequences of dataset3
if filterindex == 0
else
    for i=1:numel(AllNtd3)
        Seq3{i}= [AllSeq3{i} '      ' '0' '      ' '0' '      ' num2str(AllOccur3(i)) '      ' AllNtd3{i}];
    end;
    display('Finished dataset 3');
end;



%% PRINTING THE FILES

fh = fopen(fullfile(indir1, outdir,[outname '_stats.txt']),'a');
fprintf(fh, '\r\n\r\ndifferent sequences\r\n');
fprintf(fh, 'seq.considered 1: %d\r\nseq.considered 2: %d\r\nseq.considered 3: %d\r\n', [ K, L, M]);
fprintf(fh, '\r\ntotal sequences\r\n');
fprintf(fh, 'seq.considered 1: %d\r\nseq.considered 2: %d\r\nseq.considered 3: %d\r\n', [ I, J, N]);
fprintf(fh, '\r\nseq. exclusively in 1: %d \r\nseq. exclusively in 2: %d \r\nseq. exclusively in 3: %d \r\n\r\n', [numel(Seq1), numel(Seq2), numel(Seq3)]);
fprintf(fh, 'seq. in 1 and 2: %d \r\nseq. in 1 and 3: %d \r\nseq. in 2 and 3: %d \r\n\r\n', [numel(Seq12), numel(Seq13), numel(Seq23)]);
fprintf(fh, 'seq. in 1,2,3: %d \r\n\r\n', numel(Seq123));
fclose('all');

fh = fopen(fullfile(indir1, outdir,[outname '_Seq1.txt']),'w');
for i=1:numel(Seq1)
    fprintf(fh, '%s\r\n', Seq1{i});
end
fclose('all');

fh = fopen(fullfile(indir1, outdir,[outname '_Seq2.txt']),'w');
for i=1:numel(Seq2)
    fprintf(fh, '%s\r\n', Seq2{i});
end
fclose('all');

if filterindex == 0
else
    fh = fopen(fullfile(indir1, outdir,[outname '_Seq3.txt']),'w');
    for i=1:numel(Seq3)
        fprintf(fh, '%s\r\n', Seq3{i});
    end
    fclose('all');
end;

fh = fopen(fullfile(indir1, outdir,[outname '_Seq12.txt']),'w');
for i=1:numel(Seq12)
    fprintf(fh, '%s\r\n', Seq12{i});
end
fclose('all');

if filterindex == 0
else
    fh = fopen(fullfile(indir1, outdir,[outname '_Seq13.txt']),'w');
    for i=1:numel(Seq13)
        fprintf(fh, '%s\r\n', Seq13{i});
    end
    fclose('all');

    fh = fopen(fullfile(indir1, outdir,[outname '_Seq23.txt']),'w');
    for i=1:numel(Seq23)
        fprintf(fh, '%s\r\n', Seq23{i});
    end
    fclose('all');

    fh = fopen(fullfile(indir1, outdir,[outname '_Seq123.txt']),'w');
    for i=1:numel(Seq123)
        fprintf(fh, '%s\r\n', Seq123{i});
    end
    fclose('all');
end;

% fh = fopen(fullfile(indir1, outdir,[outname '_merged.txt']),'w');
% for i=1:numel(Seq123)
%     fprintf(fh, '%s\r\n', Seq123{i});
% end
% for i=1:numel(Seq12)
%     fprintf(fh, '%s\r\n', Seq12{i});
% end
% for i=1:numel(Seq23)
%     fprintf(fh, '%s\r\n', Seq23{i});
% end
% for i=1:numel(Seq13)
%     fprintf(fh, '%s\r\n', Seq13{i});
% end
% for i=1:numel(Seq1)
%     fprintf(fh, '%s\r\n', Seq1{i});
% end
% for i=1:numel(Seq2)
%     fprintf(fh, '%s\r\n', Seq2{i});
% end
% for i=1:numel(Seq3)
%     fprintf(fh, '%s\r\n', Seq3{i});
% end
% fclose('all');


end

