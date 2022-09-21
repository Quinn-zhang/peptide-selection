function LoopLengths(varargin)

%New version that works with matlab releases before and after R2012b (there was a change in function unique)

%Usage:
% LoopLengths
% LoopLengths('inname','filename.txt','indir','path')
% LoopLengths(...,'cutoff',n)
% LoopLengths(...,'cter','XXX')

%input (optional)
% LoopLengths('inname','filename.txt','indir','path') indicates the file 
%       name and path to the folder where it is located. If not specified, a dialog box to choose the file will open.
%       LoopLengths needs a file with data on the format: peptide seq - abundance - nucleotide seq
% LoopLengths(...,'cutoff',n), where n specifies the minimum abundance to 
%       be considered. If not specified, no cutoff is applied and all 
%       sequences are considered
% LoopLengths(...,'cter','XXX'), where XXX is the amino acid sequence found at the C-terminus of the peptide. Allows to remove frame-shifted clones.

%LoopLengths separates all sequences in different files according to the
%   peptide format (i.e. number of cysteines and loop lengths if 2, 3 or 4
%   cysteines.

%output
%   A new folder called "LoopLengths_filename" with a series of files named
%       LoopLengths_n_NUMBERcys.txt, where "n" indicates the format and
%       "NUMBER" the number of cysteines
%        Notation for "n":   
%                       3_twocys = C XXX C
%                       5_twocys = C XXXXX C
%
%                       3_threecys = 0x3 = CC XXX C
%                       300_threecys = 3x0 = C XXX CC                      
%                       304_threecys = 3x4 = C XXX C XXXX C
%                       305_threecys = 3x5 = C XXX C XXXXX C
%
%                       4_fourcys = 0x0x4 = CCC XXXX C
%                       400_fourcys = 0x4x0 = CC XXXX CC
%                       403_fourcys = 0x4x3 = CC XXXX C XXX C
%                       40000_fourcys = 4x0x0 = C XXXX CCC
%                       40302_fourcys = 4x3x2 = C XXXX C XXX C XX C
%
%       LoopLengths_stats.txt with the information about the total and
%       different sequences assigned to each category
%                       
%   


%% INPUT SECTION
inname = '';  
outdir = ''; % default save directory is a new directory inside the input directory
outname= 'LoopLengths';   % default save name
indir = '';
Cter = '';
CUTOFF = 0;

% check for input variable
if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); 
    end

    % read input variables
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'inname', inname = varargin{ni+1};
            case 'outdir', outdir = varargin{ni+1};
            case 'indir', indir=varargin{ni+1};
            case 'cutoff', CUTOFF=varargin{ni+1};
            case 'cter', Cter=varargin{ni+1};
        end
    end
end

% check if inname was defined
if strcmp(inname,'')
    [inname,indir,~] = uigetfile('*.txt','Select .txt file');
else
    [~,message] = fopen(fullfile(indir, inname));
    if strcmp(message,'') == 0
        display('File not found, a dialog box will open...');
        [inname,indir,~] = uigetfile('*.txt','Select .txt file');
    end;
end;


if strcmp(outdir,'')
    mkdir(indir, ['LoopLengths_'  strrep(inname,'.txt','')] );
    outdir = ['LoopLengths_'  strrep(inname,'.txt','')];
end;
mkdir(fullfile(indir, outdir, '4cys'));
mkdir(fullfile(indir, outdir, '3cys'));
mkdir(fullfile(indir, outdir, '2cys'));
mkdir(fullfile(indir, outdir, 'other'));


% check if Cter was specified
if strcmp(Cter,'')
    display('No C-terminus specified');
else
    display(['Considering sequences whose C-terminus is ' Cter]);
end;

% check if cutoff was specified
if CUTOFF == 0
    display('No abundance cutoff applied');
else
    display(['Considering sequences with a minimum abundance of ' num2str(CUTOFF)]);
end;

%% DATA READING
%open file and read data

file = fopen(fullfile(indir, inname));
AllVar = textscan(file, '%s %d %s %*[^\n]');
fclose('all');

AllSeq = AllVar{1}; %Sequences are stored as a cell array of strings
AllOccur = AllVar{2};
AllNtd = AllVar{3};

KEEP = find(AllOccur>=CUTOFF); % Discard sequences that appeared less than 'CUTOFF' times)
AllSeq = AllSeq(KEEP);
AllOccur = AllOccur(KEEP);
AllNtd = AllNtd(KEEP);

if strcmp(Cter,'')
else
    KEEP = ~cellfun('isempty',(strfind(AllSeq,Cter)));
    AllSeq = AllSeq(KEEP);
    AllOccur = AllOccur(KEEP);
    AllNtd = AllNtd(KEEP);
end;

if numel(AllSeq)>1 && ischar(AllSeq{1}) && isnumeric(AllOccur(1)) && ischar(AllNtd{1})           
    total_seq_considered = sum(AllOccur);
    display(['Total sequences considered = ' num2str(total_seq_considered)]);
else
    display('Error: file chosen not suitable for this script');
    display('Choose a file of the format: peptide seq - abundance - nucleotide seq');
    return
end;
        

clear('AllVar');


%% DATA ANALYSIS: group by number of cys
numelall = numel(AllSeq);
cys = strfind(AllSeq,'C');
%clasif = zeros(numelall,1);
%threecys=cell(numelall,2);

%count the cysteines
fourcys = cell(0);
threecys = cell(0);
twocys = cell(0);
ii=0;
iii=0;
iv=0;
for i=1:(size(AllSeq))
    if numel(cys{i})==1
        clasif(i)=1;
    elseif numel(cys{i})==4
        clasif(i)=4;
        a=cys{i}(2)-cys{i}(1)-1;
        b=cys{i}(3)-cys{i}(2)-1;
        c=cys{i}(4)-cys{i}(3)-1;
        iv=iv+1;
        fourcys{iv,1}=a*10000+b*100+c; %it stores in an array: 1st column = e.g. 30304->3x3x4, 2nd column = sequence, 3rd column = abundance, 4th column = ntds
        fourcys{iv,3}=AllSeq{i};
        fourcys{iv,2}=AllOccur(i);
        fourcys{iv,4}=AllNtd{i};  
    elseif numel(cys{i})==2
        clasif(i)=2;
        o=cys{i}(2)-cys{i}(1)-1;
        iii=iii+1;
        twocys{iii,1}=o;%it stores in an array: 1st column = e.g. 3->C3C, 2nd column = abundance, 3rd column = sequence
        twocys{iii,3}=AllSeq{i};
        twocys{iii,2}=AllOccur(i);
        twocys{iii,4}=AllNtd{i};        
    elseif isempty(cys{i})
        clasif(i)=10;
    elseif numel(cys{i})==3
        clasif(i)=0; % store them separately
        n=cys{i}(2)-cys{i}(1)-1;
        m=cys{i}(3)-cys{i}(2)-1;
        ii=ii+1; %counter
        threecys{ii,1}=n*100+m;
        threecys{ii,3}=AllSeq{i}; %it stores in an array: 1st column = e.g. 304->3x4, 2nd column = sequence, 3rd column = abundance
        threecys{ii,2}=AllOccur(i);
        threecys{ii,4}=AllNtd{i};
    else
        clasif(i)=5;
    end;
end;

manycys=AllSeq(clasif==5); A=numel(manycys);
fourcys2=AllSeq(clasif==4); B=numel(fourcys2);
twocys2=AllSeq(clasif==2); C=numel(twocys2);
onecys=AllSeq(clasif==1); D=numel(onecys);
nocys=AllSeq(clasif==10); E=numel(nocys);
threecys2=AllSeq(clasif==0); F=numel(threecys2);

many_cys = find(clasif==5);
many_cys_seq=AllSeq(many_cys);
many_cys_abundances=AllOccur(many_cys);
many_cys_ntd=AllNtd(many_cys);

one_cys = find(clasif==1);
one_cys_seq=AllSeq(one_cys);
one_cys_abundances=AllOccur(one_cys);
one_cys_ntd=AllNtd(one_cys);

three_cys = find(clasif==0);
% three_cys_seq=AllSeq(three_cys);
three_cys_abundances=AllOccur(three_cys);
% three_cys_ntd=AllNtd(three_cys);

no_cys = find(clasif==10);
no_cys_seq=AllSeq(no_cys);
no_cys_abundances=AllOccur(no_cys);
no_cys_ntd=AllNtd(no_cys);

two_cys = find(clasif==2);
% two_cys_seq=AllSeq(two_cys);
two_cys_abundances=AllOccur(two_cys);
% two_cys_ntd=AllNtd(two_cys);

four_cys = find(clasif==4);
% four_cys_seq=AllSeq(four_cys);
four_cys_abundances=AllOccur(four_cys);
% four_cys_ntd=AllNtd(four_cys);

%% DATA ANALYSIS: if three cysteines, group by loop length

if isempty(threecys) == 1;
else
    threecys=sortrows(threecys,[1 -2]); %sorts the rows first according to first column ascent order, then second column in descent order
    SizesMat = cell2mat(threecys(:,1)); %transform to matrix
    
    %check MatLab version to use unique
    CurrVer = version('-date');
    RefVer = 'September 9, 2012'; %corresponding to 2 days before R2012b, versions afterwards changed function 'unique'
    %change date to numbers
    CurrVer = datenum(CurrVer);
    RefVer = datenum(RefVer);
    
    if CurrVer - RefVer > 0 %newer version
        [Size,Pos,~] = unique(SizesMat,'legacy'); % find where is the first of each size
    else
        [Size,Pos,~] = unique(SizesMat);
    end;

    Loops(1,1) = Size(1);
    Loops(1,2) = Pos(1);
    for iii=2:size(Size)
        Loops(iii,1) = Size(iii);
        Loops(iii,2) = Pos(iii)-Pos(iii-1);
    end;
end;

%% DATA ANALYSIS: if two cysteines, group by loop length

if isempty(twocys) == 1;
else
    
    twocys=sortrows(twocys,[1 -2]); %sorts the rows first according to first column ascent order, then second column in descent order
    SizesMat = cell2mat(twocys(:,1)); %transform to matrix

    if CurrVer - RefVer > 0 %newer version
        [SizeT,PosT,~] = unique(SizesMat,'legacy'); % find where is the first of each size
    else
        [SizeT,PosT,~] = unique(SizesMat);
    end;
    
    Loop(1,1) = SizeT(1);
    Loop(1,2) = PosT(1);
    for iii=2:size(SizeT)
        Loop(iii,1) = SizeT(iii);
        Loop(iii,2) = PosT(iii)-PosT(iii-1);
    end;
end;

%% DATA ANALYSIS: if 4 cysteines, group by loop length

if isempty(fourcys) == 1;
else 
    fourcys=sortrows(fourcys,[1 -2]); %sorts the rows first according to first column ascent order, then second column in descent order
    SizesMat = cell2mat(fourcys(:,1)); %transform to matrix

    if CurrVer - RefVer > 0 %newer version
        [SizeF,PosF,~] = unique(SizesMat,'legacy'); % find where is the first of each size
    else
        [SizeF,PosF,~] = unique(SizesMat);
    end;    
    
    Loopss(1,1) = SizeF(1);
    Loopss(1,2) = PosF(1);
    for iii=2:size(SizeF)
        Loopss(iii,1) = SizeF(iii);
        Loopss(iii,2) = PosF(iii)-PosF(iii-1);
    end;
end;

%% WRITE FILES from adaptor sequences analysis

%Printing the not-3cys file
fh = fopen(fullfile(indir,outdir,'other',[outname '_many_cys.txt']),'w');
for i=1:numel(many_cys)
    fprintf(fh, '%s\r\n', [ many_cys_seq{i} '     ' num2str(many_cys_abundances(i)) '    ' many_cys_ntd{i}]);
end;
fprintf(fh, '%d\r\n', sum(many_cys_abundances));
fclose('all');

%Printing the 4cys files:
if isempty(fourcys) == 1;
else
    fh = fopen(fullfile(indir,outdir,'4cys',[outname '_' num2str(Loopss(1,1)) '_fourcys.txt']),'w');
    % fprintf(fh, '%d\r\n', sum(threecys{,2}));
    sumcum4 = zeros();
    for v=1:PosF(1)
        fprintf(fh, '%s\r\n', [num2str(fourcys{v,3}) '     ' num2str(fourcys{v,2}) '    '  fourcys{v,4}]);
        sumcum4(1) = sumcum4(1) + fourcys{v,2};
    end;
    fprintf(fh, '%d\r\n', sumcum4(1));
    fclose('all');

    for v=2:size(SizeF)
        fh = fopen(fullfile(indir,outdir,'4cys',[outname '_' num2str(Loopss(v,1)) '_fourcys.txt']),'w'); 
        sumcum4(v) = 0;
        for iv=PosF(v-1)+1:PosF(v)
            fprintf(fh, '%s\r\n', [fourcys{iv,3} '     ' num2str(fourcys{iv,2}) '    ' fourcys{iv,4}]);
            sumcum4(v) = sumcum4(v) + fourcys{iv,2};
        end;
        fprintf(fh, '%d\r\n', sumcum4(v));
        fclose('all');
    end;
end;

%Printing the 1cys file
fh = fopen(fullfile(indir,outdir,'other',[outname '_one_cys.txt']),'w');
for i=1:numel(one_cys)
    fprintf(fh, '%s\r\n', [ one_cys_seq{i} '     ' num2str(one_cys_abundances(i)) '    '  one_cys_ntd{i}]);
end;
fprintf(fh, '%d\r\n', sum(one_cys_abundances));
fclose('all');

%Printing the 3cys files:
if isempty(threecys) == 1;
else
    fh = fopen(fullfile(indir,outdir,'3cys',[outname '_' num2str(Loops(1,1)) '_threecys.txt']),'w');
    % fprintf(fh, '%d\r\n', sum(threecys{,2}));
    sumcum3 = zeros();
    for v=1:Pos(1)
        fprintf(fh, '%s\r\n', [num2str(threecys{v,3}) '     ' num2str(threecys{v,2}) '    '  threecys{v,4}]);
        sumcum3(1) = sumcum3(1) + threecys{v,2};
    end;
    fprintf(fh, '%d\r\n', sumcum3(1));
    fclose('all');

    for v=2:size(Size)
        fh = fopen(fullfile(indir,outdir,'3cys',[outname '_' num2str(Loops(v,1)) '_threecys.txt']),'w'); 
        sumcum3(v) = 0;
        for iv=Pos(v-1)+1:Pos(v)
            fprintf(fh, '%s\r\n', [threecys{iv,3} '     ' num2str(threecys{iv,2}) '    ' threecys{iv,4}]);
            sumcum3(v) = sumcum3(v) + threecys{iv,2};
        end;
        fprintf(fh, '%d\r\n', sumcum3(v));
        fclose('all');
    end;
end;

%Printing the 2cys files:
if isempty(twocys) == 1;
else
    fh = fopen(fullfile(indir,outdir,'2cys',[outname '_' num2str(Loop(1,1)) '_twocys.txt']),'w');
    % fprintf(fh, '%d\r\n', sum(threecys{,2}));
    sumcum2 = zeros();
    for v=1:PosT(1)
        fprintf(fh, '%s\r\n', [num2str(twocys{v,3}) '     ' num2str(twocys{v,2}) '    ' twocys{v,4}]);
        sumcum2(1) = sumcum2(1) + twocys{v,2};
    end;
    fprintf(fh, '%d\r\n', sumcum2(1));
    fclose('all');

    for v=2:size(SizeT)
        fh = fopen(fullfile(indir,outdir,'2cys',[outname '_' num2str(Loop(v,1)) '_twocys.txt']),'w'); 
        sumcum2(v) = 0;
        for iv=PosT(v-1)+1:PosT(v)
            fprintf(fh, '%s\r\n', [twocys{iv,3} '     ' num2str(twocys{iv,2}) '    ' twocys{iv,4}]);
            sumcum2(v) = sumcum2(v) + twocys{iv,2};
        end;
        fprintf(fh, '%d\r\n', sumcum2(v));
        fclose('all');
    end;
end;

%Printing the 0cys file
fh = fopen(fullfile(indir,outdir,'other',[outname '_no_cys.txt']),'w');
for i=1:numel(no_cys)
    fprintf(fh, '%s\r\n', [ no_cys_seq{i} '     ' num2str(no_cys_abundances(i)) '    '  no_cys_ntd{i}]);
end;
fprintf(fh, '%d\r\n', sum(no_cys_abundances));
fclose('all');

%Printing stats files
fh = fopen(fullfile(indir,outdir,'other', [outname '_stats.txt']),'a');
fprintf(fh, 'Distribution by number of cysteine residues: \r\ndif.seq\r\n\r\n');
fprintf(fh, 'many \t %d \r\n4cys \t %d \r\n3cys \t %d \r\n2cys \t %d \r\n1cys \t %d \r\n0cys \t %d \r\ntot \t %d \r\n\r\n', [ A, B, F, C, D, E, numelall]);

fprintf(fh, 'Distribution by number of cysteine residues: \r\ntotal seq\r\n\r\n');
fprintf(fh, 'many \t %d \r\n4cys \t %d \r\n3cys \t %d \r\n2cys \t %d \r\n1cys \t %d \r\n0cys \t %d \r\ntotal seq \t %d \r\n\r\n', [ sum(many_cys_abundances), sum(four_cys_abundances), sum(three_cys_abundances), sum(two_cys_abundances), sum(one_cys_abundances), sum(no_cys_abundances), total_seq_considered]);

if isempty(threecys) == 1;
else
fprintf(fh, 'Distribution of loop lengths among the 3-cysteine clones: dif.seq. \t total seq. \r\n\r\n');
for i=1:size(Loops)
    fprintf(fh, '%s\r\n', [num2str(Loops(i,1)) '    ' num2str(Loops(i,2)) '     ' num2str(sumcum3(i))]);
end;
end;

if isempty(twocys) == 1;
else
fprintf(fh, 'Distribution of loop lengths among the 2-cysteine clones: dif.seq. \t total seq. \r\n\r\n');
for i=1:size(Loop)
    fprintf(fh, '%s\r\n', [num2str(Loop(i,1)) '    ' num2str(Loop(i,2)) '     ' num2str(sumcum2(i))]);
end;
end;

if isempty(fourcys) == 1;
else
fprintf(fh, 'Distribution of loop lengths among the 4-cysteine clones: dif.seq. \t total seq. \r\n\r\n');
for i=1:size(Loopss)
    fprintf(fh, '%s\r\n', [num2str(Loopss(i,1)) '    ' num2str(Loopss(i,2)) '     ' num2str(sumcum4(i))]);
end;
fclose('all');
end;


    






    



    









    
    