function FindSeq(varargin)


%Usage:
%   FindSeq('seq','XXX')
%   FindSeq('inname','filenamt.txt','indir','path')
%   FindSeq(...,'cter','XXXX','cutoff',n)


CUTOFF = 0; %discard sequences that appeared less than CUTOFF times

%% INPUT SECTION
inname = '';  
outdir = ''; % default save directory is the same as input directory
indir = '';
cter = '';
seq = '';
% check for input variable
if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); 
    end

    % read input variables
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'inname', inname = varargin{ni+1};
            case 'seq', seq = varargin{ni+1};
            case 'indir', indir=varargin{ni+1};
            case 'cter', cter=varargin{ni+1};
            case 'cutoff', CUTOFF=varargin{ni+1};
        end
    end
end


% check whether inname was defined
if strcmp(inname,'')
    [inname,indir,~] = uigetfile('*.txt','Select text file');
else
    [~,message] = fopen(fullfile(indir, inname));
    if strcmp(message,'') == 0
        display('File not found, a dialog box will open...');
        [inname,indir,~] = uigetfile('*.txt','Select text file');
    end;
end;

if strcmp(seq,'')
    'you must enter a regular expression as a sequence';
    return
else
    name = strrep(seq,'?','_');
    name = strrep(name,'.','X');
    name = strrep(name,'*','_');
    name = strrep(name,'\','_');
    name = strrep(name,'/','_');
    name = strrep(name,'%','_');
    name = strrep(name,':','_');
    name = strrep(name,'"','_');
    name = strrep(name,'<','_');
    name = strrep(name,'>','_');
    name = strrep(name,'|','_');
    name = strrep(name,' ','_');
    outname= strcat('Seq_',name);   % default save name
    display(['Looking for motif: ' seq]);
end;

if strcmp(cter,'')
    display('considering all sequences');
else
    display(['considering all sequences whose C-ter is ' cter]);
end;

if CUTOFF == 0
    display('No minimum abundance specified');
else
    display(['Considering sequences whose abundance is  ' num2str(CUTOFF) ' or higher']);
end;


if strcmp(outdir,'')
    mkdir(indir, 'Seq');
    outdir = 'Seq';
end

%% DATA READING
%open file and read data

file = fopen(fullfile(indir, inname));
AllVar = textscan(file, '%s %d %s %*[^\n]');
fclose('all');

AllSeq = AllVar{1}; %Sequences are stored as a cell array of strings
AllOccur = AllVar{2};
AllNtd = AllVar{3};

KEEP = find(AllOccur>=CUTOFF); %I discard sequences that appeared less than 'CUTOFF' times (for naive libraries set CUTOFF to 1!!)
AllSeq = AllSeq(KEEP);
AllOccur = AllOccur(KEEP);
AllNtd = AllNtd(KEEP);

if strcmp(cter,'') == 1
else
    KEEP = ~cellfun('isempty',(strfind(AllSeq,cter)));
    AllSeq = AllSeq(KEEP);
    AllOccur = AllOccur(KEEP);
    AllNtd = AllNtd(KEEP);
end;

A = sum(AllOccur); %total sequences considered
B = numel(AllSeq); %total different sequences

clear('AllVar');

%% DATA ANALYSIS: finding sequences containing 'seq'

MatchMat = regexp(AllSeq,seq);
NoMatchMat = cellfun('isempty',MatchMat);
MatchMat = ~cellfun('isempty',MatchMat);

MatchSeq = AllSeq(MatchMat);
MatchOccur = AllOccur(MatchMat);
MatchNtd = AllNtd(MatchMat);
C = sum(MatchOccur);
D = numel(MatchOccur);

NoMatchSeq = AllSeq(NoMatchMat);
NoMatchOccur = AllOccur(NoMatchMat);
NoMatchNtd = AllNtd(NoMatchMat);
E = sum(NoMatchOccur);
F = numel(NoMatchOccur);


%% DATA PRINTING:

%Printing stats files
fh = fopen(fullfile(indir,outdir,[outname '_stats.txt']),'a');
fprintf(fh, 'Sequences considered:\r\n');
fprintf(fh, 'total \t %d \r\ndif. \t %d \r\n\r\n\r\n', [ A, B]);

fprintf(fh, 'Distribution:\r\n');
fprintf(fh, 'Match total \t %d \r\nMatch dif \t %d \r\nOther total \t %d \r\nOther dif \t %d \r\n\r\n', [ C, D, E, F]);

fclose('all');

%Printing the match file
fh = fopen(fullfile(indir,outdir,[outname '_match.txt']),'w');
for i=1:numel(MatchSeq)
    fprintf(fh, '%s\r\n', [MatchSeq{i} '     ' num2str(MatchOccur(i)) '     ' MatchNtd{i} ]);
end;
fclose('all');

%Printing the 4cys file
fh = fopen(fullfile(indir,outdir,[outname '_nomatch.txt']),'w');
for i=1:numel(NoMatchSeq)
    fprintf(fh, '%s\r\n', [NoMatchSeq{i} '     ' num2str(NoMatchOccur(i)) '     ' NoMatchNtd{i} ]);
end;
fclose('all');





