function Step2(varargin)

%Usage:
% Step2
% Step2('inname','filename.txt','indir','path')
% Step2(...,'badmax',n,'q','quality value')
% Step2(...,'uplimit',m,'downlimit',o,'midlimit',p)
% Step2(...,'start','ATGGC','end','GCTGAAAC')
% Step2(...,'fixerr',n,'cter','XXX')
% Step2(...,'translateall','on')
%
% Input (optional):
%   Step2('inname','filename.txt','indir','path') indicates the file name and path to the folder where it is located. If not specified, a dialog box to choose the file will open.
%       Step2 can read the output file from Step1 (BC1.txt, BC2.txt...)
%   Step2('badmax',n) where n is the maximum number of bases below the
%       quality threshold allowed. If not specified, badmax = 3
%   Step2('q',Q) where Q indicates the quality threshold (18 for Q18, 20
%       for Q20, etc...). If not specified, Q = 18.
%   Step2('uplimit',m,'downlimit',o,'midlimit',p) specifies the maximum
%       (m) and minimum (o) peptide length (in residues). Additionally, an
%       internal limit (=<) can be indicated (p).
%   Step2('start',NNNNN,'end',NNNNN), specify start and end of the DNA
%       region of interest. If not specified, it ises the ones described in
%       this manuscript.
%   Step2('fixerr',n,'cter','XXX') allows to fix sequencing errors: it merges together
%       sequencing with only 1 or 2 differences in the DNA sequence. 
%       It only does so with the "n" most abundant sequences. For more
%       information, see fixingerrors.m. 'cter' indicates a constant region
%       in the C-terminus of the peptide and helps fixing errors correctly
%       but it is not necessary.
%
% Output:
%   Command window: information about parameters used for quality filter
%
%   QF_BCn folder containing intermediate quality files:
%       QF_BCn_longGOOD.txt and QF_BC3_longBAD.txt contain good and bad
%           quality reads, respectively, whose length is longer than the
%           midlimit specified
%       QF_BCn_shortGOOD.txt and QF_BC3_shorBAD.txt, IDEM with reads whose
%           length is shorter or equal to the midlimt specified.
%       QF_BCn_toolong.txt and QF_BCn_tooshort.txt contain reads whose
%           length was longer or shorter than the maximum or minimum length
%           specified, respectively.
%       QF_BCn_NOLIM.txt contain sequences for which either the start or
%           the end of the DNA region of interest was not found.
%       QF_BCn_stats contain the number of reads assigned to each file
%
%   Translation_BCn folder containing translation files:
%       Translated_BCn_longGOOD.txt contain translated DNA sequences from
%           good quality reads whose length is longer than the midlimit
%           specified
%       Translated_BCn_shortGOOD.txt contain translated DNA sequences from
%           good quality reads whose length is shorter than the midlimit
%           specified
%       Translated_BCn_stats contain the information about the total number
%           of reads and the number of different reads.
%       fixerrTranslated_BCn...txt contains translated DNA sequences after
%           fixing sequencing errors (if this option was activated)
%
%       Additional files folder
%           Translated_toolong.txt and QF_BCn_tooshort.txt contain reads whose
%               length was longer or shorter than the maximum or minimum length
%               specified, respectively.
%           Translated_BCn_NOLIM.txt contain sequences for which either the start or
%               the end of the DNA region of interest was not found.
%           Translated_BCn_longBAD.txt contain translated DNA sequences from
%               bad quality reads whose length is longer than the midlimit
%               specified
%           Translated_BCn_shortBAD.txt contain translated DNA sequences from
%               bad quality reads whose length is shorter than the midlimit
%               specified
%   
%       




%% INPUT SECTION
inname = '';  
indir = '';
outdir = '';
badmax = '';
q = '';
pepstart = '';
adaptorR = '';
SHORTLIMIT = '';
UPlimit = '';
DOWNlimit = '';
fixerr = '';
cter = '';
translateall = '';
keepQF = '';

% check for input variable
if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); 
    end

    % read input variables
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'inname', inname = varargin{ni+1};
            case 'indir', indir=varargin{ni+1};
            case 'badmax', badmax=varargin{ni+1};
            case 'q', q=varargin{ni+1};
            case 'midlimit', SHORTLIMIT=varargin{ni+1}*3+1;
            case 'uplimit', UPlimit=varargin{ni+1}*3+1;
            case 'downlimit', DOWNlimit=varargin{ni+1}*3-1;
            case 'start', pepstart=varargin{ni+1};
            case 'end', adaptorR=varargin{ni+1};
            case 'fixerr', fixerr=varargin{ni+1};
            case 'cter', cter=varargin{ni+1};
            case 'translateall', translateall=varargin{ni+1};
            case 'keepqf', keepQF=varargin{ni+1};
        end
    end
end

% check whether inname was defined
if strcmp(inname,'')
    [inname,indir,~] = uigetfile('*.txt','Select BCn.txt file');
else
    [~,message] = fopen(fullfile(indir, inname));
    if strcmp(message,'') == 0
        display('File not found, a dialog box will open...');
        [inname,indir,~] = uigetfile('*.fastq','Select .fastq file');
    end;
end;

% check quality parameters
if strcmp(badmax,'') == 1
    badmax = 3;
end;
display(['Maximum number of bases below quality accepted = ' num2str(badmax)]);

if strcmp(q,'') == 1
    q = 18;
end;
Q = char(q+33);
p = 10^(-q/10);
display(['Quality threshold = ' num2str(q) ' , meaning p = ' num2str(p,3)]);

% check start and end of the peptide sequence
if strcmp(pepstart,'') == 1
    pepstart = 'ATGGC';
    display('Peptide start not specified, using default start ATGGC');
else
    display(['Peptide start specifcied: ' pepstart])
end;
if strcmp(adaptorR,'') == 1
    adaptorR = 'GCTGAAAC';
    display('Peptide end not specified, using default end GCTGAAAC')
else
    display(['Peptide end specified: ' adaptorR])
end;

% check length parameters
if strcmp(DOWNlimit,'') == 1
    DOWNlimit = 0;
    display('No minimum length of the peptide');
else
    display(['Minimum length of the peptide in bases: ' num2str(DOWNlimit)]);
end;
if strcmp(UPlimit,'') == 1
    UPlimit = Inf;
    display('No maximum length of the peptide');
else
    display(['Maximum length of the peptide in bases: ' num2str(UPlimit)]);
end;
if strcmp(SHORTLIMIT,'') == 1
    SHORTLIMIT = DOWNlimit;
    display('No intermediate limit');
else
    display(['Intermediate limit in bases (<=): ' num2str(SHORTLIMIT)]);
end;

% check fixerr
if strcmp(fixerr,'') == 1
    display('No fixing sequencing errors');
elseif strcmp(fixerr,'') == 0
    if isnumeric(fixerr) == 0
        display('ERROR: fixerr should be a number > 0')
        return
    else
        display(['Fixing sequencing errors: correcting sequencing errors of the ' num2str(fixerr) ' most abundant sequences']);
    end;
end;


%% Step 2.1: Quality filter
tic

%% DATA READING
%open file and read data

file = fopen(fullfile(indir, inname));
AllVar = textscan(file, '%s %s %*[^\n]');
fclose('all');

AllSeq = AllVar{1}; %Sequences are stored as a cell array of strings
AllQual = AllVar{2};
clear('AllVar');


%% DATA ANALYSIS
%Search for adapter rev
posR = strfind(AllSeq, adaptorR); %posR = position adaptorR, but these are cell arrays
posF = strfind(AllSeq, pepstart); %where the peptide starts


%posR is a cell array, if the adaptor sequence was not found, it
%leaves an empty cell. Transform empty cells to zeroes:
for ii=1:size(posR,1)
    if isempty(posR{ii})
        posR{ii}=0;
    elseif numel(posR{ii})>1 %it could in principle appear several times
        posR{ii}=max(posR{ii}); %use the maximum value in case it appears in the random region too
    end  
end

for ii=1:size(posF,1)
    if isempty(posF{ii})
        posF{ii}=0;
    elseif numel(posF{ii})>1 %it could in principle appear several times
        posF{ii}=min(posF{ii}); %use the minimum value in case it appears in the random region too
    end  
end


%now transform them into a matrix
MposR = cell2mat(posR);
MposF = cell2mat(posF);

%now look for quality lower than the threshold in the region of interest
clasif = zeros(numel(AllSeq),1);
for i=1:size(AllSeq)
    if MposF(i)==0
        clasif(i) = 4;
    elseif MposR(i)==0
        clasif(i) = 4;
    else
        fragmentlength = MposR(i) - MposF(i);
        if fragmentlength > UPlimit
            clasif(i) = 5; %too long sequences
        elseif fragmentlength < DOWNlimit
            clasif(i) = 6; % too short sequences     
        elseif fragmentlength <= SHORTLIMIT
        fragment = AllQual{i}(MposF(i):MposR(i)); %region of interest in the quality data
        bad = find(fragment<Q); %finds positions with low quality
            if numel(bad) > badmax;
            clasif(i) = 12; %sequences shorter than the middle limit but low quality
            else
            clasif(i) = 13; %sequences shorter than the middle limit and good quality
            end;
        else
        fragment = AllQual{i}(MposF(i):MposR(i)); %region of interest in the quality data
        bad = find(fragment<Q); %finds positions with low quality
            if numel(bad) > badmax;
            clasif(i) = 2; %sequences longer than the middle limit but low quality
            else
            clasif(i) = 3; %sequences longer than the middle limit and good quality
            end;
        end;
    end;
end;
good=find(clasif==3);
shortgood=find(clasif==13);
shortbad=find(clasif==12);
lowqual=find(clasif==2);
nolimits= find(clasif==4);
toolong= find(clasif==5);
tooshort=find(clasif==6);


   
%% WRITE FILES from adaptor sequences analysis

inname = strrep(inname, '.txt', '');
outname = ['QF_' inname];   % default save name

if strcmp(outdir,'')
    outdir = ['QF_' inname];
    mkdir(indir, outdir);
end;

fh = fopen(fullfile(indir, outdir,[outname '_stats.txt']),'a');
fprintf(fh, '\r\n\r\nPARAMETERS \r\n %s\r\n' , ['Q' num2str(q) '   Max.bases ' num2str(badmax)]);
fprintf(fh, '\nLIMITS \n %s \n', ['min. ' num2str(DOWNlimit) '   max. ' num2str(UPlimit) '   middle ' num2str(SHORTLIMIT)]);
fprintf(fh, '\n NUMBER OF READS \n');
fprintf(fh, 'longgood  \tshortgood \tlongbad \tshortbad \tnolimits \ttotal \r\n%d\t %d\t %d\t %d\t %d\t %d\r\n', [ numel(good), numel(shortgood), numel(lowqual), numel(shortbad), numel(nolimits) numel(AllSeq) ]);
fprintf(fh, 'toolong  \ttooshort \r\n%d\t %d\t \r\n', [ numel(toolong), numel(tooshort)]);

fclose('all');

GOODSeq = AllSeq(good);
GOODMposR = MposR(good);
GOODMposF = MposF(good);
fh = fopen(fullfile(indir, outdir,[outname '_longGOOD.txt']),'a');
for i=1:size(GOODSeq)
    fprintf(fh, '%s\r\n', GOODSeq{i}(GOODMposF(i):GOODMposR(i))); %I only store the info of the wanted region
end
fclose('all');

GOODShortSeq = AllSeq(shortgood);
GOODShortMposR = MposR(shortgood);
GOODShortMposF = MposF(shortgood);
fh = fopen(fullfile(indir, outdir,[outname '_shortGOOD.txt']),'a');
for i=1:size(GOODShortSeq)
    fprintf(fh, '%s\r\n', GOODShortSeq{i}(GOODShortMposF(i):GOODShortMposR(i))); %I only store the info of the wanted region
end
fclose('all');

BADShortSeq = AllSeq(shortbad);
BADShortMposR = MposR(shortbad);
BADShortMposF = MposF(shortbad);
fh = fopen(fullfile(indir, outdir,[outname '_shortBAD.txt']),'a');
for i=1:size(BADShortSeq)
    fprintf(fh, '%s\r\n', BADShortSeq{i}(BADShortMposF(i):BADShortMposR(i))); %I only store the info of the wanted region
end
fclose('all');

LOWQUALSeq = AllSeq(lowqual);
LOWQUALposR = MposR(lowqual);
LOWQUALposF = MposF(lowqual);
fh = fopen(fullfile(indir, outdir,[outname '_longBAD.txt']),'a');
for i=1:size(LOWQUALSeq)
    fprintf(fh, '%s\r\n', LOWQUALSeq{i}(LOWQUALposF(i):LOWQUALposR(i))); %I store the info of the wanted region
end
fclose('all');

NOLIMSeq = AllSeq(nolimits);
fh = fopen(fullfile(indir, outdir,[outname '_NOLIM.txt']),'a');
for i=1:size(NOLIMSeq)
    fprintf(fh, '%s\r\n', NOLIMSeq{i}); %I store the whole sequence
end
fclose('all');

toolongSeq = AllSeq(toolong);
toolongMposR = MposR(toolong);
toolongMposF = MposF(toolong);
fh = fopen(fullfile(indir, outdir,[outname '_toolong.txt']),'a');
for i=1:size(toolongSeq)
    fprintf(fh, '%s\r\n', toolongSeq{i}(toolongMposF(i):toolongMposR(i))); %I only store the info of the wanted region
end
fclose('all');

tooshortSeq = AllSeq(tooshort);
tooshortMposR = MposR(tooshort);
tooshortMposF = MposF(tooshort);
fh = fopen(fullfile(indir, outdir,[outname '_tooshort.txt']),'a');
for i=1:size(tooshortSeq)
    fprintf(fh, '%s\r\n', tooshortSeq{i}(tooshortMposF(i):tooshortMposR(i))); %I only store the info of the wanted region
end;
fclose('all');

display(['Quality filter completed in ' num2str(toc) ' sec']);

%If no middle limit was indicated, rename the files and delete short:
if SHORTLIMIT == DOWNlimit
    %current folder = QF
    movefile(fullfile(indir,outdir,[outname '_longGOOD.txt']),fullfile(indir,outdir,[outname '_GOOD.txt']));
    movefile(fullfile(indir,outdir,[outname '_longBAD.txt']),fullfile(indir,outdir,[outname '_BAD.txt']));
    delete(fullfile(indir,outdir,[outname '_shortGOOD.txt']));
    delete(fullfile(indir,outdir,[outname '_shortBAD.txt']));
end;

%% Step 2.2: Translation
outdir = fullfile(indir, outdir);
transdir = fullfile(indir, ['Translation_' inname]); %default save directory is translation
if ~isdir(transdir)
    mkdir(transdir);    
end;

% If all need to be translated an additional files folder is created
if strcmp(translateall,'on') == 1
    adddir = fullfile(transdir, 'additional files');
    if ~isdir(adddir)
        mkdir(adddir);    
    end;
end;

if SHORTLIMIT == DOWNlimit
    Translate('inname',[outname '_GOOD.txt'],'indir',outdir,'fixerr',fixerr,'outdir',transdir,'cter',cter);
    if strcmp(translateall,'on') == 1
        Translate('inname',[outname '_BAD.txt'],'indir',outdir,'fixerr',fixerr,'outdir',adddir,'cter',cter);        
    end;
else
    Translate('inname',[outname '_longGOOD.txt'],'indir',outdir,'fixerr',fixerr,'outdir',transdir,'cter',cter);
    Translate('inname',[outname '_shortGOOD.txt'],'indir',outdir,'fixerr',fixerr,'outdir',transdir,'cter',cter);
    if strcmp(translateall,'on') == 1
        Translate('inname',[outname '_longBAD.txt'],'indir',outdir,'outdir',adddir,'cter',cter);
        Translate('inname',[outname '_shortBAD.txt'],'indir',outdir,'outdir',adddir,'cter',cter);   
    end;
end;

if strcmp(translateall,'on') == 1
    Translate('inname',[outname '_toolong.txt'],'indir',outdir,'outdir',adddir,'cter',cter);
    Translate('inname',[outname '_tooshort.txt'],'indir',outdir,'outdir',adddir,'cter',cter);
end;
%Translate('inname',[outname '_NOLIM.txt'],'indir',outdir,'outdir',adddir,'cter','NOLIM');
    
% Unless specified otherwise, remocing QF folder
if strcmp(keepQF,'on') == 1
else
    rmdir(outdir,'s')
end;





    



    









    
    