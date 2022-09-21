function Translate(varargin)

%Usage:

%   Translate('inname','filename.txt,'indir,'path','outdir','path')
%   Translate(...,'cter','XXX')
%   Translate(...,'fixerr',n)

% Input:

%   Translate('inname','filename.txt,'indir,'path','outdir','path'), where
%   filename.txt is the name of the file (QF_BCn_xxx.txt file, containing
%   nucleotide sequences), and its folder 'indir' is specified. The folder
%   where the output files will be generated is specified as 'outdir'.

%   Translate(...,'cter','XXX') indicates a constant peptide sequence at the
%   end of the random region, such as the linker to the phage, allowing to
%   identify frame-shifts. For example for the libraries described in this 
%   study it would be 'GGSG' for Library A and B, and 'SHS' for 3x3 and 4x4

%   Translate(...,'fixerr',n) indicates how many different positions are
%   allowed to be considered a sequencing error



%% INPUT SECTION
inname = '';  
outdir = ''; % default save directory is the same as input directory
indir = '';
cter = '';
fixerr = '';

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
            case 'cter', cter=varargin{ni+1};
            case 'fixerr', fixerr=varargin{ni+1};
        end
    end
end

%check if parameters was defined

if strcmp(inname,'')
    inname = uigetfile;
end

if strcmp(indir,'')
    indir = uigetdir;
end

if strcmp(outdir,'')
    outdir = indir;
end

%% DATA READING
%open file and read data

tic

file = fopen(fullfile(indir, inname));
AllVar = textscan(file, '%s %*[^\n]');
fclose('all');

AllSeq = AllVar{1}; %Sequences are stored as a cell array of strings
clear('AllVar');


%% DATA ANALYSIS: Search for duplicates
numelall = numel(AllSeq);

%sort the sequences
Sorted = sort(AllSeq);
%find unique sequences and the position of the first one
[seq,M,~] = unique(Sorted,'first');
numelunique = numel(seq);
M2 = [M(2:end); numel(Sorted)+1];
Occur = M2 - M;
[Occur,IX] = sort(Occur, 'descend');

%make an array with sequences and occurences
UniqueArraySorted = seq(IX);


%% DATA ANALYSIS: translation
GeneticCode = geneticcode(1);
AmberQ = setfield(GeneticCode,'TAG','Q'); %change amber codon to glutamine
aaseq = nt2aa(UniqueArraySorted,'GeneticCode',AmberQ);

   
%% WRITE FILES from adaptor sequences analysis

inname = regexprep(inname,'.txt','');
inname = regexprep(inname,'QF_','');
outname= ['Translated_' inname];   % default save name
fh = fopen(fullfile(outdir,'Translated_stats.txt'),'a');
fprintf(fh, '%s\n different\t max.abundance\t total\t \r\n', inname);
fprintf(fh, ' \t %d\t %d\t %d\r\n', [ numelunique, max(Occur), numelall ]);
fclose('all');

fh = fopen(fullfile(outdir,[outname '.txt']),'w');
for i=1:size(UniqueArraySorted)
    fprintf(fh, '%s\r\n', [aaseq{i} '    ' num2str(Occur(i)) '  ' UniqueArraySorted{i} ]);
end
fclose('all');

display(['Translation completed in ' num2str(toc) ' sec']);


%% OPTIONAL: removing sequencing errors

if strcmp(fixerr,'') == 1
else
    fixingerrors('inname',[ outname '.txt'],'indir',outdir, 'limit', fixerr);
end;

    






    



    









    
    