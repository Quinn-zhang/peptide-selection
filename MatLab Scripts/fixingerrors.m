function fixingerrors(varargin)

%Usage:

% fixingerrors('inname','filename.txt','indir','path','outdir','path')
% fixingerrors(..., 'nb_mut',n)
% fixingerrors(..., 'cter','XXX')
% fixingerrors(..., 'limit',m)

%Input:

% fixingerrors('inname','filename.txt','indir','path'),
% where filename.txt within the folder specified as 'indir' is the name of the file (output of Step2.m) having the
% format: peptide sequence - abundance - nucleotide sequence. If it is not
% indicated, a dialog box will open

% fixingerrors(..., 'nb_mut',n) specifies the number of different positions
% allowed to merge two nucleotide sequences as the same parent sequence


% fixingerrors(..., 'Cter','XXX') specifies the constant residues at the
% Cter of the random region. When two sequences are merged, it will take
% the most abundant as correct, unless:
%       - it is less than 4x the less abundant, and
%       - the most abundant does not have the constant residues at the
%       Cterminus and the less abundant has them.
% If these two conditions are fulfilled, the less abundant sequence will be
% considered correct.

% fixingerrors(..., 'limit',m) specifies how many sequences will be checked
% for mutations. If not specified, it checks the top200 abundat ones.
% fixingerrors(..., 'limit','all') will fix all sequences.



%% INPUT SECTION
tic
nb_mut = 2.5;
Cter = '';
limit = 200;
inname = '';

% check for input variable
if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); 
    end;

    % read input variables
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'inname', inname = varargin{ni+1};
            case 'indir', indir=varargin{ni+1};
            case 'nb_mut', nb_mut=varargin{ni+1}+0.5;
            case 'cter', Cter=varargin{ni+1};
            case 'limit', limit=varargin{ni+1};
        end;
    end;
end;

% check whether inname was defined
if strcmp(inname,'')
    [inname,indir,~] = uigetfile('*.txt','Select text file with this format: peptide seq - abundance - nucleotide seq');
else
    [~,message] = fopen(fullfile(indir, inname));
    if strcmp(message,'') == 0
        display('File not found, a dialog box will open...');
        [inname,indir,~] = uigetfile('*.txt','Select text file');
    end;
end;

%make outdir
outdir = fullfile(indir, ['correctiondata' inname]); % default save directory for fraction files
if ~isdir(outdir)
    mkdir(outdir);    
end;


%% DATA READING
%open file and read data

file = fopen(fullfile(indir, inname));
AllVar = textscan(file, '%s %d %s %*[^\n]');
fclose('all');

AllSeq = AllVar{1}; %Sequences are stored as a cell array of strings
AllOccur = AllVar{2};
AllNtd = AllVar{3};

A=sum(AllOccur);

clear('AllVar');

% Do not continue if it is an empty file
if numel(AllOccur) == 0
    return
end;


%% DATA ANALYSIS: fixing errors


seq=cell(2,1);
total = numel(AllNtd);

if strcmp(limit,'all')
    limit = total;
elseif limit > total
    limit = total;
end;

display(['Fixing errors of top ' num2str(limit) ' abundant sequences']);

for i=1:limit(1)
    display(['Looking for sequencing errors of sequence ' num2str(i) ]);
    InitialAbun = AllOccur(i);
    if strcmp(AllSeq{i},'') == 1
        continue %if this was removed already, skip it
    else
        nb_ntd = numel(AllNtd{i});
        for ii=i+1:total
            if strcmp(AllSeq{ii},'') == 1
                continue %if this was removed already, skip it
            else
                seq{1}=AllNtd{i};
                seq{2}=AllNtd{ii};
                dist=nb_ntd*seqpdist(seq,'method','p-distance');
                if dist<nb_mut
                    
                    % I check difference in abundances
                    if AllOccur(i) > 4*AllOccur(ii) % if it is more than 4x the other clone ok
                    
                        % I store the correction data in a file
                        fh = fopen(fullfile(outdir, 'Correctiondata.txt'),'a');
                        fprintf(fh, '%s\r\n', [num2str(AllOccur(i)) '  ' AllNtd{i} '   ' AllSeq{i}]);
                        fprintf(fh, '%s\r\n', [num2str(AllOccur(ii)) '   ' AllNtd{ii} '   ' AllSeq{ii}]);
                        fprintf(fh, '%d\r\n \r\n', round(dist));
                        fclose('all');

                        % Keep the first clone
                        AllOccur(i)= AllOccur(ii)+ AllOccur(i);
                        AllOccur(ii)=0;
                        AllNtd{ii}='';
                        AllSeq{ii}='';
                    
                    else % the most abundant clone is less than 4x the other clone
                        
                        % I store the correction data in another file
                        fh = fopen(fullfile(outdir, 'Conflict.txt'),'a');
                        fprintf(fh, '%s\r\n', [num2str(AllOccur(i)) '  ' AllNtd{i} '   ' AllSeq{i}]);
                        fprintf(fh, '%s\r\n', [num2str(AllOccur(ii)) '   ' AllNtd{ii} '   ' AllSeq{ii}]);
                        fprintf(fh, '%d\r\n \r\n', round(dist));
                        fclose('all');

                        % Check if there is a framshift in the first clone and
                        % not in the second one, then accept the second one
                        if isempty(strfind(AllSeq{i},Cter)) == 1 %frameshift
                            if isempty(strfind(AllSeq{ii},Cter)) == 1 %frameshift also in the second too - keep first
                                AllOccur(i)= AllOccur(ii)+ AllOccur(i);
                                AllOccur(ii)=0;
                                AllNtd{ii}='';
                                AllSeq{ii}='';
                            elseif isempty(strfind(AllSeq{ii},Cter)) == 0 %not framshift in the second
                                    AllSeq{i} = AllSeq{ii};
                                    AllOccur(i)= AllOccur(ii)+ AllOccur(i);
                                    AllNtd{i} = AllNtd{ii};

                                    AllOccur(ii)=0;
                                    AllNtd{ii}='';
                                    AllSeq{ii}='';
                             end;
                        else %no frameshif (or no "Cter" specified) - keep first
                            AllOccur(i)= AllOccur(ii)+ AllOccur(i);
                            AllOccur(ii)=0;
                            AllNtd{ii}='';
                            AllSeq{ii}='';
                        end;
                    end;
                end;
            end;
        end;
    end;
    fh = fopen(fullfile(outdir, 'ErrorRates.txt'),'a');
    fprintf(fh, '%s\r\n', [num2str(num2str((AllOccur(i) - InitialAbun)*100/AllOccur(i))) '    ' num2str(InitialAbun) '    ' num2str(AllOccur(i))]);
    fclose('all');                    
end;

% Sort the resulting data
[AllOccur,IX] = sort(AllOccur, 'descend');
AllSeq = AllSeq(IX);
AllNtd = AllNtd(IX);

B=sum(AllOccur);

% Check
display(['Number of sequences before = ' num2str(A)]);
display(['Number of sequences after = ' num2str(B)]);

display(['Fixing completed in ' num2str(toc) ' sec']);
   
%% WRITE FILES

outname = ['fixerr' inname];
fh = fopen(fullfile(indir,outname),'w');
for i=1:size(AllNtd)
    if AllOccur(i) == 0
        return
    else
    fprintf(fh, '%s\r\n', [AllSeq{i} '    ' num2str(AllOccur(i)) '  ' AllNtd{i} ]);
    end;
end
fclose('all');


    






    



    









    
    