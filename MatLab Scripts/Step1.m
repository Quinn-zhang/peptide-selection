function Step1(varargin)
% Usage: 
% Step1
% Step1('inname','filename.fastq','indir','path')
% Step1(...,'bc',{'AAAAAA';'TTTTTT';'GGGGGG';'...'})
% Step1(...,'indelmut','on')
% 
% input: (optional)

%   Step1('inname','filename.fastq','indir','path') indicates the file name and path to the folder where it is located. If not specified, a dialog box to choose the file will open.
%   Step1('indelmut','on') allows one insertion, deletion or mutation in the barcodes. If not specified, it is off.
%   Step1('bc',{'AAAAAA','TTTTTT','GGGGGG','...'}) indicates the barcodes used. They must be separated by comma. If not specified, it uses the ones described in Rentero Rebollo et al. 2014.

% It reads the initial filename.fastq file and generates files containing reads according to their barcode (named BC1.txt, BC2.txt... BCNOT.txt), and saves them in a separate folder within the input folder (named "filename_BC"). BCNOT.txt contains reads whose barcode did not correspond to any of the identified barcodes.

% output:
%   A new folder called "filename_BC" with a series of files named BC1.txt, BC2.txt... BCN.txt, containing the reads corresponding to the first, second,... nth barcode respectively. Reads whose barcode did not correspond to any of the identified barcodes are stored in BCNOT.txt.
%   The chip-specific code (shown in the command window in MatLab)
%   BC_stats.txt file, containing information about how many reads where found per barcode.


tic

%% Reading input
inname = '';  
chunk = 200000;     % default value 200000, MUST BE LESS THAN SEQUENCES
outname= 'Fraction';   % default save name for temporary files, will have a number appended to it
indir = '';
indelmut = '';
BC = '';


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
            case 'indelmut', indelmut = varargin{ni+1};
            case 'bc', BC=varargin{ni+1};
        end
    end
end

% check whether inname was defined
if strcmp(inname,'')
    [inname,indir,~] = uigetfile('*.fastq','Select .fastq file');
else
    [~,message] = fopen(fullfile(indir, inname));
    if strcmp(message,'') == 0
        display('File not found, a dialog box will open...');
        [inname,indir,~] = uigetfile('*.fastq','Select .fastq file');
    end;
end;

if strcmp(BC,'') %barcodes described in Rentero Rebollo I. et al. 2014
    BC{1} = 'GCATAG';
    BC{2} = 'CGTATC';
    BC{3} = 'ATCGCA';
    BC{4} = 'ACGATA';
    BC{5} = 'AGACTC';
    BC{6} = 'GATACA';
    BC{7} = 'CATCTC';
    BC{8} = 'GTTCAG';
    BC{9} = 'TACCAG';
    BC{10} = 'ATGGAG';
    BC{11} = 'AGTTAC';
    BC{12} = 'GGTGAA';
end;

% check if BC was entered correctly
if iscell(BC) == 0
    display(['invalid BC format input, please indicate it as: Step1(' char(39) 'bc' char(39) ',{' char(39) 'AAAAAA' char(39) ',' char(39) 'TTTTTT' char(39) ',' char(39) 'GGGGGG' char(39) ',' char(39) '...' char(39) '})']);
    return
end;
BC_sizes = cellfun('length',BC);
if max(BC_sizes) == min(BC_sizes)
else
    display('All barcodes must have the same length');
    return
end;
display(BC);

% check if indelmut is active
if strcmp(indelmut,'on')
    display('Considering barcodes having 1 in-del-mut');
else
    display('Condireing barcodes with perfect match');
end

%make outdir
outdir = fullfile(indir, [inname(1:end-6) '_RAW']); % default save directory for fraction files
if ~isdir(outdir)
    mkdir(outdir);    
end;

BCdir = fullfile(indir, [inname(1:end-6) '_BC']); %default save directory for BC files
if ~isdir(BCdir)
    mkdir(BCdir);    
end;


%% Step 1.1: breaking the file in smaller temporary files called FractionN.txt

FID = fopen(fullfile(indir, inname));

% look for the code
FirstLine = textscan(FID, '%s %*[^\n]', 1);
FirstLine = FirstLine{1}{1};
End = strfind(FirstLine,':');
End = End(1);
code = textscan(FirstLine(1:End-1), '%s');
code = code{1}{1};

%it will display the code and number of barcodes in the command window
display(['Chip-specific code = ' code]);
SizeBC = size(BC);
if SizeBC(1) > 1
    display('Invalid barcode input')
    return
end;
display(['Number of barcodes = ' num2str(SizeBC(2))]);

% Breaking the file
currL=cell(4*chunk,1);
i=0;
ii=0;
n=0;

while ~feof(FID) %feof() test for end-of-file
    i=i+1;
    ii = ii+1;
    currL{ii} = fgetl(FID); %fgetl() returns the next line of the file removing newline character
                            %It separates the file in a list of lines,
                            %stored in a cell array
    
    if mod(i,4*chunk) == 0 %mod(X,X) = 0, when we get to the 4* number of sequences
        A1 = strncmp(code,currL,length(code)); %find all strings that start from @NAME,
         %currL is an array, compares each element; A1 is an array too
         %the number must correspond to the number of characters to compare
         %gives 1 if TRUE, 0 if FALSE
                                          
              
        Slines = find(A1)+1;               %shift the line number
        AllSeq = currL(Slines);
            %Makes an array with sequences
        
        Qlines = find(A1)+3;               %shift the line number
        AllQual = currL(Qlines); 
            %Makes an array with qualities
        
        n=n+1; % counter for the files
        fh = fopen(fullfile(outdir, [outname num2str(n) '.txt']), 'w');
        
        %write the data
        for jj=1:size(AllSeq,1)
          fprintf(fh, '%s\r\n', [ AllSeq{jj} '  ' AllQual{jj}]);
        end;
        
        disp('reading file...');
        
        % clear the variable, restart the indices
        currL=cell(4*chunk,1);
        fclose(fh);
        ii=0;
    end;
    
    if feof(FID) == 1 %mod(X,X) = 0, when we get to the 4* number of sequences

        A1 = strncmp(code,currL,length(code)); %find all strings that start from @NAME,
         %currL is an array, compares each element; A1 is an array too
         %the number must correspond to the number of characters to compare
         %gives 1 if TRUE, 0 if FALSE
                                          
              
        Slines = find(A1)+1;               %shift the line number
        AllSeq = currL(Slines);
        
        Qlines = find(A1)+3;               %shift the line number
        AllQual = currL(Qlines); 

        n=n+1; % counter for the files
        fh = fopen(fullfile(outdir, [outname num2str(n) '.txt']), 'w');
        
        %write the data
        for jj=1:size(AllSeq,1)
          fprintf(fh, '%s\r\n', [ AllSeq{jj} '  ' AllQual{jj}]);
        end;

        disp('reading file...');
        
        % clear the variable, restart the indices
        currL=cell(4*chunk,1);
        fclose(fh);
        ii=0;
    end;


end;

fclose('all');

display(['File read in ' num2str(toc) ' sec']);


%% Step 1.2: assigning barcodes

tic

if strcmp(indelmut,'on') == 1
    for i=1:n
        filtering_byBC('inname',[outname num2str(i) '.txt'],'indir',outdir,'outdir',BCdir,'bc',BC,'indelmut',indelmut);
        display('Sorting by barcode...');
    end;
else
    for i=1:n
        filtering_byBC('inname',[outname num2str(i) '.txt'],'indir',outdir,'outdir',BCdir,'bc',BC);
        display('Sorting by barcode...')
    end;
end;

% Remove temporary files
rmdir(outdir,'s');

% Sum the number of reads of each fraction in the stats file
SizeBC = size(BC);

scanformat = '';
for n=1:SizeBC(2)+1
    scanformat = [scanformat '%d '];
end;
scanformat = [scanformat '%*[^\n]'];
    
fh = fopen(fullfile(BCdir,'statsBC.txt'));
    Number_reads = textscan(fh, scanformat);
fclose('all');

Number_reads = cell2mat(Number_reads);
Number_reads = sum(Number_reads,1);

fh = fopen(fullfile(BCdir,'statsBC.txt'),'w');
for i=1:SizeBC(2)
    fprintf(fh, '%s\n',[ 'BC' num2str(i) '  ' BC{i} '  ' num2str(Number_reads(i))]);
end;
    fprintf(fh, '%s\n',[ 'BCNOT  ' num2str(Number_reads(SizeBC(2)+1))]);
fclose('all');

display(['Barcodes asigned in ' num2str(toc) ' sec']);
    
    


