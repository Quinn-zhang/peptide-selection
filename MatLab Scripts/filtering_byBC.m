function filtering_byBC(varargin)

%Usage:
% filtering_byBC('inname','BCn.txt','indir','input directory','bc',BC)
% filtering_byBC('inname','BCn.txt','indir','input directory','bc',BC,'indelmut','on')

% filtering_byBC distributes the sequences in different files according to
% the barcode (PERFECT MATCH or INDELMUT MATCH if 'indelmut' is 'on')

%Output:
% - n files called 'BCn.txt', n being the number of the barcode
% - stats file, with the number of sequences assigned to the different
% barcodes (last column is no barcode assigned)


%% INPUT SECTION
indelmut = '';
inname = '';  
indir = '';

% check for input variable
if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); 
    end

    % read input variables
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'inname', inname = varargin{ni+1};
            case 'bc', BC=varargin{ni+1};
            case 'indir', indir=varargin{ni+1};
            case 'indelmut', indelmut=varargin{ni+1};
            case 'outdir', outdir=varargin{ni+1};
        end
    end
end


%% DATA READING
%open file and read data

file = fopen(fullfile(indir, inname));
AllVar = textscan(file, '%s %s %*[^\n]');
fclose('all');

AllSeq = AllVar{1}; %Sequences are stored as a cell array of strings
AllQual = AllVar{2};
clear('AllVar');

% Calculating number of barcodes
SizeBC = size(BC);


%% DATA ANALYSIS: 'fixing' in-del-mut barcodes

if strcmp(indelmut,'on') == 1
    BClength = length(BC{1});
    %Generating barcodes
    indelmutBC = cell(SizeBC(2),BClength);
    for i=1:SizeBC(2)
        indelmutBC{i,1} = ['^((.)?|' BC{i}(1) '(.)?)' BC{i}(2:BClength) ];
        for ii=2:BClength-1
            indelmutBC{i,ii} = ['^' BC{i}(1:ii-1) '((.)?|' BC{i}(ii) '(.)?)' BC{i}(ii+1:BClength)];
        end;
        indelmutBC{i,BClength} = ['^' BC{i}(1:BClength-1)];
    end;

    %Finding & replacing
    for i=1:SizeBC(2)
        for ii=1:BClength
            AllSeq = regexprep(AllSeq,indelmutBC{i,ii},BC{i});
        end;
    end; 
else
end;


%% DATA ANALYSIS: search for barcodes

%creating a matrix NxM (M barcodes), containing '1' when the barcode is present and '0' if it is not present or in a position other than 1
MposBC = zeros(size(AllSeq,1),SizeBC(2));
for i=1:SizeBC(2)
   posBC=strfind(AllSeq,BC{i}); %posBC is a column of cells
   for ii=1:size(posBC,1)
    if isempty(posBC{ii})
        posBC{ii}=0;
    elseif posBC{ii}(1)>1 %Don't consider the barcodes found out of the first position (barcode comes right after the trimmed key sequence, i.e. position 1)
        posBC{ii}=0;
    elseif posBC{ii}(1)==1 
        posBC{ii}=posBC{ii}(1); 
    end 
   end
   MposBCtemp = cell2mat(posBC);
   MposBC(:,i) = MposBCtemp;
end

%now we find which rows we need to take into the different groups - sum of all elements in a row (without barcodes means sum of elements of a row ==0)
SUMrows = sum(MposBC,2);

BCrow = cell(SizeBC(2),1);
for i=1:SizeBC(2)
    BCrow{i} = find(MposBC(:,i));
end;

BCNOTrow = find(SUMrows==0);
BCNOTnumel = numel(BCNOTrow);

%if we sum all the elements in a column, we have the number of sequences in which a certain barcode was present
SUM = sum(MposBC); 
    
    
%% WRITE FILES from adaptor sequences analysis

% Printing stats file
fh = fopen(fullfile(outdir,'statsBC.txt'),'a');
fprintf(fh, '%d ', [ SUM, BCNOTnumel ]);
fprintf(fh, '\r\n');
fclose('all');

% Printing barcode files
for i=1:SizeBC(2)
    AllBCSeq = AllSeq(BCrow{i});
    AllBCQual = AllQual(BCrow{i});
    fh = fopen(fullfile(outdir,['BC' num2str(i) '.txt']),'a');
    for ii=1:size(BCrow{i})
        fprintf(fh, '%s\r\n', [ AllBCSeq{ii} '  ' AllBCQual{ii}]);
    end
fclose('all');
end;

% Printing BCNOT file
AllBCNOTSeq = AllSeq(BCNOTrow);
AllBCNOTQual = AllQual(BCNOTrow);
fh = fopen(fullfile(outdir,'BCNOT.txt'),'a');
for i=1:size(BCNOTrow)
    fprintf(fh, '%s\r\n', [ AllBCNOTSeq{i} '  ' AllBCNOTQual{i}]);
end
fclose('all');

end




    






    



    









    
    