function [a,b,c] = MakeCDJVFilter(request,DN)

% MakeCDJVFilter --  Set up filters for CDJV Wavelet Transform by reading
% in the text files in the directory CDJV_Filters
%  Usage
%    [a,b,c] = MakeCDJVFilter(request,degree)
%  Inputs
%    request  string: 'HighPass', 'LowPass', 'Precondition', 'Postcondition'
%    degree   integer: 2 or 3 (number of vanishing moments)
%  Outputs
%    a,b,c    filter, left edge filter, right edge filter
%             ('HighPass', 'LowPass')
%    a        conditioning matrix ('Precondition', 'Postcondition')
%
%  Description
%    CDJV have developed an algorithm for wavelets on the interval which
%    preserves the orthogonality, vanishing moments, smoothness, and compact
%    support of Daubechies wavelets on the line.
%
%    The algorithm for wavelets on the interval of CDJV involves four objects
%    not present in the usual periodized algorithm: right edge filters, left
%    edge filters, and pre- and post- conditioning operators.
%
%    These objects are supplied by appropriate requests to MakeCDJVFilter.
%
%  References
%    Cohen, Daubechies, Jawerth and Vial, 1992.
%

cur_dir =  fileparts(mfilename('fullpath'));

if (strcmp(request,'HighPass')||strcmp(request,'LowPass'))
    
    %Internal Filter
    filename = fullfile(cur_dir, 'CDJV_Filters', sprintf('FILT%d.txt',  DN));
    fid = fopen(filename);
    B=textscan(fid,'%s',inf,'delimiter','\n');
    fclose(fid);

    data = B{1,1};
    str = sprintf(' %s,',data{:});
    data_array = sprintf('[%s]', str);

    internal_filters = eval(data_array);
    
    s_f = sum(internal_filters);
    if (s_f - sqrt(2) > 0.00001)
        internal_filters =  (sqrt(2)/s_f) *internal_filters;
    end


    %%%Left Filter
    filename = fullfile(cur_dir, 'CDJV_Filters', sprintf('LFILT%d.txt',  DN));
    fid = fopen(filename);
    A = fscanf(fid, '%e %e', [2 inf]);
    fclose(fid);
    left_LP = A(1,:);
    left_HP = A(2,:);

    %%%Right Filter

    filename = fullfile(cur_dir, 'CDJV_Filters', sprintf('RFILT%d.txt',  DN));
    fid = fopen(filename);
    A = fscanf(fid, '%e %e', [2 inf]);
    fclose(fid);
    right_LP = A(1,:);
    right_HP = A(2,:);



    if strcmp(request,'HighPass'),
        LEHI = zeros(DN,3*DN-1); k=1; l=1;
        for j = DN+1:2:3*DN-1
            LEHI(k, 1:j) = left_HP(l:l+j-1);
            l=l+j; k=k+1;
        end

        REHI = zeros(DN,3*DN-1); k=1; l=1;
        for j = DN+1:2:3*DN-1
            REHI(k, 1:j) = right_HP(l:l+j-1);
            l=l+j; k=k+1;
        end
        
        a = reverse(MirrorFilt(internal_filters)); b = LEHI; c = REHI;
    end
    if strcmp(request, 'LowPass')
        LELO = zeros(DN,3*DN-1); k=1; l=1;
        for j = DN+1:2:3*DN-1
            LELO(k, 1:j) = left_LP(l:l+j-1);
            l=l+j; k=k+1;
        end

        RELO = zeros(DN,3*DN-1); k=1; l=1;
        for j = DN+1:2:3*DN-1
            RELO(k, 1:j) = right_LP(l:l+j-1);
            l=l+j; k=k+1;
        end
        
        a = internal_filters; b = LELO; c = RELO;
    end
end
    if (strcmp(request,'PreCondition')||strcmp(request,'PostCondition'))

    filename = fullfile(cur_dir, 'CDJV_Filters', sprintf('LCOND%d.txt',  DN));
    fid = fopen(filename);
    %fid = fopen(sprintf('%s/CDJV_Filters/LCOND%d.txt', cur_dir, DN));
    A = fscanf(fid, '%f %f', [2 inf]);
    fclose(fid);
    left_precond = A(1,:);
    left_postcond = A(2,:);


    LPREMAT = reshape(left_precond, sqrt(length(left_precond)), sqrt(length(left_precond)));
    LPOSTMAT = reshape(left_postcond, sqrt(length(left_postcond)), sqrt(length(left_postcond)));


    filename = fullfile(cur_dir, 'CDJV_Filters', sprintf('RCOND%d.txt',  DN));
    fid = fopen(filename);
    %fid = fopen(sprintf('%s/CDJV_Filters/RCOND%d.txt', cur_dir, DN));
    A = fscanf(fid, '%f %f', [2 inf]);
    fclose(fid);
    right_precond = A(1,:);
    right_postcond = A(2,:);

    RPREMAT = reshape(right_precond, sqrt(length(right_precond)), sqrt(length(right_precond)));
    RPOSTMAT = reshape(right_postcond, sqrt(length(right_postcond)), sqrt(length(right_postcond)));


    if strcmp(request,'PreCondition'),
            a = (LPREMAT)'; b = rot90(RPREMAT,2)'; c = [];
    end

    if strcmp(request,'PostCondition'),
        a = (LPOSTMAT)'; b = rot90(RPOSTMAT,2)'; c = [];
    end
    
    end

end


% This file was created by Clarice Poon (cmhsp2@cam.ac.uk)
% It is based on the MakeCDJVFilter.m file from Wavelab 850
% The user may directly replace MakeCDJVFilter.m from Wavelab 850 if he/she
% require CDJV filters of 4 vanishing moments or more.

% 2018: This file have been further edited by Vegard Antun, to fix some path issues.

