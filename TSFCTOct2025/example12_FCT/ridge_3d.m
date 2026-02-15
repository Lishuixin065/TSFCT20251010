function [IF, CR] = ridge_3d(Tx, n, fa, ca, gd)
% Input
%   Tx : 3D space, for example, Tx(chirp-rate, frequency, time)
%    n : Number of ridges to extract
%   fa : Maximum allowable frequency variation for ridge tracing
%   ca : Maximum allowable chirp-rate variation for ridge tracing
%   gd : Adjustable parameter for mode reconstruction

% OUTPUTS:  
%   CR: Matrix of instantaneous frequencies extracted along the ridges
%   IF : Matrix of instantaneous chirprates extracted along the ridges

Et = (abs(Tx) + eps).^2; % Lower the energy size

[Nc, Nf, Nt] =size(Tx); 
IF = zeros( Nt,n); % Instantaneous frequency
CR = zeros( Nt,n); % Instantaneous chirp-rate

e = 10e-8; % Threshold

% The algorithm
for r = 1:n
   
    [~, idx1] = max(abs(Et(:)));
%最大值和线性索引
    [idc, idf, idt]= ind2sub(size(Et),idx1);
    IF( idt,r) = idf;
    CR(idt,r) = idc;

    % Forward search
    for c = idt + 1:Nt
        a =max(1, idc - ca): min(Nc, idc + ca);
        b = max(1, idf - fa): min(Nf, idf + fa); 
        E1 = squeeze(Et(a, b, c));      
        [M, I] = max(E1(:));

        if M > e
            [idcs, idfs] = ind2sub([length(a), length(b)], I);            
            idc = a(idcs);
            idf = b(idfs);
          
            if idc > Nc
             idc = Nc;
            end
             if idf > Nf
             idf = Nf;
            end

            IF( c,r) = idf;
            CR( c,r) = idc;

        else
            IF( c,r) = IF( c - 1,r);
            CR( c,r) = CR( c - 1,r);
        end
    end

    % Backward search
    idf = IF( idt,r);
    idc = CR( idt,r);

    for c = idt - 1:-1:1
        b = max(1, idf - fa): min(Nf, idf + fa);
        a = max(1, idc - ca): min(Nc, idc + ca);

        E2 = squeeze(Et(a, b, c));
        [M, I] = max(E2(:));       
        if M > e
            [idcs, idfs] = ind2sub([length(a), length(b)], I);
            idf = b(idfs);
            idc = a(idcs);
            IF( c,r) = idf;
            CR( c,r) = idc;
             if idf > Nf
             idf = Nf;
            end

            if idc > Nc
             idc = Nc;
            end

        else
            IF( c,r) = IF( c + 1,r);
            CR( c,r) = CR(c + 1,r);
        end
    end

% Reconstruct the mode & Remove the curve
for c = 1:Nt
   

        c_index_low = max(CR(c,r)-gd,1);
        c_index_up = min(CR(c,r)+gd,Nc);
        f_index_low = max(IF(c,r)-gd,1);
        f_index_up = min(IF(c,r)+gd,Nf);
        Et( c_index_low:c_index_up,f_index_low:f_index_up,c)=0;
end
   
end



end