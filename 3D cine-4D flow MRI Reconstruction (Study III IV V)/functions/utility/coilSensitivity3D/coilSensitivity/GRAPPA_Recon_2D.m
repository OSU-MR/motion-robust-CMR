

function Reduced_Kspace = GRAPPA_Recon_2D(Reduced_Kspace, Coef, option)

s_0 = size(Reduced_Kspace);

KernelSize = option.KernelSize;
Pattern_pe = option.KernelPattern;
OutPattern = option.OutPattern;

Start_fe = (KernelSize(2)-1)/2 ;
End_fe = (KernelSize(2)-1)/2 ;

Pattern_fe = -Start_fe:End_fe;

Start_pe = -min(Pattern_pe) ;
End_pe = max(Pattern_pe) ;


CoefSize = prod(KernelSize)*s_0(3);

MisSize = length(OutPattern)*s_0(3) ;
OutSize = length(OutPattern);
SamplingPattern = double(option.SamplingPattern) ;

if length(size(KernelSize))~= 2, 
    return, 
end

Data_Acq = zeros( s_0(1)*length(SamplingPattern), prod(KernelSize)*s_0(3), 'single' );
Data_Acq = complex( Data_Acq, Data_Acq );

counter  = 0;
for index_fe = 1:s_0(1)
    for index_pe = 1:length(SamplingPattern)
        counter  = counter  + 1 ;
        k_fe = mod(index_fe+Pattern_fe-1, s_0(1) ) + 1 ;
        k_pe = mod(SamplingPattern(index_pe) + Pattern_pe - 1, s_0(2) ) + 1;
        temp = Reduced_Kspace( k_fe, k_pe, :) ;
        Data_Acq(counter, :) = temp(:);
%        Data_Mis(counter, :) = Ref_Img( index_fe, SamplingPattern(index_pe) + OutPattern, :);
    end
end

% Reconstruction
Data_Mis = single(Data_Acq * Coef);

% Reshape to fill k-space
counter = 0;
for index_fe = 1:s_0(1)
    for index_pe = 1:length(SamplingPattern)
        counter  = counter  + 1 ;
        k_pe = mod( SamplingPattern(index_pe) + OutPattern -1, s_0(2)) + 1;
        Reduced_Kspace( index_fe, k_pe, :) = single(reshape( Data_Mis(counter, :), OutSize, s_0(3) ));
    end
end



