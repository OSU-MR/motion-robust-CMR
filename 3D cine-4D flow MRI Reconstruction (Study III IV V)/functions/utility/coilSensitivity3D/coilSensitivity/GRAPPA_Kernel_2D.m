

function [Coef, Data_Acq, Data_Mis] = GRAPPA_Kernel_2D(Ref_Img, option)
% function [Coef, Data_Acq, Data_Mis] = GRAPPA_Kernel_2D(Ref_Img, option)
% Acc = 4;
% option.KernelSize = [4 5];
% option.KernelPattern = [-Acc : Acc : 2*Acc];
% option.OutPattern = [1:Acc-1];
% tic
s_0 = size(Ref_Img);
Coef = 0;

% make sure Kernel Size Consistent
KernelSize = option.KernelSize;
Pattern_pe = option.KernelPattern;
OutPattern = option.OutPattern;

if mod(KernelSize(2), 2) == 1
    Start_fe = (KernelSize(2)-1)/2 ;
    End_fe = (KernelSize(2)-1)/2 ;
else
    Start_fe = (KernelSize(2))/2 ;
    End_fe = (KernelSize(2))/2 ;
end

Pattern_fe = -Start_fe:End_fe;

Start_pe = -min(Pattern_pe) ;
End_pe = max(Pattern_pe) ;


CoefSize = prod(KernelSize)*s_0(3);

MisSize = length(OutPattern)*s_0(3) ;
%SamplingPattern = option.SamplingPattern ;

if length(size(KernelSize))~= 2, 
    return, 
else
    Coef = zeros( CoefSize, MisSize );
end

%{'1', toc}, tic,
counter  = 0;
if length(s_0) == 3
    Data_Acq = zeros( (s_0(1) - Start_fe - End_fe)*(s_0(2) - Start_pe - End_pe ), prod(KernelSize)*s_0(3), 'single' );
    Data_Acq = complex( Data_Acq, Data_Acq );
    Data_Mis = zeros( (s_0(1) - Start_fe - End_fe)*(s_0(2) - Start_pe - End_pe ), length(OutPattern)*s_0(3), 'single' );
    Data_Mis = complex( Data_Mis, Data_Mis );
    for index_fe = Start_fe+1:s_0(1)-End_fe
        for index_pe = Start_pe+1:s_0(2)-End_pe
            counter  = counter  + 1 ;
            temp = Ref_Img( index_fe+Pattern_fe, index_pe+Pattern_pe, :) ;
            %keyboard
            Data_Acq(counter, :) = temp(:);
            %size(Data_Mis(counter, :)), size(Ref_Img( index_fe, index_pe+OutPattern, :))
            Data_Mis(counter, :) = reshape( Ref_Img( index_fe, index_pe+OutPattern, :), 1, MisSize);
        end
    end
elseif length(s_0) == 4
    Data_Acq = zeros( (s_0(1) - Start_fe - End_fe)*(s_0(2) - Start_pe - End_pe )*s_0(4), prod(KernelSize)*s_0(3), 'single' );
    Data_Acq = complex( Data_Acq, Data_Acq );
    Data_Mis = zeros( (s_0(1) - Start_fe - End_fe)*(s_0(2) - Start_pe - End_pe )*s_0(4), length(OutPattern)*s_0(3), 'single' );
    Data_Mis = complex( Data_Mis, Data_Mis );
    for index_fr = 1:s_0(4)
        for index_fe = Start_fe+1:s_0(1)-End_fe
            for index_pe = Start_pe+1:s_0(2)-End_pe
                counter  = counter  + 1 ;
                temp = Ref_Img( index_fe+Pattern_fe, index_pe+Pattern_pe, :, index_fr) ;
                %size(temp), size(Data_Acq(counter, :))
                %keyboard
                Data_Acq(counter, :) = temp(:);
                %size(Data_Mis(counter, :)), size(Ref_Img( index_fe, index_pe+OutPattern, :))
                Data_Mis(counter, :) = reshape( Ref_Img( index_fe, index_pe+OutPattern, :, index_fr), 1, MisSize);
            end
        end
    end
else
    'Error! The Ref image MUST be 3D or 4D!'
    Coef = 0;
    Data_Acq = 0;
    Data_Mis = 0;
    return
end

counter  = counter;
%{'2',toc}, tic, 
s_A = size(Data_Acq); 
Asquare = (Data_Acq'*Data_Acq); 
Tr_A = 1e-5*trace(Asquare)/CoefSize ; %%% Original Ding's code has 1e-4
% 0.0001*Tr_A/CoefSize;
%cond(Asquare)
%max((Asquare(:))), min((Asquare(:))), % test code Ding 2011-07-26
%E = eig(Asquare)/1; {'Max', max(E), 'Min', min(E), 'Tr', Tr_A, 'Thresh', 0.0001*Tr_A/CoefSize }, % test code Ding 2012-10-10
%figure, semilogy(E, 'o'), title('Eigenvalue'), grid on
%figure, hist(E, 32)  , title('Histogram')
%{'3',toc}, tic, 
%Coef = inv(Asquare + 0.0001*Tr_A/CoefSize*eye(CoefSize) ) * (Data_Acq' * Data_Mis ) ; 
if isfield(option, 'IcePAT')||isfield(option, 'Threshold')
    if isfield(option, 'IcePAT') && option.IcePAT
        Coef = (Asquare + Tr_A*eye(CoefSize) ) \ (Data_Acq' * Data_Mis ) ;
    elseif isfield(option, 'Threshold') && (~isfield(option, 'Truncation'))
        %{'NewThresh', s_A(1)*option.Threshold^2}, test code
        Coef = (Asquare + 1*s_A(1)*option.Threshold^2*eye(CoefSize) ) \ (Data_Acq' * Data_Mis ) ;%'Threshold',
    elseif isfield(option, 'Threshold') && isfield(option, 'Truncation')
        Coef = Truncate_Pinv(Asquare, option.Threshold^2 ) * (Data_Acq' * Data_Mis ) ;
        {'SVD', 'Threshold', option.Threshold^2}
    end
else % IcePAT is the default
    Coef = (Asquare + Tr_A*eye(CoefSize) ) \ (Data_Acq' * Data_Mis ) ;
%     {'Default', 'Threshold', Tr_A}
end
%{'4',toc},
%figure(1), semilogy(1:s_A(2), eig(Asquare) ), %hold on
%{'condition #', 'Orig', cond(Asquare)}


% % add pre-conditioner
% temp = 0;
% for i = 1:s_A(1)
%     temp = max(abs(Data_Acq(i,:)));
%     if temp > 0
%         Data_Acq(i,:) = Data_Acq(i,:)/temp;
%         Data_Mis(i,:) = Data_Mis(i,:)/temp;
%     end
% end
% 
% %preconM = diag(precon/max(precon)); {min(preconM(:)), max(preconM(:))}
% %Data_Acq = preconM*Data_Acq;
% Asquare = (Data_Acq'*Data_Acq); 
% Tr_A = trace(Asquare) ;
% 0.0001*Tr_A/CoefSize;
% Coef = inv(Asquare + 0.0001*Tr_A/CoefSize*eye(CoefSize) ) * (Data_Acq' * Data_Mis ) ;
% {'condition #', 'pre-cond-ed', cond( Asquare )}
% figure(1), semilogy(1:s_A(2), eig(Asquare), 'o'), hold off
% 









