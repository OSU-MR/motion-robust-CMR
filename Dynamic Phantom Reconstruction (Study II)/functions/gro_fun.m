function[PEInd, samp] = gro_fun(param)
% Author: Rizwan Ahmad (ahmad.46@osu.edu)

% Output
% ky(i,e) = PEInd(i,e), where 'e' is encoding and 'i' is the order of acquisition
% t(i) = FRInd(i)
% samp = binary mask in ky-t-encoding domain

n   = param.n;   % Number of phase encoding (PE) lines per frame
FR  = param.FR;   % Frames
N   = param.PE;  % Size of of PE grid
E   = param.E;    % Number of encoding, E=1 for cine, E=2 for flow
tau = param.tau;
PF  = param.PF;
s   = param.s;
a   = param.alph;
dsp = param.dsp;
M   = param.M;

gr = (1+sqrt(5))/2; % golden ratio
gr = 1/(gr+tau-1); % golden angle, sqrt(2) works equally well
% R  = N/(n+PF); % Acceleration
% s  = max(1, s*(R).^(1/3)); % This is to adapt the value of's' based on acceleration. This line can be commented out.


%% Size of the smaller pseudo-grid which after stretching gives the true grid size
Ns = ceil(N * 1/s); % Size of shrunk PE grid
k = (N/2-Ns/2)/((Ns/2)^a); %(1/2*(n-Ns)/max(tmp)); % location specific displacement

samp  = zeros(N, FR, E); % sampling on k-t grid
PEInd = zeros((n-PF)*FR, E); % The ordered sequence of PE indices and encoding
FRInd = zeros((n-PF)*FR, 1); % The ordered sequence of frame indices

% figure;
v0 = (1/2+1e-10:Ns/(n+PF):Ns+1/2-1e-10); % Start with uniform sampling for each frame
for e=1:E
    v0 = v0 + 1/E*Ns/(n+PF); % Start with uniform sampling for each frame
    kk=E+1-e;
    for j=1:FR
        v = rem((v0 + (j-1)*Ns/(n+PF)*gr)-1, Ns) + 1; % In each frame, shift by golden shift of PES/TR*gr
        v = v - Ns.*(v>=(Ns+0.5));

        if rem(N,2)==0 % if even, shift by 1/2 pixel
            vC = v - k*sign((Ns/2+1/2)-v).*(abs((Ns/2+1/2)-v)).^a + (N-Ns)/2 + 1/2;%(ctrn - ctrnS);
            vC = vC - N.*(vC>=(N+0.5));
        elseif rem(N,2)==1 % if odd don't shift
            vC = v - k*sign((Ns/2+1/2)-v).*(abs((Ns/2+1/2)-v)).^a + (N-Ns)/2;%(ctrn - ctrnS);
        end
        vC = round(sort(vC));
        vC(1:PF) = [];

        if rem(j,2) == 1
            PEInd((j-1)*n + 1 : j*n, e) = vC;
        else
            PEInd((j-1)*n + 1 : j*n, e) = flip(vC);
        end
        FRInd((j-1)*n + 1 : j*n) = j;

        samp(vC, j, e) = samp(vC, j, e) + kk;
    end
end

if dsp == 1
    tiFont = 20; % title font
    axFont = 14; % axis font
    laFont = 18; % label font
    figure;
    subNum = E + double((E/2) >= 1); % Number of subplots
    tiledlayout(1,subNum,'TileSpacing','compact', 'Padding', 'compact')
    for e = 1 : subNum
        nexttile;
%         subplot(1,subNum,e);
        if e == 1
            imagesc(samp(:,:,e),[0,E]); axis('image'); colormap(hot);
            set(gca, 'FontSize', axFont, 'FontName','times');
            xlabel('$t$', 'FontSize', laFont,'Interpreter','latex'); 
            ylabel('$k_y$', 'FontSize', laFont,'Interpreter','latex');   
            if E==1 % do nothing
            elseif E>1, title(['Encoding ' num2str(e)], 'FontSize', tiFont,'Interpreter','latex'); end
        elseif e < (E+1)
            imagesc(samp(:,:,e),[0,E]); axis('image'); colormap(hot);
            set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]);
            xlabel('$t$', 'FontSize', laFont,'Interpreter','latex');    
            title(['Encoding ' num2str(e)], 'FontSize', tiFont,'Interpreter','latex');
        else
            imagesc(max(samp,[],3)); axis('image'); colormap(hot);
            set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]);
            xlabel('$t$', 'FontSize', laFont,'Interpreter','latex');  
            title('Superimposed', 'FontSize', tiFont,'Interpreter','latex');
        end
    end
    set(gcf,'color','w','units','points','position',[10,10,450,350]); %export_fig('gro',gcf,'-m4','-png');
    
    figure;
    %len = min([length(PEInd), 120]);
    len=length(PEInd);
    plot(1:len, PEInd(1:len,1), '-'); hold on;
    plot(1:len, PEInd(1:len,1), '.','MarkerSize', 12);
    set(gca, 'FontSize', axFont, 'FontName','times');
    xlabel('Acquisition Order', 'FontSize', laFont,'Interpreter','latex');
    ylabel('$k_y$', 'FontSize', laFont,'Interpreter','latex');
    title('PE Index vs. Acquisition Order', 'FontSize', tiFont,'Interpreter','latex');
    set(gcf,'color','w','units','points','position',[10,10,600,350]); %export_fig('gro-index',gcf,'-m4','-png');
end

samp = logical(samp); % convert to logical before returning
