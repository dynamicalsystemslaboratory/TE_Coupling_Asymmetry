%%% This code plots percentage incorrect inference by ---------------------
%%% different measures - TE, TE_C1, TE_C2 ---------------------------------
%%% for multi-node networks -----------------------------------------------

clear all; format short; close all; 
FTsz = 22;

set(groot,'defaulttextFontName','Arial');
set(groot,'defaultLegendFontName','Arial');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

save_fig_dir = './figures4/';

fpos = [100 100 540 420]; % Size and position of contour plots

% Parameters
sigw=1; beta=0;
bst=-1.0; Nb=2001; bend=1.0; db=0.001;
ast=-1.0; Na=2001; aend=1.0; da=0.001;

Nyarr = [1 2 4 7 9]; 
netw = 1; % 1: fully connected network, -1: star network

thresh2 = 1e-5;

load(['./data/Incorr_perc_Netw_TE_TEs_TEs2_eps0_02_Nab',num2str(Nb),'.mat']);

hc = zeros(size(Nyarr));
colorarr = ["#FF0000", "#0000FF", "#77AC30", "#7E2F8E", "#4DBEEE", "#FF00FF", "#EDB120","#D95319"];
if (length(colorarr)<length(Nyarr))
    disp('Add more colors!'); return;
end
thresh=1e-4; % Threshold to decide if we can take the 
FTsz = 15; FTszlabel = 20;        

for iNy=1:length(Nyarr)

    Ny=Nyarr(iNy);

    TEincorr_perc = TEincorr_percarr(iNy,:); 
    TEsincorr_perc = TEsincorr_percarr(iNy,:); 
    TEs2incorr_perc = TEs2incorr_percarr(iNy,:); 

    f1=figure(30);
    semilogx(alphparr,TEincorr_perc,'-.',...
                   'Linewidth',2,'Color',colorarr(iNy)); hold on;
    semilogx(alphparr,TEsincorr_perc,'s','Linewidth',2,...
        'Color',colorarr(iNy),'MarkerSize',10,'MarkerFaceColor',colorarr(iNy)); hold on;
    hc(iNy) = semilogx(alphparr,TEs2incorr_perc,'-','Linewidth',2,...
        'Color',colorarr(iNy)); hold on;
    
end
ylabel("Incorrect inference (%)"); %,'FontSize',FTszlabel,'Interpreter','latex');%,'FontSize',FTszlabel-3);
xlabel("\alpha'"); %xlabel("$\alpha$",'FontSize',FTszlabel,'Interpreter','latex');
lh1 = legend(hc, string(Nyarr+1)); title(lh1,'N');

figure(1);
hc1 = semilogx(alphparr,TEincorr_perc,'-.k',...
           'Linewidth',2); hold on;
hc2 = semilogx(alphparr,TEsincorr_perc,'s','Linewidth',2,...
        'Color','k','MarkerSize',10,'MarkerFaceColor','k'); hold on;
hc3 = semilogx(alphparr,TEs2incorr_perc,'-k','Linewidth',2); hold on;
lh2 = legend([hc1 hc2 hc3], ...
    'net TE','net TE^{C1}',...
    'net TE^{C2}');
saveas(gcf,[save_fig_dir 'incorr_inf_perc_netw'],'epsc');


