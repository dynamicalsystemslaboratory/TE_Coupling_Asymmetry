%%% This code plots percentage incorrect inference by ---------------------
%%% different measures - TE, TE_C1, TE_C2 ---------------------------------

clear all; format short; close all; 
FTsz = 17;

set(groot,'defaulttextFontName','Arial');
set(groot,'defaultLegendFontName','Arial');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

save_fig_dir = './figures/';

fpos = [100 100 350 290]; % Size and position of contour plots

% Parameters
sigw=1; beta=0;
bst=-1.0; Nb=2001; bend=1.0; db=0.001;
ast=-1.0; Na=2001; aend=1.0; da=0.001;

Ny=1; 
netw = 1; % 1: fully connected network, -1: star network

thresh2 = 1e-5;

load(['./data/Incorr_perc_TE_TEs_TEs2_Nab',num2str(Nb),'.mat']);

epsarr = [epsarr(1:3) epsarr(6:end)];
TEincorr_percarr = [TEincorr_percarr(1:3,:); TEincorr_percarr(6:end,:)];
TEsincorr_percarr = [TEsincorr_percarr(1:3,:); TEsincorr_percarr(6:end,:)];
TEs2incorr_percarr = [TEs2incorr_percarr(1:3,:); TEs2incorr_percarr(6:end,:)];

hc = zeros(size(epsarr));
colorarr = ["#FF0000", "#0000FF", "#77AC30", "#7E2F8E", "#4DBEEE", "#FF00FF", "#EDB120","#D95319"];
if (length(colorarr)<length(epsarr))
    disp('Add more colors!'); return;
end
thresh=1e-4; 
FTszlabel = FTsz;       

for ieps=1:length(epsarr)

    eps=epsarr(ieps);
    TEincorr_perc = TEincorr_percarr(ieps,:); 
    TEsincorr_perc = TEsincorr_percarr(ieps,:); 
    TEs2incorr_perc = TEs2incorr_percarr(ieps,:); 

    f1=figure(30);
    semilogx(alpharr,TEincorr_perc,'-.',...
                   'Linewidth',2,'Color',colorarr(ieps)); hold on;
    semilogx(alpharr,TEsincorr_perc,'s','Linewidth',2,...
        'Color',colorarr(ieps),'MarkerSize',10,'MarkerFaceColor',colorarr(ieps)); hold on;
    hc(ieps) = semilogx(alpharr,TEs2incorr_perc,'-','Linewidth',2,...
        'Color',colorarr(ieps)); hold on;
    
end
ylabel("Incorrect inference (%)"); 
xlabel("\alpha"); 
lh1 = legend(hc, string(epsarr)); title(lh1,'\epsilon');
f1.Position = fpos;
saveas(f1,[save_fig_dir 'incorr_inf_perc'],'epsc');

f2=figure(1);
hc1 = semilogx(alpharr,TEincorr_perc,'-.k',...
           'Linewidth',2); hold on;
hc2 = semilogx(alpharr,TEsincorr_perc,'s','Linewidth',2,...
        'Color','k','MarkerSize',10,'MarkerFaceColor','k'); hold on;
hc3 = semilogx(alpharr,TEs2incorr_perc,'-k','Linewidth',2); hold on;
lh2 = legend([hc1 hc2 hc3], ...
    'net TE','net TE^{C1}',...
    'net TE^{C2}');
f2.Position = fpos;


