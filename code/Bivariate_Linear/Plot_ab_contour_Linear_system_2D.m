%%% This code calculates and plots the TE, TE_C1, TE_C2 for ---------------
%%% bivariate linear system at specific values of epsilon & alpha ---------
%%% in the entire (a,b) domain that satisfies stationarity condition ------

clear all; format short; close all; 

FTsz = 17;

set(groot,'defaulttextFontName','Arial');
set(groot,'defaultLegendFontName','Arial');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

save_fig_dir = './figures/';
fstr = 'alph5_eps0_1';

% multishade contour colormap - custom
rgb = [1 0 0; 1 0.8346 0.8346; 1.0000 0.8567 0.8567; 0.8567 0.8567 1.0000; 0.8346 0.8345 1; 0 0 1];

% Size and position of contour plots
fpos = [900 400 350 290];

% Parameters
sigw=1; beta=0;eps=0.1;
bst=-1.0; Nb=2001; bend=1.0; db=0.001; 
ast=-1.0; Na=2001; aend=1.0; da=0.001;

Nyarr = 1;
alphparr = 5; 
figure; 

for iNy=1:length(Nyarr)

    Ny=Nyarr(iNy);

    % Covariance matrix of noise W^(n)
    Q = zeros(Ny+1,Ny+1);
    Q(1,1)=sigw^2;
    for i=2:1:Ny+1
        Q(i,i)=(sigw + (beta*eps))^2;
    end
    vecQ = reshape(Q,[size(Q,1)*size(Q,2),1]);
    
    flag1 = zeros(size(alphparr));
    percneg = zeros(size(alphparr));
    
    for ialph=1:length(alphparr)
        alph=alphparr(ialph)*Ny;
        Iden = eye(Ny+1,Ny+1);
        IXI= kron(Iden,Iden);

        %% Find the max a,b of the realizable domain
        ib=0; 
        rhoAarr=zeros(Na,Nb); 
        netnormTExtoy = zeros(Na,Nb); netnormTExtoy(:,:)=nan;
        netTExtoy=zeros(Na,Nb); netTExtoy(:,:)=nan;
        netTExtoystarold=zeros(Na,Nb); netTExtoystarold(:,:)=nan;
        netTExtoystar2old=zeros(Na,Nb); netTExtoystar2old(:,:)=nan;
        netnormTExtoystar=zeros(Na,Nb); netnormTExtoystar(:,:)=nan;
        netnormTExtoystar2=zeros(Na,Nb); netnormTExtoystar2(:,:)=nan;

        for b=bst:db:bend
            ib = ib + 1; 
            ia=0;
            for a=ast:da:aend
                ia = ia + 1;
                %% Setting up the system: X^(n+1) = A X^(n) + W^(n)
                % where X^(n) = [x^(n) y1^(n) y2^(n) .... yN^(n)]
                
                % Matrix A
                A = ones(Ny+1,Ny+1); 
                A=A*b; %%%%%% If y1,y2,y3... are all connected with each oher by "b"
                A(1,1)=a;
                for j=2:1:Ny+1
                    A(1,j)=b;
                end
                for i=2:1:Ny+1
                    A(i,1)= b + eps;
                    A(i,i)= a + alph*eps;
                end
                AXA = kron(A,A);
                
                % Spectral Radius should be <1
                rhoA = max(abs(eig(A)));
                rhoAarr(ia,ib) = rhoA;
                rhoAXA = max(abs(eig(AXA)));
%                 if (max(rhoA,rhoAXA)>=1)
                IA = IXI-AXA;
                if ((rhoA>=1) || (abs(det(IA))<1e-10)) 
                    continue;
                end
                
                
                %% Calculate Covariance matrices
                
                % C^(1)(X) ---------------------->
                % C_n,n
                C11vec = (IA)\vecQ; 
                C11 = reshape(C11vec,[round(sqrt(size(C11vec,1))),round(sqrt(size(C11vec,1)))]);
                
                % C^(2)(X) ---------------------->
                % C_n,n+1
                C12 = C11*transpose(A);
                C21 = transpose(C12);
                C22 = C11;
                C2X = [C11, C12;
                       C21, C22];
            
                % C^(3)(X) ----------------------->
                D = C21*inv(C11);
                Cc = C22 - (C21*inv(C11)*C12);
                C13 = C11*transpose(D*D);
                C23 = D*C13 + Cc*transpose(D);
                C33 = D*C23 + Cc;
                C31 = transpose(C13);
                C32 = transpose(C23);
                
                C3X = [C11, C12, C13;
                       C21, C22, C23;
                       C31, C32, C33];
                
                %% Calculate TE b/w x and any y_i
                
                x_idx=1; x1_idx=(1+Ny)+1;
                i=1; yi_idx=1+i; yi1_idx=(1+Ny)+1+i;
                
                % TE_yi_x
                Hxx1 = entropy([x_idx,x1_idx],C2X); % H(x^(n),x^(n+1))
                Hxy = entropy([x_idx,yi_idx],C2X); % H(x^(n),yi^(n))
                Hxyx1 = entropy([x_idx,yi_idx,x1_idx],C2X);  % H(x^(n),yi^(n),x^(n+1))
                Hx = entropy(x_idx,C2X);  % H(x^(n))
                TEytox = Hxx1 + Hxy - Hx - Hxyx1;
                
                % TE_x_yi
                Hyy1 = entropy([yi_idx,yi1_idx],C2X); % H(yi^(n),yi^(n+1))
                Hxy = entropy([x_idx,yi_idx],C2X); % H(x^(n),yi^(n))
                Hxyy1 = entropy([x_idx,yi_idx,yi1_idx],C2X);  % H(x^(n),yi^(n),yi^(n+1))
                Hy = entropy(yi_idx,C2X);  % H(y^(n))
                TExtoy = Hyy1 + Hxy - Hy - Hxyy1;
                
                netTExtoy(ia,ib) = TExtoy - TEytox;

                Hcondx = Hxx1 - Hx;
                Hcondy = Hyy1 - Hy;
                
                netnormTExtoy(ia,ib) = TExtoy/Hcondy - TEytox/Hcondx;
                
        
                %% TE*_{x->y} = I(Xn+1;Yn+2|Yn+1,Xn)
                x_idx=1; x1_idx=(1+Ny)+1; x2_idx=2*(1+Ny)+1;
                i=1; yi_idx=1+i; yi1_idx=(1+Ny)+1+i; yi2_idx=2*(1+Ny)+1+i;
            
                % TE*_x_yi
                Hxy1y2 = entropy([x_idx,yi1_idx,yi2_idx],C3X); 
                Hxy1 = entropy([x_idx,yi1_idx],C3X); 
                Hxx1y1y2 = entropy([x_idx,x1_idx,yi1_idx,yi2_idx],C3X);  
                Hxx1y1 = entropy([x_idx,x1_idx,yi1_idx],C3X); 
                TExtoystarold = Hxy1y2 - Hxy1 - Hxx1y1y2 + Hxx1y1;
                % TE*_yi_x
                Hyx1x2 = entropy([yi_idx,x1_idx,x2_idx],C3X); 
                Hyx1 = entropy([yi_idx,x1_idx],C3X);
                Hyx1y1x2 = entropy([yi_idx,x1_idx,yi1_idx,x2_idx],C3X);
                Hyx1y1 = entropy([yi_idx,x1_idx,yi1_idx],C3X);
                TEytoxstarold = Hyx1x2 - Hyx1 - Hyx1y1x2 + Hyx1y1;
            
                netTExtoystarold(ia,ib)= TExtoystarold - TEytoxstarold;

                Hcondx2 = Hyx1x2 - Hyx1;
                Hcondy2 = Hxy1y2 - Hxy1;
                
                netnormTExtoystar(ia,ib) = TExtoystarold/Hcondy2 - TEytoxstarold/Hcondx2;

                %% TE**_{x->y} = I(Xn+1;Yn+2|Yn+1,Xn,Yn)
                x_idx=1; x1_idx=(1+Ny)+1; x2_idx=2*(1+Ny)+1;
                i=1; yi_idx=1+i; yi1_idx=(1+Ny)+1+i; yi2_idx=2*(1+Ny)+1+i;
            
                % TE**_x_yi
                Hxyy1y2 = entropy([x_idx,yi_idx,yi1_idx,yi2_idx],C3X); 
                Hxyy1 = entropy([x_idx,yi_idx,yi1_idx],C3X); 
                Hxyx1y1y2 = entropy([x_idx,yi_idx,x1_idx,yi1_idx,yi2_idx],C3X);  
                Hxyx1y1 = entropy([x_idx,yi_idx,x1_idx,yi1_idx],C3X); 
                TExtoystar2old = Hxyy1y2 - Hxyy1 - Hxyx1y1y2 + Hxyx1y1;
                % TE**_yi_x
                Hxyx1x2 = entropy([x_idx,yi_idx,x1_idx,x2_idx],C3X); 
                Hxyx1 = entropy([x_idx,yi_idx,x1_idx],C3X);
                Hxyx1y1x2 = entropy([x_idx,yi_idx,x1_idx,yi1_idx,x2_idx],C3X);
                Hxyx1y1 = entropy([x_idx,yi_idx,x1_idx,yi1_idx],C3X);
                TEytoxstar2old = Hxyx1x2 - Hxyx1 - Hxyx1y1x2 + Hxyx1y1;
            
                netTExtoystar2old(ia,ib)= TExtoystar2old - TEytoxstar2old;

                Hcondx3 = Hxyx1x2 - Hxyx1;
                Hcondy3 = Hxyy1y2 - Hxyy1;
                
                netnormTExtoystar2(ia,ib) = TExtoystar2old/Hcondy3 - TEytoxstar2old/Hcondx3;

            end
        end
        flag1(ialph)= (min(netTExtoy(:))<0); % disp('netTE < 0 at some a,b')
    
        negcount=sum(netTExtoy(:)<0);
        totcount=sum(~isnan(netTExtoy(:))); % total (a,b) points where rho(A)<1
        percneg(ialph) = 100*negcount/totcount; %percneg(ialph)
    
        negcountarr(ialph,iNy)=negcount;
        totcountarr(ialph,iNy)=totcount;
        percnegarr(ialph,iNy)=percneg(ialph);

        xlims = [-1,0.5]; ylims = [-1,1]; 
 
        [a_arr,b_arr] = meshgrid(ast:da:aend,bst:db:bend);
        
        if (abs(alph)<1e-5)
            idx1 = (netTExtoy>0.1); idx2 = (netTExtoy<-0.1);
            netTExtoy2 = netTExtoy; 
            netTExtoy2(idx1) = 0.0; netTExtoy2(idx2) = 0.0;
    
            f1=figure; 
            contourf(a_arr',b_arr',netTExtoy2,20,'LineColor','none'); hold on;
            [M,c]=contour(a_arr',b_arr',netTExtoy,[0.0,0.0],'ShowText','on');
            c.LineWidth=2;c.LineColor='k';
            colormap(multigradient(rgb, [1 19 20 22 24 41]));
            xlabel('a');
            ylabel('b');%title(['netTE_{x->y_i}, N=' num2str(Ny+1) ', alpha=' num2str(alph)]);
            xlim(xlims); ylim(ylims);
            c = colorbar; title(c,'net TE_{X \rightarrow Y}','Fontsize',FTsz);
            c.Title.Position = [-61 192 0];
            f1.Position = fpos;
            caxis([-0.1 0.1]);
            saveas(f1,[save_fig_dir 'netTE_' fstr],'epsc');
        else
            f1=figure; 
            contourf(a_arr',b_arr',netTExtoy,20,'LineColor','none'); hold on;
            [M,c]=contour(a_arr',b_arr',netTExtoy,[0.0,0.0],'ShowText','on');
            c.LineWidth=2;c.LineColor='k';
            colormap(multigradient(rgb, [1 19 20 22 24 41]));
            xlabel('a');
            ylabel('b');
            xlim(xlims); ylim(ylims);
            c = colorbar; title(c,'net TE_{X \rightarrow Y}','Fontsize',FTsz);
            c.Title.Position = [-65 192 0]; caxis([-0.2 0.2]);
            f1.Position = fpos;
            saveas(f1,[save_fig_dir 'netTE_' fstr],'epsc');
        end

        f2=figure; 
        contourf(a_arr',b_arr',netnormTExtoy,20,'LineColor','none'); colorbar; hold on;
        [M,c]=contour(a_arr',b_arr',netnormTExtoy,[0.0,0.0],'ShowText','on');
        c.LineWidth=2;c.LineColor='k';
        colormap(multigradient(rgb, [1 19 20 22 24 41]));
        xlabel('a');
        ylabel('b');
        xlim(xlims); ylim(ylims);
        c = colorbar; title(c,"net TE^{'}_{X \rightarrow Y}",'Fontsize',FTsz);
        caxis([-0.05 0.05]);
        c.Title.Position = [-61 192 0];
        f2.Position = fpos;
        saveas(f2,[save_fig_dir 'netnormTE_' fstr],'epsc');

        f3=figure; 
        contourf(a_arr',b_arr',netTExtoystarold,20,'LineColor','none'); colorbar; hold on;
        [M,c]=contour(a_arr',b_arr',netTExtoystarold,[0.0,0.0],'ShowText','on');
        c.LineWidth=2;c.LineColor='k';
        colormap(multigradient(rgb, [1 19 20 22 24 41]));
        xlabel('a');
        ylabel('b');
        xlim(xlims); ylim(ylims);
        c = colorbar; title(c,'net TE^{C1}_{X \rightarrow Y}','Fontsize',FTsz);
        c.Title.Position = [-65 192 0]; caxis([-0.06 0.06]);
        f3.Position = fpos;
        saveas(f3,[save_fig_dir 'netTE_C1_' fstr],'epsc');

        f5=figure; 
        contourf(a_arr',b_arr',netTExtoystar2old,20,'LineColor','none'); colorbar; hold on;
        [M,c]=contour(a_arr',b_arr',netTExtoystar2old,[0.0,0.0],'ShowText','on');
        c.LineWidth=2;c.LineColor='k';
        colormap(multigradient(rgb, [1 19 20 22 24 41]));
        xlabel('a');
        ylabel('b');       
        xlim(xlims); ylim(ylims);
        c = colorbar; title(c,'net TE^{C2}_{X \rightarrow Y}','Fontsize',FTsz);
        c.Title.Position = [-65 192 0];
        f5.Position = fpos; caxis([-0.06 0.06]);
        saveas(f5,[save_fig_dir 'netTE_C2_' fstr],'epsc');

        return
    
    end

end

%% Functions

% Calc entropy/joint entropy ----->
function H = entropy(idx,C)
    nidx = length(idx);
    if (nidx==1)
        cov = C(idx(1),idx(1));
    elseif (nidx==2)
        cov = [C(idx(1),idx(1)), C(idx(1),idx(2));
               C(idx(2),idx(1)), C(idx(2),idx(2))];
    elseif (nidx==3)
        cov = [C(idx(1),idx(1)), C(idx(1),idx(2)), C(idx(1),idx(3));
               C(idx(2),idx(1)), C(idx(2),idx(2)), C(idx(2),idx(3));
               C(idx(3),idx(1)), C(idx(3),idx(2)), C(idx(3),idx(3))];
    elseif (nidx==4)
        cov = [C(idx(1),idx(1)), C(idx(1),idx(2)), C(idx(1),idx(3)), C(idx(1),idx(4));
               C(idx(2),idx(1)), C(idx(2),idx(2)), C(idx(2),idx(3)), C(idx(2),idx(4));
               C(idx(3),idx(1)), C(idx(3),idx(2)), C(idx(3),idx(3)), C(idx(3),idx(4));
               C(idx(4),idx(1)), C(idx(4),idx(2)), C(idx(4),idx(3)), C(idx(4),idx(4))];
    elseif (nidx==5)
        cov = [C(idx(1),idx(1)), C(idx(1),idx(2)), C(idx(1),idx(3)), C(idx(1),idx(4)), C(idx(1),idx(5));
               C(idx(2),idx(1)), C(idx(2),idx(2)), C(idx(2),idx(3)), C(idx(2),idx(4)), C(idx(2),idx(5));
               C(idx(3),idx(1)), C(idx(3),idx(2)), C(idx(3),idx(3)), C(idx(3),idx(4)), C(idx(3),idx(5));
               C(idx(4),idx(1)), C(idx(4),idx(2)), C(idx(4),idx(3)), C(idx(4),idx(4)), C(idx(4),idx(5));
               C(idx(5),idx(1)), C(idx(5),idx(2)), C(idx(5),idx(3)), C(idx(5),idx(4)), C(idx(5),idx(5))];
    else
        disp('Error: we can only calculate H for upto 5 timeseries!');
    end
    H = 0.5*log( ((2*pi*exp(1))^nidx)*det(cov) );
end