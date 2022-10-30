%%% This code calculates the TE, TE_C1, TE_C2 for -------------------------
%%% multivariate linear network systems -----------------------------------
%%% for a range of values of epsilon & alpha ------------------------------
%%% and compares with true dominant coupling direction to -----------------
%%% obtain the percentage of incorrect inferences out of the --------------
%%% entire (a,b) domain that satisfies stationarity condition -------------
%%% alpha scaled by Ny ----------------------------------------------------

clear all; format short; %clc; 
close all; 
% Parameters
sigw=1; beta=0; eps=0.02;
bst=-1.0; Nb=2001; bend=1.0; db=0.001;
ast=-1.0; Na=2001; aend=1.0; da=0.001;
thresh2 = 1e-5;

netw = 1; % 1: fully connected network, -1: star network
Nyarr = [1 2 4 7 9];
hc = zeros(size(Nyarr));
colorarr = ["#FF0000", "#0000FF", "#77AC30", "#7E2F8E", "#4DBEEE", "#EDB120", "#FF00FF"];
thresh=1e-4; % Threshold to decide if we can take the 
FTsz = 15; FTszlabel = 20;

alphparr = [0.1 0.2 0.5 1 2 5 10]; % alpha' = alpha/Ny

%% Correct sign of info flow from equations
corrIF = zeros(Na,Nb); corrIF(:,:)=nan;
ib=0;
for b=bst:db:bend
    ib = ib + 1; 
    ia=0;
    for a=ast:da:aend
        ia = ia + 1;
        if (abs(abs(b+eps)-abs(b))<thresh2) % no net X Y flow
            corrIF(ia,ib) = 0;
        elseif (abs(b+eps)>abs(b)) % X->Y
            corrIF(ia,ib) = 1;
        else                    % Y->X
            corrIF(ia,ib) = -1;
        end
    end
end

for iNy=1:length(Nyarr)

    Ny=Nyarr(iNy)

    % Covariance matrix of noise W^(n)
    Q = zeros(Ny+1,Ny+1);
    Q(1,1)=sigw^2;
    for i=2:1:Ny+1
        Q(i,i)=(sigw + (beta*eps))^2;
    end
    vecQ = reshape(Q,[size(Q,1)*size(Q,2),1]);

    TEincorr_perc = zeros(size(alphparr));
    TEsincorr_perc = zeros(size(alphparr));
    TEs2incorr_perc = zeros(size(alphparr));

    for ialph=1:length(alphparr)
        alph=alphparr(ialph)*Ny;
        Iden = eye(Ny+1,Ny+1);
        IXI= kron(Iden,Iden);

        %% Find the max a,b of the realizable domain
        ib=0; 
        rhoAarr=zeros(Na,Nb); 
        netTExtoy=zeros(Na,Nb); netTExtoy(:,:)=nan;
        netTExtoystarold=zeros(Na,Nb); netTExtoystarold(:,:)=nan;
        netTExtoystar2old=zeros(Na,Nb); netTExtoystar2old(:,:)=nan;

        for b=bst:db:bend
            ib = ib + 1; 
            ia=0;
            for a=ast:da:aend
                ia = ia + 1;
                %% Setting up the system: X^(n+1) = A X^(n) + W^(n)
                % where X^(n) = [x^(n) y1^(n) y2^(n) .... yN^(n)]
                
                % Matrix A
                if (netw > 0) %%%%%% If y1,y2,y3... are all connected with each oher by "b"
                    A = ones(Ny+1,Ny+1); 
                    A=A*b; 
                else             %%%%%% If y1,y2,y3... are disconnected wrt each oher
                    A = zeros(Ny+1,Ny+1);
                end 
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
                rhoA = max(abs(eig(A))); rhoAarr(ia,ib) = rhoA;
                rhoAXA = max(abs(eig(AXA)));
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

            end
        end

        TEtotcount=sum(~isnan(netTExtoy(:))); % total (a,b) points where rho(A)<1
        TEincorr_count = TEtotcount - numcorr(netTExtoy,corrIF,thresh,bst,db,bend,ast,da,aend);
        TEincorr_perc(ialph) = 100*TEincorr_count/TEtotcount; % percneg(ialph)

        TEincorr_countarr(iNy,ialph)=TEincorr_count;
        TEtotcountarr(iNy,ialph)=TEtotcount;
        TEincorr_percarr(iNy,ialph)=TEincorr_perc(ialph);

        TEstotcount=sum(~isnan(netTExtoystarold(:))); % total (a,b) points where rho(A)<1
        TEsincorr_count = TEstotcount - numcorr(netTExtoystarold,corrIF,thresh,bst,db,bend,ast,da,aend);
        TEsincorr_perc(ialph) = 100*TEsincorr_count/TEstotcount; % percneg(ialph)
        
        TEsincorr_countarr(iNy,ialph)=TEsincorr_count;
        TEstotcountarr(iNy,ialph)=TEstotcount;
        TEsincorr_percarr(iNy,ialph)=TEsincorr_perc(ialph);

        TEs2totcount=sum(~isnan(netTExtoystar2old(:))); % total (a,b) points where rho(A)<1
        TEs2incorr_count = TEs2totcount - numcorr(netTExtoystar2old,corrIF,thresh,bst,db,bend,ast,da,aend);
        TEs2incorr_perc(ialph) = 100*TEs2incorr_count/TEs2totcount; % percneg(ialph)
        
        TEs2incorr_countarr(iNy,ialph)=TEs2incorr_count;
        TEs2totcountarr(iNy,ialph)=TEs2totcount;
        TEs2incorr_percarr(iNy,ialph)=TEs2incorr_perc(ialph);

    end

end

save(['./data/Incorr_perc_Netw_TE_TEs_TEs2_eps0_02_Nab',num2str(Nb),'.mat'], ...
    'Nyarr', 'alphparr', 'TEincorr_percarr', 'TEsincorr_percarr', 'TEs2incorr_percarr');

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

function corr = numcorr(X,corrIF,thresh,bst,db,bend,ast,da,aend)
    ib=0;corr=0;
    for b=bst:db:bend
        ib = ib + 1; 
        ia=0;
        for a=ast:da:aend
            ia = ia + 1;
            if (sign(X(ia,ib)) == sign(corrIF(ia,ib)))
                corr = corr + 1;
            elseif ( (abs(X(ia,ib))<thresh) && (corrIF(ia,ib)==0) )
                corr = corr + 1;
            end
        end
    end

end