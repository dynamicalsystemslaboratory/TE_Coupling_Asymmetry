%%% This code calculates and plots the TE, TE_C1, TE_C2 for ---------------
%%% bivariate linear system for a given set of epsilon,a,b values ---------
%%% for a range of values of alpha ----------------------------------------

clear all; format short; close all; 
FTsz = 17; 

set(groot,'defaulttextFontName','Arial');
set(groot,'defaultLegendFontName','Arial');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

save_fig_dir = './figures/';
fstr = 'a0_15_b0_45_eps0_1';

% Parameters
sigw=1; beta=0;
eps=0.1;
Ny = 1;
alpharr = -5:0.1:5;  
dalph= 4; 
xlimits = [alpharr(1) alpharr(end)];
a=0.15; b=0.45;

bst=0; bend=1.0; db=0.001; Nb=ceil((bend-bst)/db)+1;
ast=0; aend=1.0; da=0.001; Na=ceil((aend-ast)/da)+1;

FTsz = 15; FTszlabel = 20;        

Iden = eye(Ny+1,Ny+1);
IXI= kron(Iden,Iden);

for ialph = 1:length(alpharr)
    alph = alpharr(ialph); %alphparr(ialph)*Ny;

    %% Find the max a,b of the realizable domain
    ib=0; 
    atmp=0;
    for btmp=bst:db:bend
        ib = ib + 1;
        %% Setting up the system: X^(n+1) = A X^(n) + W^(n)
        % where X^(n) = [x^(n) y1^(n) y2^(n) .... yN^(n)]
        
        % Matrix A
        A = ones(Ny+1,Ny+1); 
        A=A*btmp; %%%%%% If y1,y2,y3... are all connected with each oher by "b"
        A(1,1)=atmp;
        for j=2:1:Ny+1
            A(1,j)=btmp;
        end
        for i=2:1:Ny+1
            A(i,1)= btmp + eps;
            A(i,i)= atmp + alph*eps;
        end
        
        % Spectral Radius should be <1
        rhoA = max(abs(eig(A))); 
        if (rhoA>=1) 
            bmax = btmp-db;
            break;
        end
    end
    ia=0; 
    btmp=0;
    for atmp=ast:da:aend
        ia = ia + 1;
        %% Setting up the system: X^(n+1) = A X^(n) + W^(n)
        % where X^(n) = [x^(n) y1^(n) y2^(n) .... yN^(n)]
        
        % Matrix A
        A = ones(Ny+1,Ny+1); 
        A=A*btmp; %%%%%% If y1,y2,y3... are all connected with each oher by "b"
        A(1,1)=atmp;
        for j=2:1:Ny+1
            A(1,j)=btmp;
        end
        for i=2:1:Ny+1
            A(i,1)= btmp + eps;
            A(i,i)= atmp + alph*eps;
        end
        
        % Spectral Radius should be <1
        rhoA = max(abs(eig(A)));
        if (rhoA>=1) 
            amax = atmp-da;
            break;
        end
    end

    % Check if our a,b is within this
    if (a > amax) && (b > bmax)
        disp('Error, a,b out of rho(A)<1 domain!');
        return;
    end

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
    if ((rhoA>=1) || (abs(det(IXI-AXA))<1e-10))
        netTExtoy(ialph) = nan; 
        netTExtoystar(ialph) = nan; 
        netTExtoystar2(ialph) = nan; 
%         discarded=discarded+1;
        break;
    end
    
    % Covariance matrix of noise W^(n)
    Q = zeros(Ny+1,Ny+1);
    Q(1,1)=sigw^2;
    for i=2:1:Ny+1
        Q(i,i)=(sigw + (beta*eps))^2;
    end
    
    %% Calculate Covariance matrices
    
    % C^(1)(X) ---------------------->
    % C_n,n
    vecQ = reshape(Q,[size(Q,1)*size(Q,2),1]);
    C11vec = (IXI-AXA)\vecQ; % IAinv = inv(IXI-AXA); C11vec = IAinv*vecQ;
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
    TEytox(ialph) = Hxx1 + Hxy - Hx - Hxyx1;
    
    % TE_x_yi
    Hyy1 = entropy([yi_idx,yi1_idx],C2X); % H(yi^(n),yi^(n+1))
    Hxy = entropy([x_idx,yi_idx],C2X); % H(x^(n),yi^(n))
    Hxyy1 = entropy([x_idx,yi_idx,yi1_idx],C2X);  % H(x^(n),yi^(n),yi^(n+1))
    Hy = entropy(yi_idx,C2X);  % H(yi^(n))
    TExtoy(ialph) = Hyy1 + Hxy - Hy - Hxyy1;
    
    netTExtoy(ialph) = TExtoy(ialph) - TEytox(ialph);
    
    % Individual dynamics
    Hx1y = entropy([x1_idx,yi_idx],C2X); % H(x^(n+1)|yi^(n))
    Hx1condy(ialph) = Hx1y-Hy;
    Hy1x = entropy([yi1_idx,x_idx],C2X); % H(yi^(n+1)|x^(n))
    Hy1condx(ialph) = Hy1x-Hx;

    Hx1condxy = Hxyx1 - Hxy;
    Hy1condxy = Hxyy1 - Hxy;
    indivx(ialph) = Hx1condy(ialph)-Hx1condxy; % H(x^(n+1)|yi^(n))-H(x^(n+1)|yi^(n),x^(n))
    indivy(ialph) = Hy1condx(ialph)-Hy1condxy;

    % Individual dynamics 2
    Hx1 = entropy(x1_idx,C2X);  % H(x^(n+1))
    indivx2(ialph) = Hx1 - Hxx1 + Hx;
    Hy1 = entropy(yi1_idx,C2X);  % H(y^(n+1))
    indivy2(ialph) = Hy1 - Hyy1 + Hy;

    % MI
    Hx1 = entropy(x1_idx,C2X);  % H(x^(n+1))
    Hy1 = entropy(yi1_idx,C2X);  % H(yi^(n+1))
    MIxtoy(ialph) = Hy1+Hx-Hy1x;
    MIytox(ialph) = Hx1+Hy-Hx1y;

    %% TE'_{x->y} = normalized TE
    Hcondx = Hxx1 - Hx;
    Hcondy = Hyy1 - Hy;
    normTExtoy(ialph) = TExtoy(ialph)/Hcondy;
    normTEytox(ialph) = TEytox(ialph)/Hcondx;
    netnormTExtoy(ialph) = TExtoy(ialph)/Hcondy - TEytox(ialph)/Hcondx;

    %% TE*_{x->y} = I(Xn+1;Yn+2|Yn+1,Xn)
    x_idx=1; x1_idx=(1+Ny)+1; x2_idx=2*(1+Ny)+1;
    i=1; yi_idx=1+i; yi1_idx=(1+Ny)+1+i; yi2_idx=2*(1+Ny)+1+i;

    % TE*_x_yi
    Hxy1y2 = entropy([x_idx,yi1_idx,yi2_idx],C3X); 
    Hxy1 = entropy([x_idx,yi1_idx],C3X); 
    Hxx1y1y2 = entropy([x_idx,x1_idx,yi1_idx,yi2_idx],C3X);  
    Hxx1y1 = entropy([x_idx,x1_idx,yi1_idx],C3X); 
    TExtoystar(ialph) = Hxy1y2 - Hxy1 - Hxx1y1y2 + Hxx1y1;
    % TE*_yi_x
    Hyx1x2 = entropy([yi_idx,x1_idx,x2_idx],C3X); 
    Hyx1 = entropy([yi_idx,x1_idx],C3X);
    Hyx1y1x2 = entropy([yi_idx,x1_idx,yi1_idx,x2_idx],C3X);
    Hyx1y1 = entropy([yi_idx,x1_idx,yi1_idx],C3X);
    TEytoxstar(ialph) = Hyx1x2 - Hyx1 - Hyx1y1x2 + Hyx1y1;

    netTExtoystar(ialph)= TExtoystar(ialph) - TEytoxstar(ialph);

    %% TE**_{x->y} = I(Xn+1;Yn+2|Yn+1,Xn,Yn)
    x_idx=1; x1_idx=(1+Ny)+1; x2_idx=2*(1+Ny)+1;
    i=1; yi_idx=1+i; yi1_idx=(1+Ny)+1+i; yi2_idx=2*(1+Ny)+1+i;

    % TE**_x_yi
    Hxyy1y2 = entropy([x_idx,yi_idx,yi1_idx,yi2_idx],C3X); 
    Hxyy1 = entropy([x_idx,yi_idx,yi1_idx],C3X); 
    Hxyx1y1y2 = entropy([x_idx,yi_idx,x1_idx,yi1_idx,yi2_idx],C3X);  
    Hxyx1y1 = entropy([x_idx,yi_idx,x1_idx,yi1_idx],C3X); 
    TExtoystar2(ialph) = Hxyy1y2 - Hxyy1 - Hxyx1y1y2 + Hxyx1y1;
    % TE**_yi_x
    Hxyx1x2 = entropy([x_idx,yi_idx,x1_idx,x2_idx],C3X); 
    Hxyx1 = entropy([x_idx,yi_idx,x1_idx],C3X);
    Hxyx1y1x2 = entropy([x_idx,yi_idx,x1_idx,yi1_idx,x2_idx],C3X);
    Hxyx1y1 = entropy([x_idx,yi_idx,x1_idx,yi1_idx],C3X);
    TEytoxstar2(ialph) = Hxyx1x2 - Hxyx1 - Hxyx1y1x2 + Hxyx1y1;

    netTExtoystar2(ialph)= TExtoystar2(ialph) - TEytoxstar2(ialph);

end

f1=figure;
plot(alpharr,TExtoy,'-k',LineWidth=2); hold on;
plot(alpharr,TEytox,'-.k',LineWidth=2); 
% xlim(xlimits); ylim([-0.005 0.015]);
xlabel('\alpha');ylabel('TE');
legend('TE_{X \rightarrow Y}','TE_{Y \rightarrow X}', ...
    'Location','north'); legend boxoff;
f1.Position = [100 100 260 240];

f2=figure;
plot(alpharr,netTExtoy,'-k',LineWidth=2);  hold on;
plot(alpharr,zeros(size(alpharr)),'--k',Linewidth=1.2);
% xlim(xlimits); 
xlabel('\alpha');ylabel('net TE_{X \rightarrow Y}');
f2.Position = [500 100 300 240];

f3=figure;
plot(alpharr,indivx,'-b',LineWidth=2); hold on;
plot(alpharr,indivy,'-.r',LineWidth=2); 
% xlim(xlimits); 
xlabel('\alpha');
ylabel('SE');
legend('SE_{X|Y}','SE_{Y|X}', ...
    'Location','northwest'); legend boxoff;
f3.Position = [900 100 300 240];

f4 = figure;
plot(alpharr,normTExtoy,'-','Color','#7E2F8E',LineWidth=2); hold on;
plot(alpharr,normTEytox,'-.','Color','#7E2F8E',LineWidth=2); 
xlabel('\alpha');ylabel("TE'");
% xlim(xlimits); 
legend("TE'_{X \rightarrow Y}","TE'_{Y \rightarrow X}", ...
    'Location','northwest'); legend boxoff;
f4.Position = [400 700 300 240];

f5 = figure;
plot(alpharr,TExtoystar,'-','Color','#D95319',LineWidth=2); hold on;
plot(alpharr,TEytoxstar,'-.','Color','#D95319',LineWidth=2); 
xlabel('\alpha');ylabel('TE^{C1}');
legend('TE^{C1}_{X \rightarrow Y}',...
    'TE^{C1}_{Y \rightarrow X}'); legend boxoff;
ylo = min(TEytoxstar(1),TExtoystar(1)); yup = max(TEytoxstar(1),TExtoystar(1)); ydel = yup - ylo;
dyt = 0.005; 
ylo = floor((ylo - ydel)/dyt)*dyt; yup = ceil((yup + ydel)/dyt)*dyt; 
% xlim(xlimits); ylim([0.1,0.3]);
f5.Position = [100 400 300 240];

f6=figure;
plot(alpharr,TExtoystar2,'-','Color','#4DBEEE',LineWidth=2); hold on;
plot(alpharr,TEytoxstar2,'-.','Color','#4DBEEE',LineWidth=2); 
xlabel('\alpha');ylabel('TE^{C2}');
legend('TE^{C2}_{X \rightarrow Y}',...
    'TE^{C2}_{Y \rightarrow X}'); legend boxoff; 
% xlim(xlimits); ylim([0.05,0.25]);
f6.Position = [500 400 300 240];

f7=figure;
plot(alpharr,netTExtoy,'-k',LineWidth=1.5); hold on;
plot(alpharr,netTExtoystar,'-.','Color','#D95319',LineWidth=2); hold on;
plot(alpharr,netTExtoystar2,'--','Color','#4DBEEE',LineWidth=2); hold on
plot(alpharr,zeros(size(alpharr)),'--k',Linewidth=1.2);
xlabel('\alpha');ylabel('net TE_{X \rightarrow Y}');
% xlim(xlimits); ylim([-0.02 0.02]);
legend('net TE_{X \rightarrow Y}', ...%     "net TE'_{X \rightarrow Y}",...
    'net TE^{C1}_{X \rightarrow Y}','net TE^{C2}_{X \rightarrow Y}', ...
    'Location','south'); legend boxoff;
f7.Position = [900 400 300 240];

return

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
    H = 0.5*log2( ((2*pi*exp(1))^nidx)*det(cov) );
end