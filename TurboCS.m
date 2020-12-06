%**************************************************************************
% This code reproduces the results in the following paper

% "Turbo compressed sensing with partial DFT sensing matrix",
% available at http://arxiv.org/abs/1408.3904

% If you have used this code, please cite the corresponding paper. Thanks.

% Junjie Ma, 26 Augest, 2014.
%**************************************************************************

clear all;

lambda = 0.4;                       % sparsity level

N = 32768;                           % number of variables.
% It is set to a very large number to see the agreement with evolution. You can try a smaller one.

M = N * 0.55;                        % number of measurements
M = fix(M);

NSIM = 10;                        % number of simulations. To save time, you can try a much smaller number, say 10.

SNRdB = [50:4:50];

Iteration = 70;                     % number of iterations

MSE_simulation = zeros(Iteration,length(SNRdB));
MSE_evolution  = zeros(Iteration,length(SNRdB));
MSE_GAMP_cal=zeros(Iteration,length(SNRdB));
MSE_GAMP_mea=zeros(Iteration,length(SNRdB));
EXP_MAX = 50;
EXP_MIN = -50;

LOG_MAX = 1e50;
LOG_MIN = 1e-50;

load Table_sparsity04;              % For evolution, the function mmse(eta) (see Section III-A) is saved as a table. Only for lambda = 0.4
SNR_Table = Table_sparsity04(:,1);
VAR_POST_Table = Table_sparsity04(:,2);

for iii = 1:length(SNRdB)
    
    sigma2 = 10^(-SNRdB(iii)/10);
    
    rand('state',0);
    randn('state',0);
    
    %%  Evolution, only for the case sparsity = 0.4
    
    v_A_pri_evo = 1;
    
    for it = 1:Iteration
        
        v_A2B_ext_evo = N/M * (v_A_pri_evo + sigma2) - v_A_pri_evo;
        
        v_B_pri_evo = v_A2B_ext_evo;
        
        v_B_post_evo = interp1q(SNR_Table,VAR_POST_Table,1/v_B_pri_evo);
        
        v_B_ext_evo = 1/(1/v_B_post_evo-1/v_B_pri_evo);
        
        v_A_pri_evo = v_B_ext_evo;
        
        MSE_evolution(it, iii) = v_B_post_evo;
    end
    
    %% Simulation
    
    for nsim = 1:NSIM
        
        pattern = random('Binomial',1,lambda,N,1);
        x = sqrt(1/lambda) * (randn(N,1)+1i*randn(N,1))/sqrt(2);
        x = x .* pattern;
        
        permuation = randperm(N);
        index_selection = permuation(1:M);
        
        z = fft(x)/sqrt(N);
        y = z(index_selection) + sqrt(sigma2)*(randn(M,1)+1i*randn(M,1))/sqrt(2);
        
        % Signal recovery
        
        v_A_pri = 1;
        z_A_pri = zeros(N,1);
        
        % The implementations are slightly different from that shown in Fig. 1, see the discussions in Section II-B-(4)
        for it = 1:Iteration
            
            % compute a posteriori
            z_A_post = zeros(N,1);
            z_A_post(index_selection) = v_A_pri / ( v_A_pri + sigma2 ) * ( y - z_A_pri(index_selection) );
            z_A_post = z_A_post + z_A_pri;
            v_A_post = v_A_pri - M/N * v_A_pri^2 / ( v_A_pri + sigma2 );
            
            % update extrinsic
            %             v_A2B_ext = N/M * ( v_A_pri + sigma2 ) - v_A_pri;             % equal to 1/(1/v_A_post- 1/v_A_pri)
            v_A2B_ext = 1/( 1 / v_A_post- 1 / v_A_pri );
            z_A2B_ext = ( v_A2B_ext / v_A_post ) * z_A_post - ( v_A2B_ext / v_A_pri) * z_A_pri;
            
            % tranform to the x domain
            x_B_pri = ifft(z_A2B_ext) * sqrt(N);
            v_B_pri = v_A2B_ext;
            qqplot(real(x_B_pri-x))
            % MMSE sparsity combiner (also called Baye optimal denoiser in the literature)
            
            exponent =  - abs(x_B_pri).^2 * 1/lambda / ( v_B_pri * (v_B_pri+1/lambda));
            
            exponent(exponent > EXP_MAX) = EXP_MAX;
            exponent(exponent < EXP_MIN) = EXP_MIN;
            
            den = 1 + ( v_B_pri + 1/lambda ) / v_B_pri * (1-lambda) / lambda * exp(exponent);
            C = 1./den;
            
            % compute a posteriori
            
            VAR = 1/lambda * v_B_pri / ( v_B_pri + 1/lambda ) * C + abs( 1/lambda / ( 1/lambda + v_B_pri ) * x_B_pri ).^2 .* C .* (1-C);
            v_B_post = mean(VAR);
            
            x_B_post = ( 1/lambda ) / ( 1/lambda + v_B_pri ) * x_B_pri .* C;
            qqplot(real(x_B_post-x))
            skdl=norm(x_B_post-x)^2/N;
            % update extrinsic
            
            v_B2A_ext = 1/(1/v_B_post-1/v_B_pri);
            x_B2A_ext = ( v_B2A_ext / v_B_post ) * x_B_post - (v_B2A_ext / v_B_pri) * x_B_pri;
            testx=norm(x_B2A_ext-x)^2/N;
            % transorm to the z domain
            
            z_A_pri = fft(x_B2A_ext) / sqrt(N);
            v_A_pri = v_B2A_ext;
            
            % measure MSE
            
            MSE_simulation(it,iii) = MSE_simulation(it,iii) + norm(x_B_post - x)^2/(N);
            
            
        end
        
        u_hat = zeros(M,1);
        DIFF = 0;
        x_hat = zeros(N,1);
        for it = 1:Iteration
            Hxtmp=fft(x_hat)/sqrt(N);
            Hx=Hxtmp(index_selection);
            u_hat =N/M*( y - Hx  + u_hat * DIFF);
            
            v_u = 1/N * norm(u_hat)^2;
            
            Htu=zeros(N,1);
            Htu(index_selection)=u_hat;
            Htu=ifft(Htu)*sqrt(N);
            z_PRI = x_hat + Htu;
            %qqplot(real(z_PRI-x))
            Var_PRI = v_u;
            
            %------------------ MMSE with Gauss-Bernoulli prior------------
            % see note "MMSE with Gauss-Bernoulli sparse prior"
            % Define p0 := Pr(b=0|y) and p1 := Pr(b=1|y)
            exponent =  - abs(z_PRI).^2 * 1/lambda ./ ( Var_PRI .* (Var_PRI+1/lambda));
            exponent(exponent > EXP_MAX) = EXP_MAX;
            exponent(exponent < EXP_MIN) = EXP_MIN;
            den = 1 + (Var_PRI+1/lambda)./Var_PRI * (1-lambda)/lambda .* exp(exponent);
            C = 1./den;
            
            z_POST = ( 1/lambda ) ./ ( 1/lambda + Var_PRI ) .* z_PRI .* C;
            
            VAR = 1/lambda.*Var_PRI./(Var_PRI+1/lambda) .* C + abs(1/lambda./(1./lambda+Var_PRI) .* z_PRI).^2 .* C .* (1-C);
            
            
            %------------------------------
            x_POST = z_POST;
            
            Var_POST = mean(VAR);
            
            x_hat = x_POST;
            DIFF = Var_POST / Var_PRI;
            
            
            MSE_GAMP_mea(it,iii) =  MSE_GAMP_mea(it,iii)+norm(x_POST - x)^2/(N);
            MSE_GAMP_cal(it,iii) = MSE_GAMP_cal(it,iii)+mean(VAR);
            
        end
        
        if(mod(nsim,10)==0)
            fprintf('======================================================\n');
            fprintf('snr = %d dB, nsim = %d\n',SNRdB(iii),nsim);
            fprintf('MSE_sim = %e, MSE_evo = %e\n',MSE_simulation(Iteration,iii)/nsim,MSE_evolution(Iteration,iii));
        end
        
    end
    
    MSE_simulation(:,iii) = MSE_simulation(:,iii)/NSIM;
    MSE_GAMP_cal(:,iii) = MSE_GAMP_cal(:,iii)/NSIM;
    MSE_GAMP_mea(:,iii) = MSE_GAMP_mea(:,iii)/NSIM;
    
end

% save;
semilogy(1:Iteration,MSE_simulation,'r--',1:Iteration,MSE_evolution,'bp',1:Iteration,MSE_GAMP_cal);
legend('Turbo-CS, partial DFT, simulation','Turbo-CS, partial DFT, evolution','AMP, partial DFT, simulation');
xlabel('Iteration');
ylabel('MSE');
for it = 1:Iteration
    semilogy(SNRdB,MSE_simulation(it,:),'r--',SNRdB,MSE_evolution(it,:),'bp');
    hold on;
end
legend('simulation','evolution');
xlabel('SNR');
ylabel('MSE');
saveas(gcf,'MSE.fig');
saveas(gcf,'MSE.png');
close;


figure;
for iii = 1:length(SNRdB)
    semilogy(1:Iteration,MSE_simulation(:,iii),'r--',1:Iteration,MSE_evolution(:,iii),'bp');
    hold on;
end
legend('simulation','evolution');
xlabel('it');
ylabel('MSE');
saveas(gcf,'MSE2.fig');
saveas(gcf,'MSE2.png');
close;


