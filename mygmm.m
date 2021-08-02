function [gmm_mean, gmm_vcov, gmm_conf, gmm_hat, all_mu, spp] = mygmm(x,y,k,T)
disp('MCMC : Starts');
%% initiate variables
n = size(x,1); % Number of Observations
p = size(x,2); % Number of Dimensions In Data

z = zeros(n,1,T); % latent variable
q = zeros(k,1,T);
mu = zeros(k*p,1,T); % collection of mu vectors
R = zeros(k*p,p,T); % collection of R vectors

tau = zeros(k*p,1,2);
omega = zeros(k*p,p,2);

% variables : updated at each iteration
n_count = zeros(k,1);
x_bar = zeros(k*p,1);
S = zeros(k*p,p);
p_vec = zeros(n,k);
%% Initialization with k-means
disp(['Iteration: ',num2str(1)])
z(:,:,1) = kmeans(x,k)-1;
choto = 1e-3; %%%%% can be changed %%%%%
for i = 0:(k-1)
    gablu = (z(:,:,1)==i); % index
    hablu = ((i*p)+1):(p*(i+1)); % group index
    
    mu(hablu,:,1) = mean(x(gablu,:))';
    R(hablu,:,1) = inv(cov(x(gablu,:)));
    q(i+1,:,1) = sum(gablu)/n;
    
    tau(hablu,1,1) = zeros(p,1);
    omega(hablu,:,1) = choto*eye(p);
end
%% MCMC
for i = 2:T
    disp(['Iteration: ',num2str(i)])
    loop = i;
    for j = 0:(k-1)
        temp = (z(:,:,i-1)==j);
        n_count(j+1,1) = sum(temp);
        x_bar(((j*p)+1):(p*(j+1))) = mean(x(temp,:),1)';   
    end
    
    q(:,:,i) = drchrnd(ones(1,k) + n_count',1);
    
    for jj = 0:(k-1)
        
        temp1 = (z(:,:,i-1)==jj);
        
        temp2 = (1+(jj*p)):(p*(jj+1));
        
        omega(temp2,:,2)= n_count(jj+1)*R(temp2,:,i-1)...
            + omega(temp2,:,1);
        
        tau(temp2,:,2) = inv(omega(temp2,:,2))*(n_count(jj+1)*R(temp2,:,i-1)*x_bar(temp2)...
            + omega(temp2,:,1)*tau(temp2,:,1));
        
        mu(temp2,:,i) = mvnrnd(tau(temp2,:,2),inv(omega(temp2,:,2)),1);
        
        R(temp2,:,i) = wishrnd(inv(getS(x(temp1,:),mu(temp2,:,i))),n_count(jj+1));
    end
    %% create probability vector and z vector
    for ii = 1:n
        for u = 0:(k-1)
            temp3 = (1+(u*p)):(p*(u+1));
            p_vec(ii,u+1) = mvnpdf(x(ii,:)',mu(temp3,:,i),inv(R(temp3,:,i)));
        end
    end
    mamu(:,:,i) = q(:,:,i)'.*p_vec;
    pp(:,:,i) = mamu(:,:,i)./sum(mamu(:,:,i),2);
    z(:,:,i)=getz(mamu(:,:,i));
    n_count = zeros(k,1);
    x_bar = zeros(k*p,1);
    omega(:,:,1) = omega(:,:,2);
    omega(:,:,2) = zeros(k*p,p,1);
    tau(:,:,1) = tau(:,:,2);
    tau(:,:,2) = zeros(k*p,1,1);
end
%%
disp('MCMC: Final Calculation');
tt = floor(T/2);
gmm_mean = mean(mu(:,:,tt:end),3);
gmm_vcov = getinvR(mean(R(:,:,tt:end),3));
gmm_hat = mode(z(:,:,tt:end),3);
gmm_conf = confusionmat(y,gmm_hat);
spp = mean(pp(:,:,tt:end),3);
all_mu = mu;
disp('MCMC : Ends');
end