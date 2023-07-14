categorical.baseline <- "
functions{
    matrix cov_GPL2(matrix x, real a, real b, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;

        for (i in 1:(N-1)) {
          K[i, i] = a + delta;
          for (j in (i + 1):N) {
            K[i, j] = a * exp(- b * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = a + delta;
        return K;
    }
    real partial_sum(int[] y_slice,
                   int start, int end,
                   int N, int M, vector alpha, vector beta, vector gamma, vector phi, int[ , ] Y, row_vector X1, real minp ) {
        real lp = 0.0;
        real p;
        vector[M+1] prob;

        for( i in start:end){
            for( j in 1:N ){
            
                p = exp(- alpha[i] - 0.5 * gamma[i] * pow(X1[j] - beta[i], 2));
            
                prob[1] = 1 - phi[1]*p - minp*(M+1);
                       
                for (k in 2:M){
                   prob[k]  =  phi[k-1]*p - phi[k]*p + minp;
                }
                       
                prob[M+1]  = phi[M]*p + minp;

                lp += categorical_lpmf( Y[i, j] | prob);
            }
        }
        return lp;
    }
}
data{
    int N;
    int L;
    int M;
    real minp;
    int Y[L, N];
    int indices[L];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    real<lower=0, upper=1> zphi[M];
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_b;
    real<lower=0> etasq_b;
    real<lower=0> rhosq_b;
    real<lower=0> sigma_g;
    real<lower=0> etasq_g;
    real<lower=0> rhosq_g;
}
transformed parameters{
    vector[L] alpha;
    vector[L] beta;
    vector[L] gamma;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;

    vector[M] phi;
    phi[1] = zphi[1];
    for (i in 2:M){
       phi[i] = zphi[i] * phi[i - 1];
    }

    alpha = exp(zalpha * sigma_a + alpha_bar);

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
}
model{
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    rhosq_g ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    zgamma ~ std_normal();
    zbeta ~ std_normal();
    zphi ~ beta(1, 1);
    
    int grainsize = 1;

    target += reduce_sum(partial_sum, indices,
                     grainsize,
                     N, M, alpha, beta, gamma, phi, Y, X1, minp);
}
generated quantities{
    vector[L*N] log_lik;
    int k;
    real p;
    vector[M+1] prob;

    k = 1;
    for ( i in 1:L ){
        for (j in 1:N){
           p = exp(-alpha[i] - 0.5 * gamma[i] * pow(X1[j] - beta[i], 2));

           prob[1] = 1 - phi[1]*p - minp*(M+1);

           for (l in 2:M){
              prob[l]  =  phi[l-1]*p - phi[l]*p + minp;
           }

           prob[M+1]  = phi[M]*p + minp;

           log_lik[k] = categorical_lpmf( Y[i, j] | prob);
           k = k + 1;
        }
    }
}
"

categorical.fat.tailed.skewed <- "
functions{
    matrix cov_GPL2(matrix x, real a, real b, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;

        for (i in 1:(N-1)) {
          K[i, i] = a + delta;
          for (j in (i + 1):N) {
            K[i, j] = a * exp(- b * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = a + delta;
        return K;
    }
    int sgn(real x) {
        return x < 0 ? -1 : x > 0;
    }
    real partial_sum(int[] y_slice,
                   int start, int end,
                   int N, int M, vector lambda, vector alpha, vector beta, vector gamma, vector nu, vector phi, int[ , ] Y, row_vector X1, real minp ) {
        real lp = 0.0;
        real m;
        real v;
        real p;
        vector[M+1] prob;

        for( i in start:end){
            v = sqrt(gamma[i] * ((pi()*(1+3*pow(lambda[i],2))*tgamma(3.0/nu[i])-pow(16,(1.0/nu[i]))*pow(lambda[i],2)*tgamma(0.5+1.0/nu[i])*tgamma(0.5+1.0/nu[i])*tgamma(1.0/nu[i]))/(pi()*tgamma(1.0/nu[i]))));
            m = beta[i] - pow(2, 2.0/nu[i])*lambda[i]*tgamma(0.5+1.0/nu[i])/(sqrt(pi())*v);
            for( j in 1:N ){
            
                p = exp(-alpha[i] - pow(v * fabs(X1[j] - m)/(1+lambda[i]*sgn(X1[j] - m)), nu[i]));
            
                prob[1] = 1 - phi[1]*p - minp*(M+1);
                       
                for (k in 2:M){
                   prob[k]  =  phi[k-1]*p - phi[k]*p + minp;
                }
                       
                prob[M+1]  = phi[M]*p + minp;

                lp += categorical_lpmf( Y[i, j] | prob);
            }
        }
        return lp;
    }
}
data{
    int N;
    int L;
    int M;
    real minp;
    int Y[L, N];
    int indices[L];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    real<lower=0, upper=1> zphi[M];
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    vector[L] znu;
    vector[L] zlambda;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real nu_bar;
    real lambda_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_n;
    real<lower=0> sigma_l;
    real<lower=0> sigma_b;
    real<lower=0> etasq_b;
    real<lower=0> rhosq_b;
    real<lower=0> sigma_g;
    real<lower=0> etasq_g;
    real<lower=0> rhosq_g;
}
transformed parameters{
    vector[L] lambda;
    vector[L] alpha;
    vector[L] beta;
    vector[L] gamma;
    vector[L] nu;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;
    
    vector[M] phi;
    phi[1] = zphi[1];
    for (i in 2:M){
       phi[i] = zphi[i] * phi[i - 1];
    }

    alpha = exp(zalpha * sigma_a + alpha_bar);
    nu = exp(znu * sigma_n + nu_bar)+1;
    lambda = inv_logit(zlambda*sigma_l + lambda_bar)*2-1;

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
}
model{
    sigma_a ~ exponential( 1 );
    sigma_n ~ exponential( 1 );
    sigma_l ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    rhosq_g ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    lambda_bar ~ std_normal();
    nu_bar ~ std_normal();
    zalpha ~ std_normal();
    zgamma ~ std_normal();
    zbeta ~ std_normal();
    znu ~ std_normal();
    zlambda ~ std_normal();
    zphi ~ beta(1, 1);
    
    int grainsize = 1;

    target += reduce_sum(partial_sum, indices,
                     grainsize,
                     N, M, lambda, alpha, beta, gamma, nu, phi, Y, X1, minp);
}
generated quantities{
    vector[L*N] log_lik;
    vector[L] betam;
    vector[L] gammav;
    int k;
    real p;
    vector[M+1] prob;

    for (i in 1:L){
        gammav[i] = sqrt(gamma[i] * ((pi()*(1+3*pow(lambda[i],2))*tgamma(3.0/nu[i])-pow(16,(1.0/nu[i]))*pow(lambda[i],2)*tgamma(0.5+1.0/nu[i])*tgamma(0.5+1.0/nu[i])*tgamma(1.0/nu[i]))/(pi()*tgamma(1.0/nu[i]))));
        betam[i] = beta[i] - pow(2, 2.0/nu[i])*lambda[i]*tgamma(0.5+1.0/nu[i])/(sqrt(pi())*gammav[i]);
    }

    k = 1;
    for ( i in 1:L ){
        for (j in 1:N){
           p = exp(-alpha[i] - pow(gammav[i] * fabs(X1[j] - betam[i])/(1+lambda[i]*sgn(X1[j] - betam[i])), nu[i]));

           prob[1] = 1 - phi[1]*p - minp*(M+1);

           for (l in 2:M){
              prob[l]  =  phi[l-1]*p - phi[l]*p + minp;
           }

           prob[M+1]  = phi[M]*p + minp;

           log_lik[k] = categorical_lpmf( Y[i, j] | prob);
           k = k + 1;
        }
    }
}
"

