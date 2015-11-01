% Mehmet Gonen (mehmet.gonen@gmail.com)
% Helsinki Institute for Information Technology HIIT
% Department of Information and Computer Science
% Aalto University School of Science

function state = bayesian_multitask_multiple_kernel_learning_train(Km, y, parameters)
    rand('state', parameters.seed); %#ok<RAND>
    randn('state', parameters.seed); %#ok<RAND>

    T = length(Km);
    D = zeros(T, 1);
    N = zeros(T, 1);
    for o = 1:T
        D(o) = size(Km{o}, 1);
        N(o) = size(Km{o}, 2);
    end
    P = size(Km{1}, 3);
    
    log2pi = log(2 * pi);

    lambda = cell(1, T);
    for o = 1:T
        lambda{o}.shape = (parameters.alpha_lambda + 0.5) * ones(D(o), 1);
        lambda{o}.scale = parameters.beta_lambda * ones(D(o), 1);
    end
    upsilon.shape = (parameters.alpha_upsilon + 0.5 * N * P) .* ones(T, 1);
    upsilon.scale = parameters.beta_upsilon * ones(T, 1);
    a = cell(1, T);    
    for o = 1:T
        a{o}.mean = randn(D(o), 1);
        a{o}.covariance = eye(D(o), D(o));
    end
    G = cell(1, T);
    for o = 1:T
        G{o}.mean = randn(P, N(o));
        G{o}.covariance = eye(P, P);
    end
    gamma.shape = (parameters.alpha_gamma + 0.5) * ones(T, 1);
    gamma.scale = parameters.beta_gamma * ones(T, 1);    
    omega.shape = (parameters.alpha_omega + 0.5) * ones(P, 1);
    omega.scale = parameters.beta_omega * ones(P, 1);
    epsilon.shape = (parameters.alpha_epsilon + 0.5 * N) .* ones(T, 1);
    epsilon.scale = parameters.beta_epsilon * ones(T, 1);
    be.mean = [zeros(T, 1); ones(P, 1)];
    be.covariance = eye(T + P, T + P);

    KmKm = cell(1, T);
    for o = 1:T
        KmKm{o} = zeros(D(o), D(o));
        for m = 1:P
            KmKm{o} = KmKm{o} + Km{o}(:, :, m) * Km{o}(:, :, m)';
        end
        Km{o} = reshape(Km{o}, [D(o), N(o) * P]);
    end
    
    if parameters.progress == 1
        bounds = zeros(parameters.iteration, 1);
    end

    atimesaT = cell(1, T);
    for o = 1:T
        atimesaT{o}.mean = a{o}.mean * a{o}.mean' + a{o}.covariance;
    end
    GtimesGT = cell(1, T);
    for o = 1:T
        GtimesGT{o}.mean = G{o}.mean * G{o}.mean' + N(o) * G{o}.covariance;
    end
    btimesbT.mean = be.mean(1:T) * be.mean(1:T)' + be.covariance(1:T, 1:T);
    etimeseT.mean = be.mean(T + 1:T + P) * be.mean(T + 1:T + P)' + be.covariance(T + 1:T + P, T + 1:T + P);
    etimesb.mean = zeros(P, T);
    for o = 1:T
        etimesb.mean(:, o) = be.mean(T + 1:T + P) * be.mean(o) + be.covariance(T + 1:T + P, o);
    end
    KmtimesGT = cell(1, T);
    for o = 1:T
        KmtimesGT{o}.mean = Km{o} * reshape(G{o}.mean', N(o) * P, 1);
    end
    for iter = 1:parameters.iteration
        if mod(iter, 1) == 0
            fprintf(1, '.');
        end
        if mod(iter, 10) == 0
            fprintf(1, ' %5d\n', iter);
        end

        %%%% update lambda
        for o = 1:T
            lambda{o}.scale = 1 ./ (1 / parameters.beta_lambda + 0.5 * diag(atimesaT{o}.mean));
        end
        %%%% update upsilon
        for o = 1:T
            upsilon.scale(o) = 1 / (1 / parameters.beta_upsilon + 0.5 * (sum(diag(GtimesGT{o}.mean)) ...
                                                                         - 2 * sum(sum(reshape(a{o}.mean' * Km{o}, [N(o), P])' .* G{o}.mean)) ...
                                                                         + sum(diag(KmKm{o} * atimesaT{o}.mean))));
        end
        %%%% update a
        for o = 1:T
            a{o}.covariance = (diag(lambda{o}.shape .* lambda{o}.scale) + upsilon.shape(o) * upsilon.scale(o) * KmKm{o}) \ eye(D(o), D(o));
            a{o}.mean = a{o}.covariance * (upsilon.shape(o) * upsilon.scale(o) * KmtimesGT{o}.mean);
            atimesaT{o}.mean = a{o}.mean * a{o}.mean' + a{o}.covariance;
        end
        %%%% update G        
        for o = 1:T
            G{o}.covariance = (upsilon.shape(o) * upsilon.scale(o) * eye(P, P) + epsilon.shape(o) * epsilon.scale(o) * etimeseT.mean) \ eye(P, P);
            G{o}.mean = G{o}.covariance * (upsilon.shape(o) * upsilon.scale(o) * reshape(a{o}.mean' * Km{o}, [N(o), P])' + epsilon.shape(o) * epsilon.scale(o) * (be.mean(T + 1:T + P) * y{o}' - repmat(etimesb.mean(:, o), 1, N(o))));
            GtimesGT{o}.mean = G{o}.mean * G{o}.mean' + N(o) * G{o}.covariance;
            KmtimesGT{o}.mean = Km{o} * reshape(G{o}.mean', N(o) * P, 1);
        end   
        %%%% update gamma
        gamma.scale = 1 ./ (1 / parameters.beta_gamma + 0.5 * diag(btimesbT.mean));
        %%%% update omega
        omega.scale = 1 ./ (1 / parameters.beta_omega + 0.5 * diag(etimeseT.mean));
        %%%% update epsilon
        for o = 1:T
            epsilon.scale(o) = 1 / (1 / parameters.beta_epsilon + 0.5 * (y{o}' * y{o} - 2 * y{o}' * [ones(1, N(o)); G{o}.mean]' * be.mean([o, T + 1:T + P]) ...
                                                                         + N(o) * btimesbT.mean(o, o) ...
                                                                         + sum(diag(GtimesGT{o}.mean * etimeseT.mean)) ...
                                                                         + 2 * sum(G{o}.mean, 2)' * etimesb.mean(:, o)));
        end
        %%%% update b and e
        be.covariance = [diag(gamma.shape .* gamma.scale) + diag(N .* epsilon.shape .* epsilon.scale), zeros(T, P); ...
                         zeros(P, T), diag(omega.shape .* omega.scale)];
        for o = 1:T
            be.covariance(T + 1:T + P, o) = epsilon.shape(o) * epsilon.scale(o) * sum(G{o}.mean, 2);
            be.covariance(o, T + 1:T + P) = epsilon.shape(o) * epsilon.scale(o) * sum(G{o}.mean, 2)';
            be.covariance(T + 1:T + P, T + 1:T + P) = be.covariance(T + 1:T + P, T + 1:T + P) + epsilon.shape(o) * epsilon.scale(o) * GtimesGT{o}.mean;
        end
        be.covariance = be.covariance \ eye(T + P, T + P);
        be.mean = zeros(T + P, 1);        
        for o = 1:T
            be.mean(o) = epsilon.shape(o) * epsilon.scale(o) * sum(y{o});
            be.mean(T + 1:T + P) = be.mean(T + 1:T + P) + epsilon.shape(o) * epsilon.scale(o) * G{o}.mean * y{o};
        end
        be.mean = be.covariance * be.mean;
        btimesbT.mean = be.mean(1:T) * be.mean(1:T)' + be.covariance(1:T, 1:T);
        etimeseT.mean = be.mean(T + 1:T + P) * be.mean(T + 1:T + P)' + be.covariance(T + 1:T + P, T + 1:T + P);
        for o = 1:T
            etimesb.mean(:, o) = be.mean(T + 1:T + P) * be.mean(o) + be.covariance(T + 1:T + P, o);
        end
        
        if parameters.progress == 1
            lb = 0;

            %%%% p(lambda)
            for o = 1:T
                lb = lb + sum((parameters.alpha_lambda - 1) * (psi(lambda{o}.shape) + log(lambda{o}.scale)) ...
                              - lambda{o}.shape .* lambda{o}.scale / parameters.beta_lambda ...
                              - gammaln(parameters.alpha_lambda) ...
                              - parameters.alpha_lambda * log(parameters.beta_lambda));
            end
            %%%% p(upsilon)
            lb = lb + sum((parameters.alpha_upsilon - 1) * (psi(upsilon.shape) + log(upsilon.scale)) ...
                          - upsilon.shape .* upsilon.scale / parameters.beta_upsilon ...
                          - gammaln(parameters.alpha_upsilon) ...
                          - parameters.alpha_upsilon * log(parameters.beta_upsilon));
            %%%% p(a | lambda)
            for o = 1:T
                lb = lb - 0.5 * sum(diag(diag(lambda{o}.shape .* lambda{o}.scale) * atimesaT{o}.mean)) ...
                        - 0.5 * (D(o) * log2pi - sum(log(lambda{o}.shape .* lambda{o}.scale)));
            end
            %%%% p(G | a, Km, upsilon)
            for o = 1:T
                lb = lb - 0.5 * sum(diag(GtimesGT{o}.mean)) * upsilon.shape(o) * upsilon.scale(o) ...
                        + (a{o}.mean' * KmtimesGT{o}.mean) * upsilon.shape(o) * upsilon.scale(o) ...
                        - 0.5 * sum(diag(KmKm{o} * atimesaT{o}.mean)) * upsilon.shape(o) * upsilon.scale(o) ...
                        - 0.5 * N(o) * P * (log2pi - log(upsilon.shape(o) * upsilon.scale(o)));
            end
            %%%% p(gamma)
            lb = lb + sum((parameters.alpha_gamma - 1) * (psi(gamma.shape) + log(gamma.scale)) ...
                          - gamma.shape .* gamma.scale / parameters.beta_gamma ...
                          - gammaln(parameters.alpha_gamma) ...
                          - parameters.alpha_gamma * log(parameters.beta_gamma));
            %%%% p(b | gamma)
            lb = lb - 0.5 * sum(diag(diag(gamma.shape .* gamma.scale) * btimesbT.mean)) ...
                    - 0.5 * (T * log2pi - sum(log(gamma.shape .* gamma.scale)));
            %%%% p(omega)
            lb = lb + sum((parameters.alpha_omega - 1) * (psi(omega.shape) + log(omega.scale)) ...
                          - omega.shape .* omega.scale / parameters.beta_omega ...
                          - gammaln(parameters.alpha_omega) ...
                          - parameters.alpha_omega * log(parameters.beta_omega));
            %%%% p(e | omega)
            lb = lb - 0.5 * sum(diag(diag(omega.shape .* omega.scale) * etimeseT.mean)) ...
                    - 0.5 * (P * log2pi - sum(log(omega.shape .* omega.scale)));
            %%%% p(epsilon)
            lb = lb + sum((parameters.alpha_epsilon - 1) * (psi(epsilon.shape) + log(epsilon.scale)) ...
                          - epsilon.shape .* epsilon.scale / parameters.beta_epsilon ...
                          - gammaln(parameters.alpha_epsilon) ...
                          - parameters.alpha_epsilon * log(parameters.beta_epsilon));
            %%%% p(y | b, e, G, epsilon)
            for o = 1:T
                lb = lb - 0.5 * (y{o}' * y{o}) * epsilon.shape(o) * epsilon.scale(o) ...
                        + (y{o}' * (G{o}.mean' * be.mean(T + 1:T + P))) * epsilon.shape(o) * epsilon.scale(o) ...
                        + sum(be.mean(o) * y{o}) * epsilon.shape(o) * epsilon.scale(o) ...
                        - 0.5 * sum(diag(etimeseT.mean * GtimesGT{o}.mean)) * epsilon.shape(o) * epsilon.scale(o) ...
                        - sum(G{o}.mean' * etimesb.mean(:, o)) * epsilon.shape(o) * epsilon.scale(o) ...
                        - 0.5 * N(o) * btimesbT.mean(o, o) * epsilon.shape(o) * epsilon.scale(o) ...
                        - 0.5 * N(o) * (log2pi - log(epsilon.shape(o) * epsilon.scale(o)));
            end

            %%%% q(lambda)
            for o = 1:T
                lb = lb + sum(lambda{o}.shape + log(lambda{o}.scale) + gammaln(lambda{o}.shape) + (1 - lambda{o}.shape) .* psi(lambda{o}.shape));
            end
            %%%% q(upsilon)
            lb = lb + sum(upsilon.shape + log(upsilon.scale) + gammaln(upsilon.shape) + (1 - upsilon.shape) .* psi(upsilon.shape));            
            %%%% q(a)
            for o = 1:T
                lb = lb + 0.5 * (D(o) * (log2pi + 1) + logdet(a{o}.covariance));
            end
            %%%% q(G)
            for o = 1:T
                lb = lb + 0.5 * N(o) * (P * (log2pi + 1) + logdet(G{o}.covariance));
            end
            %%%% q(gamma)
            lb = lb + sum(gamma.shape + log(gamma.scale) + gammaln(gamma.shape) + (1 - gamma.shape) .* psi(gamma.shape));
            %%%% q(omega)
            lb = lb + sum(omega.shape + log(omega.scale) + gammaln(omega.shape) + (1 - omega.shape) .* psi(omega.shape));
            %%%% q(epsilon)
            lb = lb + sum(epsilon.shape + log(epsilon.scale) + gammaln(epsilon.shape) + (1 - epsilon.shape) .* psi(epsilon.shape));
            %%%% q(b, e)
            lb = lb + 0.5 * ((T + P) * (log2pi + 1) + logdet(be.covariance));

            bounds(iter) = lb;
        end
    end

    state.lambda = lambda;
    state.upsilon = upsilon;
    state.a = a;
    state.gamma = gamma;
    state.omega = omega;
    state.epsilon = epsilon;
    state.be = be;
    if parameters.progress == 1
        state.bounds = bounds;
    end
    state.parameters = parameters;
end

function ld = logdet(Sigma)
    U = chol(Sigma);
	ld = 2 * sum(log(diag(U)));
end
