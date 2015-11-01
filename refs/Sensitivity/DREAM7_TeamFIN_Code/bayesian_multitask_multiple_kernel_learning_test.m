% Mehmet Gonen (mehmet.gonen@gmail.com)
% Helsinki Institute for Information Technology HIIT
% Department of Information and Computer Science
% Aalto University School of Science

function prediction = bayesian_multitask_multiple_kernel_learning_test(Km, state)
    T = length(Km);
    N = zeros(T, 1);
    for o = 1:T
        N(o) = size(Km{o}, 2);
    end
    P = size(Km{1}, 3);

    prediction.G = cell(1, T);
    for o = 1:T
        prediction.G{o}.mean = zeros(P, N(o));
        prediction.G{o}.covariance = zeros(P, N(o));
        for m = 1:P
            prediction.G{o}.mean(m, :) = state.a{o}.mean' * Km{o}(:, :, m);
            prediction.G{o}.covariance(m, :) = 1 / (state.upsilon.shape(o) * state.upsilon.scale(o)) + diag(Km{o}(:, :, m)' * state.a{o}.covariance * Km{o}(:, :, m));
        end
    end
    
    prediction.f = cell(1, T);
    for o = 1:T
        prediction.f{o}.mean = [ones(1, N(o)); prediction.G{o}.mean]' * state.be.mean([o, T + 1:T + P]);
        prediction.f{o}.covariance = 1 / (state.epsilon.shape(o) * state.epsilon.scale(o)) + diag([ones(1, N(o)); prediction.G{o}.mean]' * state.be.covariance([o, T + 1:T + P], [o, T + 1:T + P]) * [ones(1, N(o)); prediction.G{o}.mean]);
    end
end
