function frequencyArray = func_gParametrizedLysogeny_noshift(params, inputs)
%(TN 2023/06/01) Calculating the natural logaritm of the fraction of lysogeny
% params = 1-by-10 vector, for the g-dependent parametrization
% inputs = n-by-2 matrix, each row being [g M], where g is the growth rate and M is the bulk MOI
% This is a vectorized function, returning a column vector.

%Arguments validation
arguments
    params (1, 11) double
    inputs (:, 2) double
end

%Preparations of parameters
MOIstar = 2;
quadratic.function = @(a, b, c, x) a.*x.^2 + b.*x + c;
piecewiseQuadratic.function = @(g1, g2, x) (0) .* (x < g1)...
    - (x - g1) .* (x - g2) .* (x >= g1 & x < g2)...
    + (0) .* (x >= g2);
linear.function = @(a, b, x) a.*x + b;
piecewiseTwoSlopes.function = @(a, b, c, gstar, x) (a.*x + b) .* (x < gstar) + (c.*x + (a-c).*gstar + b) .* (x >= gstar);
piecewiseSlopeFirst.function = @(a, b, gstar, x) (a.*x + b) .* (x < gstar) + (a.*gstar + b) .* (x >= gstar);
piecewiseFlatFirst.function = @(a, b, gstar, x) (a.*gstar + b) .* (x < gstar) + (a.*x + b) .* (x >= gstar);

frequencyArray = zeros(size(inputs, 1), 1);
for i_value = 1:size(inputs, 1)
    g = inputs(i_value, 1);
    M = inputs(i_value, 2);

    %Specification of g-dependent parameters
    a = 1;%exp(quadratic.function(params(1), params(2), params(3), g));
    k = piecewiseQuadratic.function(params(4), params(5), g);
    q1 = exp(linear.function(params(6), params(7), g));
    q2 = exp(piecewiseTwoSlopes.function(params(8), params(9), params(10), params(11), g));
%     phiRecip = exp(piecewiseSlopeFirst.function(params(8), params(9), params(10), g));
%     q2 = q1 * phiRecip;
%     q1 = q2 / phiRecip;

    %Calculation of the frequency of lysogeny for this M value
    nVectorAfter = MOIstar:1:max([10 4*M]);

    %Without the repression term, with varying q
    fractionBefore = poisspdf(1, a*M) .* q1;

    %With the repression term, with q_MOI* only
    sumTerms = zeros(1, numel(nVectorAfter));
    for i_n = 1:numel(nVectorAfter)
        n = nVectorAfter(i_n);
        sumTerms(i_n) = poisspdf(n, a*M) .* q2 .* exp(-k .* (n - MOIstar));
    end
    fractionAfter = sum(sumTerms);

    frequencyArray(i_value) = log(fractionBefore + fractionAfter);
end