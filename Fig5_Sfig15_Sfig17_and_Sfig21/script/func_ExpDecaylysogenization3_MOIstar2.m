function fractionArray = func_ExpDecaylysogenization3_MOIstar2(a, q2, phi, k, Marray)
%(TN 2023/05/22) Calculating the natural logaritm of the fraction of lysogeny following the exponential-decay parametrization of the single-cell probability of lysogenization

%Preparations
MOIstar = 2;

fractionArray = zeros(size(Marray));
for i_M = 1:numel(Marray)
    M = Marray(i_M);
    nVectorAfter = MOIstar:1:max([10 4*M]);
    
    %Without the repression term, with varying q
    fractionBefore = poisspdf(1, a*M) .* phi * q2;

    %With the repression term, with q_MOI* only
    sumTerms = zeros(1, numel(nVectorAfter));
    for i_n = 1:numel(nVectorAfter)
        n = nVectorAfter(i_n);
        sumTerms(i_n) = poisspdf(n, a*M) .* q2 .* exp(-k .* (n - MOIstar));
    end
    fractionAfter = sum(sumTerms);

    fractionArray(i_M) = log(fractionBefore + fractionAfter);
end

