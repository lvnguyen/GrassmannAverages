module GrassmannAverages

function trimmed_mean(arr, N, P)
    arr_sorted = sort(arr);
    low_index = N * P / 100;
    high_index = N - N * P/100;
    return mean(arr_sorted[low_index:high_index]);
end

function grassmann_median(X, K, P)
    K >= 0 || throw(ArgumentError("Number of components must be positive"));
    (P >= 0 && P < 50) || throw(ArgumentError("Trimming average percentage must be at least 0 and less than 50 percent"));
    
    N = size(X, 1);
    D = size(X, 2);
    eps = 1e-4;  "Too small of error is problematic"
    
    basis = Array{Real}(D, K);
    for k=1:K
        "q represents the k-th principal component"
        q = randn(D);
        q /= norm(q);
        
        while true
            dotsign = sign(X * q);
            
            "Compute the weighted average of the elements. Robustness is achieved by element-wise median."            
            newQ = Array{Real}(D);
            for j=1:D
                newQ[j] = trimmed_mean(X[: , j] .* dotsign, N, P);
            end
            newQ /= norm(newQ);
            
            "If the L_1 difference between q and newQ is small, then terminate."
            max_diff = 0;
            for j=1:D
                max_diff = max(max_diff, abs(q[j] - newQ[j]));
            end
            if max_diff < eps
                break
            end
            
            q = newQ;
        end
        
        "Store the basis vector and subtract from X"
        basis[: , k] = q;
        X -= (X * q) * q';
    end
    
    return basis;
end

end # module
