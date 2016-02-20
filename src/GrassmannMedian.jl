module GrassmannAverages

function grassmann_median(X, K)
    K >= 0 || throw(ArgumentError("Number of components must be positive"));
    
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
                newQ[j] = median(X[: , j] .* dotsign);
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
