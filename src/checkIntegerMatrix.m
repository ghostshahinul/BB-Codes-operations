function checkIntegerMatrix(A, name)
% checkIntegerMatrix  Reject non-integer numeric inputs.
    if any(abs(A(:) - round(A(:))) > 1e-10)
        error('%s has non-integer entries. This code is for Z, not floating-point fields.', name);
    end
end
