function checkCSSCommutationZ(HX, HZ, label)
% checkCSSCommutationZ  Verify HX*HZ.' = 0 over integers.
    C = double(HX) * double(HZ).';
    if any(C(:) ~= 0)
        error('%s do not commute over Z. nnz=%d, max=%g', label, nnz(C), max(abs(C(:))));
    end
end
