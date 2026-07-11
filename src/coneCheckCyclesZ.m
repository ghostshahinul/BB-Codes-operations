function coneCheckCyclesZ(D0, D1)
% coneCheckCyclesZ  Verify D0*D1 = 0 over Z.
    cyc = double(D0)*double(D1);
    if any(cyc(:) ~= 0)
        error('Cycle condition failed: D0*D1 is not zero over Z. max=%g', max(abs(cyc(:))));
    end
end
