function coneCheckCoreZ(HX, F1, F0, D1)
% coneCheckCoreZ  Verify F0*D1 + HX*F1 = 0 over Z.
    core = double(F0)*double(D1) + double(HX)*double(F1);
    if any(core(:) ~= 0)
        error('Core condition failed: F0*D1 + HX*F1 is not zero over Z. max=%g', max(abs(core(:))));
    end
end
