% Verify the condition-aware BB [[72,12,6]] signed-integer result.
load('BB72_conditional_final_matrices.mat','HX_final','HZ_final','heights','ell');

HX_final = sparse(double(HX_final));
HZ_final = sparse(double(HZ_final));

C = HX_final*HZ_final.';

n  = size(HX_final,2);
nX = size(HX_final,1);
nZ = size(HZ_final,1);
wX = max(sum(HX_final~=0,2));
wZ = max(sum(HZ_final~=0,2));
qX = max(sum(HX_final~=0,1));
qZ = max(sum(HZ_final~=0,1));

fprintf('ell=%d\n',ell);
fprintf('n=%d, nX=%d, nZ=%d\n',n,nX,nZ);
fprintf('wX=%d, wZ=%d, qX=%d, qZ=%d\n',wX,wZ,qX,qZ);
fprintf('Integer commutator: nnz=%d\n',nnz(C));

assert(ell==7);
assert(isequal([n,nX,nZ,wX,wZ,qX,qZ],[3938,1886,2052,6,6,5,3]));
assert(nnz(C)==0);
assert(all(ismember(nonzeros(HX_final),[-1,1])));
assert(all(ismember(nonzeros(HZ_final),[-1,1])));

fprintf('All condition-aware BB72 checks passed.\n');
