
setup_BB_project()
HX_final = sparse(double(HX_final));
HZ_final = sparse(double(HZ_final));

opts = struct;
opts.maxTimePerSolve = 300;   % seconds for each optimization
opts.display = "off";
opts.saveWitnesses = true;

D = compute_css_distance_binary_milp(HX_final, HZ_final,'display','off','maxTimePerSolve',100,'saveWitnesses',true)