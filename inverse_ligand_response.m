function out = inverse_ligand_response(K,p)
% c = logspace(10,30,1e2); % DR=4 fold

%out = K/((p^-1)-1); 
out = K * (p / (1-p));