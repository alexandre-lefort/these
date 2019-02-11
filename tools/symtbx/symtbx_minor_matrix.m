%% Give reduced matrix 
%%

function M_minor = symtbx_minor_matrix(M, index)
	M_minor = M(1:index,1:index);