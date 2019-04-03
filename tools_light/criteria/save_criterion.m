function save_criterion(c, k, p, name)

	%% Create a new file;
	file = fopen(name,'w');
	
	kv  = fields(k);     nkv = size(kv,1); %% Structure of tunable   parameters
	pv  = fields(p);     npv = size(pv,1); %% Structure of uncertain parameters

	%% Write head of file	
	mes = '#include "ibex.h" \n \n using namespace std; \n using namespace ibex; \n \n';
	mes = strcat(mes, 'Variable x(' , num2str(nkv) , '); \nVariable p(' , num2str(npv) , '); \n \n');
	

	%% Write function
	n = size(c,2);
	
	for ii = 1:n

		mes = strcat(mes, 'const ExprNode& stab', num2str(n+1-ii), ' = ');
		expr = char(c{ii});
			
		%% Substitute generic names for parameters and gains (xi,pi) :
		
		for jj=1:nkv
			char_x = strcat('x[',num2str(jj-1),']');
			expr   = strrep(expr ,kv{jj}, char_x);
		end

		for jj=1:npv
			char_p = strcat('p[',num2str(jj-1),']');
			expr   = strrep(expr ,pv{jj}, char_p);
		end
		
		mes = strcat(mes,expr, '; \n');
	
	end

	fprintf(file,mes);
	fclose(file);