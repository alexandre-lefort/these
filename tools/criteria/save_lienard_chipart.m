function save_lienard_chipart(order, name, epsilon)

	%% Create a new file;
	file = fopen('ibex/lc.h','w');

	%% Write head of file	
	mes = '#include "tstab.h" \n \n using namespace std; \n using namespace ibex; \n \n';	
	
	coefs = build_lienard_chipart_generic(order, name);
		
	res1 = char(-coefs{1}+epsilon);
	res2 = char(-coefs{2}+epsilon);

	res = strcat('ibex::max(',res1,',',res2,')');

	%% Write function
	n = size(coefs,2);
	
	 n

	for ii =1:n
		-coefs{ii}
	end

	for ii = 3:n
		resii = char(-coefs{ii}+epsilon);
		res = strcat('ibex::max(',res,',',resii,')');
	end

	res = strcat('Function stab(x,p,',res,'); \n');

	order = order+1;
	exp_max = 10;
	for ii=exp_max:(-1):1
		for jj=1:order
			if ii>1
				char_1 = strcat('stab',num2str(jj),'^',num2str(ii));
				char_p = strcat('ibex::pow(stab',num2str(jj),',',num2str(ii),')');
				res = strrep(res, char_1, char_p);
			end
		end
	end

	mes = strcat(mes, res);

	fprintf(file,mes);
	fclose(file);