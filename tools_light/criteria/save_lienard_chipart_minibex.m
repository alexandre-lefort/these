function save_lienard_chipart_minibex(order, name, epsilon)

	%% Create a new file;
	file = fopen('ibex/lc.mbx','w');

	%% Write head of file	
	mes = '';
	coefs = build_lienard_chipart_generic(order, name);
		
	res1 = char(-coefs{1}+epsilon);
	res2 = char(-coefs{2}+epsilon);

	res = strcat('max(',res1,',',res2,')');

	%% Write function
	n = size(coefs,2);
	
	disp(n);

	for ii =1:n
		disp(-coefs{ii});
	end

	for ii = 3:n
		resii = char(-coefs{ii}+epsilon);
		res = strcat('max(',res,',',resii,')');
	end

	mes = strcat(mes, res);

	fprintf(file,mes);
	fclose(file);