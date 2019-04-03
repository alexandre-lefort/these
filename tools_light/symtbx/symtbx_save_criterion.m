function symtbx_save_criterion(c, k, p, name)

    exp_max = 20;

    file = fopen(name,'w');
	
    kv  = fields(k);  nkv = length(kv);
    pv  = fields(p);  npv = length(pv);

    %% Print parameters and gains
    disp(' '); disp('-> gains :');
    for ii=1:nkv
        disp([num2str(ii), ' : ', char(kv{ii})]);
    end

    disp(' '); disp('-> parameters :');
    for ii=1:npv
        disp([num2str(ii), ' : ', char(pv{ii})]);
    end

    %% Write head of function
    if (isa(c,'cell'))
        mes = ['function f(x[', num2str(nkv), '], p[', num2str(npv), ']) \n return('];
    else
        mes = ['function f(x[', num2str(nkv), '], p[', num2str(npv+1), ']) \n return('];
    end

    %% Write all polys
		
    n = length(c);
	
    for ii = 1:n
        
        if (isa(c,'cell'))
           expr = char(c{ii});
        else 
           expr = char(c);
        end
			
        %% Substitute generic names for parameters and gains (xi,pi) :	
        for jj=1:nkv
            char_x = ['x[', num2str(jj-1), ']'];
            expr   = strrep(expr, kv{jj}, char_x);
        end

       for jj=1:npv
            char_p = ['p[', num2str(jj-1), ']'];
            expr   = strrep(expr, pv{jj}, char_p);
       end
	
        %% Substitute w with logarithm :
        for jj=exp_max:(-1):1
            if jj==1
               char_w = 'w';
            else
                char_w = ['w^',num2str(jj)];
            end
            char_p = ['exp(ln(10)*',num2str(jj),'*p[',num2str(npv),'])'];
            expr = strrep(expr, char_w, char_p);
        end

        mes = [mes,expr];
        if ii<n
            mes = [mes,','];
        end
    end

    mes = [mes,') \n end'];
			
    fprintf(file,mes);
    fclose(file);
