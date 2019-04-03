function s_sym = symtbx_convert_struct_to_sym(s,prec)

    f = fields(s);
    n = length(f);
    
    for ii=1:n
        if (length(s.(f{ii})) <= 1)
            s_sym.(f{ii}) = symtbx_convert_double_to_sym(s.(f{ii}),prec);
        end
    end
    
    