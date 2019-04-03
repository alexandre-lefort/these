function sym_val = symtbx_convert_double_to_sym(val, prec)

    res = sprintf(['%.' num2str(prec) 'f'], val);
    while ((res(end) == '0') && (length(res) > 1))
        res = res(1:end-1);
    end
    point = strfind(res,'.');
    if (~isempty(point) && eval(res) ~= 0)
        res = erase(res,'.');
        den = 10^(length(res) - point + 1);
        if (den > 1)
            res = [res '/' num2str(den)];
        end
        while ((res(1) == '0') && (length(res) > 1))
            res = res(2:end);
        end
        if (res(1) == '-')
            while ((res(2) == '0') && (length(res) > 1))
                res = ['-' res(3:end)];
            end
        end
    end
    if (eval(res) == 0)
        res = '0';
    end
    sym_val = sym(res);

    
    