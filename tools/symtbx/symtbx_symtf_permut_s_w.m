function  res = symtbx_symtf_permut_s_w(expr)

    syms s;
    %%w = sym('r','real');
    %%res1 = simplify(subs(expr, 's',  i*w));
    res2 = simplify(subs(expr, 'w', -i*s));

    %%if (isequal(res1, expr))
        res = res2;
    %%elseif (isequal(res2, expr))
    %%    res = res1;
    %%else
    %%    res = expr;
    %%end




