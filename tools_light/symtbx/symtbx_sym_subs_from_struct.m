function  res = symtbx_sym_subs_from_struct(expr, p)


    %% Takes a symbolic expression and substitutes every
    %% parameters from p in expr with its value given in p.
    %%
    %% Inputs   : - expr : symbolic expression
    %%            - p    : one level matlab structure
    %%
    %% Outputs  : - res  : symbolic expression
    %%
    %%
    %% Example :
    %%
    %% syms a b;
    %% expr = a + b;
    %% p.a  = 1;
    %%
    %% res  = symtbx_sym_subs_from_struct(expr,p);
    %%
    %% -> res = 1 + b;


    pf   = fields(p);
    nb_p = length(pf);

    res = expr;

    for ii=1:nb_p
        res = subs(res, pf{ii}, p.(pf{ii}));
    end
