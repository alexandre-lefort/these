function symtf = symtbx_symss_to_symtf(ss_t)


    %% Return the transfert matrix 
    %% from a state space representation
    %%
    %% Inputs   : - ss_t : structure with (a,b,c,d) as
    %%                     symbolic matrices
    %%
    %% Outputs  : - symtf : symbolic transfert matrix


    syms s;

    A = ss_t.a;
    B = ss_t.b;
    C = ss_t.c;
    D = ss_t.d;

    [n,m] = size(A);
    symtf = C*inv(s*eye(n)-A)*B + D;
