function [poly, mat] = symtbx_closed_loop_poly(sys, sys_ctrl)


    %% return characteristic polynomial of a closed-loop system


    A = sys.a;
    B = sys.b;
    C = sys.c;
    %%D = sys.D;

    Ak = sys_ctrl.a;
    Bk = sys_ctrl.b;
    Ck = sys_ctrl.c;
    Dk = sys_ctrl.d;

    if (sum(sum(Ak)) == 0 && sum(sum(Ck)) == 0)
        mat = (A+B*(-Dk)*C);
    else
        mat = [(A+B*(-Dk)*C) (B*Ck) ; ...
               (-Bk*C)       (Ak)  ];
    end

    syms s;

    [n,m] = size(mat);
    poly  = simplify(det(s*eye(n) - mat));
