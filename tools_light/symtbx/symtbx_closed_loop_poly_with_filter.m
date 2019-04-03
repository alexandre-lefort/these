function [poly, mat] = symtbx_closed_loop_poly_with_filter(sys, sys_ctrl, sys_filter)


    %% return characteristic polynomial of a closed-loop system

    A = sys.a;
    B = sys.b;
    C = sys.c;
    %%D = sys.D;

    Ak = sys_ctrl.a;
    Bk = sys_ctrl.b;
    Ck = sys_ctrl.c;
    Dk = sys_ctrl.d;

    Af = sys_filter.a;
    Bf = sys_filter.b;
    Cf = sys_filter.c;

    [nf,mf] = size(Af);
    [nk,mk] = size(Ak);
    [n ,m ] = size(A);

    mat  = [A             (-B*Dk*Cf) (B*Ck)       ; ...
            Bf*C          (Af)       zeros(nf,mk) ; ...
            zeros(nk,m)   (Bk*Cf)    (Ak)        ];

    syms s;

    [n,m] = size(mat);
    poly  = simplify(det(s*eye(n) - mat));
