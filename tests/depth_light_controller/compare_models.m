model = load_model('iver2.nav');

Vs = 1.5;

r.toto = 0.0;
[A,B,C,D] = modele_lineaire(Vs, model, 'IMM',r);
        
sys_pb = ss(A,B,C,D);
        
% Augmentation vers modele [IZ Z Theta Pente Q]
A = [0 1 0 0 0 ; 0 0 0 -Vs 0 ; zeros(3,2) A];
B = [0 ; 0 ; B];
C = eye(5);
D = zeros(5,1);

ssm = ss(A,B,C,D);

model_sym = build_model_lin(1.5);


disp(model_sym.sys0.a-A(2:end,2:end));
disp(model_sym.sys0.b-B(2:end));
