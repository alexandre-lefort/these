function noise = load_noises()

deg = pi/180;

%% Modèle de bruit utilisé pour les réglages des estimateurs et des simulateurs

noise.Phi    = 0.0138;
noise.Theta  = 0.0125;
noise.Psi    = 0.0248;
noise.P      = 0.0015; 
noise.Q      = 0.0010;
noise.R      = 0.0010;
noise.NT1    = 0.01  ;
noise.U      = 0.2265;
noise.V      = 0.2268;
noise.Beta   = 0.2 * deg ;
noise.Alpha  = 0.2 * deg ;
noise.Z      = 0.06      ;
