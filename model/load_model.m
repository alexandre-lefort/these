function param = load_model(file)

% unites
deg=pi/180;
rad=180/pi;
kt=1852/3600;
%%min=60;
%%hour=60*min;

% modele de manoeuvrabilite
param  = read_DYSCO_nav(file);
param.InertieEnG = false; % Les inerties DYSCO sont donnees au centre de carene

% required for sfun_smc
optional_dyscov2p3 = {'w0', 'CzR2Q2' , 'CmR2Q2' ,'CmR3'};
for i=1:length(optional_dyscov2p3)
    if not(isfield(param,upper(optional_dyscov2p3{i})))
        param.(upper(optional_dyscov2p3{i})) = 0;
        param.(optional_dyscov2p3{i}) = 0;
    end
end

% required for sfun_smc
if not(isfield(param,'CNWP')) % bug DYSCO v2.3
    param.CNWP=param.CNVP;
    param.CnWP=param.CNVP;
end

% verrue car les majuscules et minuscules sont non significatives en DYSCO
if isfield(param,'CRemp')
	param.Cremp=param.CRemp;
else
	param.CRemp=param.Cremp;
end

% Environnement

pairs = {
    'g' 9.81000042; % (m/s2) acceleration gravit.
    'rho' 1026.0;   % (kg/m3) masse volumique eau
    'VcN' 0;        % (m/s) courant - composante north/east/down
    'VcE' 0;
    'VcD' 0
};

param.czr2q2 = 0.0;
param.cmr2q2 = 0.0;
param.cmr3   = 0.0;
param.cnwp   = 0.0;
param.xarr= param.xar;
    
for i=1:size(pairs,1)
    param.(pairs{i,1}) = pairs{i,2};
end

	
param.Beta1Max   =   20*deg ;
param.Beta1Min   =  -20*deg ;
param.Beta2Max   =   20*deg ;
param.Beta2Min   =  -20*deg ;
param.AlphaMin   =  -20*deg ;
param.AlphaMax   =   20*deg ;

