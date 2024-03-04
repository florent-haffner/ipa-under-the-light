function [tcal,pcal,btcal,bxcal, RMSEcal,SRcal,levcal,wgcal]=fwFPLS(Xcal,ycal,nbfact)

%***************************************
% correlation PLS
%***************************************
%[tcal,pcal,btcal,RMSEcal,SRcal,levcal,wgcal]=fwFPLS(Xcal,ycal,nbfact);
% tcal		= scores de calibration
% pcal		= loading calibration
% btcal		= coefficients par rapport aux scores t
% bxcal		= coefficients par rapport aux variables x originales
% RMSEcal   = RMSEP sur la base de calibration
% SRcal     = norme des spectres r�siduels
% levcal    = leverage
% wgcal		= weigth calibration

% Initialisations
nbcal=size(Xcal,1);

% on centre y et on norme 
ycalm=mean(ycal);
ycalc=ycal-ycalm;

% on centre X
Xcalm=mean(Xcal);
Xcalc=Xcal-ones(nbcal,1)*Xcalm;

[bxcal,tcal,U,C,wgcal,pcal]=Plssimfw(Xcalc,ycalc,nbfact);	

%clear BPLS U C

% r�sultats de calibration
sel=[1:nbfact];
btcal=inv(tcal(:,sel)'*tcal(:,sel))*tcal(:,sel)'*ycalc;
% ou btcal=tcal(:,sel)'*ycalc;

% r�sidus calibration
Rcal=Xcalc-tcal*wgcal';
% racine de somme des r�sidus au carr� 
SRcal=sqrt(sum(Rcal'.*Rcal'));

% leverage
% norme carr�e des tcal
stt=sum(tcal.*tcal);
% part de chaque CP sur stt
levcal=tcal.*tcal./(ones(size(tcal,1),1)*stt);
levcal=sum(levcal');

%***************************************
% pr�diction sur la base de calibration
%***************************************

yhatcal=(tcal(:,sel)*btcal)+ycalm;
% ou yhatcal=Xcalc*bxcal+ycalm;
RMSEcal=RMSEP(ycal,yhatcal);

return

