function [retour, nombre]=RMSEP(yexp,ycalc)

% RMSEP            a pour fonction de calculer la racine carr�e de 
%                   l'�cart quadratique moyen entre valeurs 
%                  exp�rimentales et valeurs calcul�es. Typiquement ce 
%                  crit�re permet de quantifier la qualit� d'un 
%                  ajustement
% OUTPUT
% retour           ='Root Mean Square Error' 
% 
% INPUT
% ycalc            = propriete calculee (vecteur) 
% yexp             = propriete experimentale (vecteur)

EP=yexp-ycalc;

nombre=max(size(EP));

if size(EP,1)==1
	retour=sqrt(EP*EP'/nombre);
else
	retour=sqrt(EP'*EP/nombre);
end
