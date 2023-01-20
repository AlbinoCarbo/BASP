% Equat2Galactic, function for selecting astrometric positions
% of satellites outside the belt of the Milky Way. Useful to
% avoid satellite astrometry in crowded star fields
% where accurate position measurements are difficult.
%
% Algoritmo: carica le posizioni astrometriche osservate, trasforma AR e DEC in coordinate galattiche Lambda e Beta,
% elimina le osservazioni che stanno dento la fascia della Via Lattea (tranne quella invernale anche se è possibile attivarla) 
% e restituisce le altre. Se tutte le osservazioni astrometriche vengono cancellate, quindi è impossibile calcolare l'orbita, 
% avvisa l'utente con un warning e non le filtra restituendo tutto l'input. 
%
% INPUT:
% anno = vettore anno di osservazione;
% mese = vettore mese di osservazione;
% giorno = vettore giorno di osservazione;
% ora = vettore ora di osservazione;
% minuto = vettore minuto di osservazione;
% secondo = vettore secondo di osservazione;
% alpha = vettore RA astrometrico (gradi);
% delta = vettore DEC astrometrico (gradi);
% norad = vettore nome satellite;
%
% OUTPUT:
% Anno = vettore anno di osservazione;
% Mese = vettore mese di osservazione;
% Giorno = vettore giorno di osservazione;
% Ora = vettore ora di osservazione;
% Minuto = vettore minuto di osservazione;
% Secondo = vettore secondo di osservazione;
% Alpha = vettore RA astrometrico fuori dalla Via Lattea (gradi);
% Delta = vettore DEC astrometrico fuori dalla Via Lattea (gradi);
% Norad = vettore nome satellite;
%
% Versione del 15 ottobre 2021
% Albino Carbognani, INAF-OAS
 
function [Anno, Mese, Giorno, Ora, Minuto, Secondo, Alpha, Delta, Norad] = equat2galactic(anno, mese, giorno, ora, minuto, secondo, alpha, delta, norad)

% Definizione costanti
alpha_NGP=192.85; % Ascensione retta polo nord galattico (J2000.0), gradi
delta_NGP=27.13;  % Declinazione polo nord galattico (J2000.0), gradi
L_NCP=122.89167;  % Longitudine galattica polo celeste nord, gradi

% Compute trigonometric combinations of coordinates
sb = sind(delta_NGP)*sind(delta) + cosd(delta_NGP)*cosd(delta).*cosd(alpha-alpha_NGP);
Y = cosd(delta).*sind(alpha-alpha_NGP);
X = cosd(delta_NGP)*sind(delta) - sind(delta_NGP)*cosd(delta).*cosd(alpha-alpha_NGP);

% Compute galactic coordinates (degrees)
lambda = L_NCP-(180/pi)*atan2(Y,X);

for i=1:length(lambda)
  if lambda(i)<0
      lambda(i)=lambda(i)+360;
  end

  if lambda(i) > 360
      lambda(i)=lambda(i)-360;
  end
end

beta=(180/pi)*asin(sb);

% Select equatorial coordinates outside Milky Way only
k=1;
for i=1:length(delta)
   
    if (lambda(i)>=0) && (lambda(i)<=109) % Via Lattea fra Cefeo e Sagittario
       if (beta(i) >13) || (beta(i) <-13)
          Alpha(k)=alpha(i); Delta(k)=delta(i);
          Anno(k)=anno(i); Mese(k)=mese(i); Giorno(k)=giorno(i); 
          Ora(k)=ora(i); Minuto(k)=minuto(i); Secondo(k)=secondo(i); 
          Norad(k)=norad(i);
          k=k+1;      
       end        
    end
    
    
    if (lambda(i)>109) && (lambda(i)<=163) % Via Lattea fra Perseo a Cassiopea
       if (beta(i) >5) || (beta(i) <-12)
          Alpha(k)=alpha(i); Delta(k)=delta(i);
          Anno(k)=anno(i); Mese(k)=mese(i); Giorno(k)=giorno(i); 
          Ora(k)=ora(i); Minuto(k)=minuto(i); Secondo(k)=secondo(i); 
          Norad(k)=norad(i);
          k=k+1; 
       end        
    end
    
%     if (lambda(i)>163) && (lambda(i)<=340) % Via Lattea (invernale) fra Auriga e Squadra 
%        if (beta(i) >8) || (beta(i) <-8)
%           Alpha(k)=alpha(i); Delta(k)=delta(i);
%           Anno(k)=anno(i); Mese(k)=mese(i); Giorno(k)=giorno(i); 
%           Ora(k)=ora(i); Minuto(k)=minuto(i); Secondo(k)=secondo(i); 
%           Norad(k)=norad(i);
%           k=k+1;  
%        end        
%     end
    
    if (lambda(i)>340) && (lambda(i)<360) % Via Lattea fra Scorpione e Ofiuco
       if (beta(i) >20) || (beta(i) <-20)
          Alpha(k)=alpha(i); Delta(k)=delta(i);
          Anno(k)=anno(i); Mese(k)=mese(i); Giorno(k)=giorno(i); 
          Ora(k)=ora(i); Minuto(k)=minuto(i); Secondo(k)=secondo(i); 
          Norad(k)=norad(i);
          k=k+1;  
       end        
    end
    
end

if k==1
   disp('     Equat2Galactic: no observations survived after Milky Way filter!');
   disp('   ')
   disp('     Equat2Galactic WARNING: Filtering canceled - manual orbit computation recommended')
   disp('   ')
   
   % Ritorna tutte le osservazioni di input in output
   Anno=anno; Mese=mese; Giorno=giorno; Ora=ora; Minuto=minuto; 
   Secondo=secondo; Alpha=alpha; Delta=delta; Norad=norad;
end

end
