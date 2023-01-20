% Time seconds formatting script
% in the format XX.XXXXX.
% Input is a number, output is a string.
%
% INPUT
% SEC = numero dei secondi con i decimali
%
% OUTPUT
% Stringa dei secondi nel formato XX.XXX
%
% Albino Carbognani, INAF-OAS
% Versione del 18 ottobre 2021


function [Time_s] = format_seconds(SEC)

SEC_int=fix(SEC);                              % Parte intera secondi
SEC_dec=fix(100000*(SEC-SEC_int));             % Parte con 5 decimali dei secondi 
str_SEC_int=num2str(SEC_int,'%02.f');          % Parte intera dei secondi con zero davanti
str_SEC_dec=num2str(SEC_dec,'%05.f');          % Parte decimale dei secondi
Time_s=strcat(str_SEC_int, '.',str_SEC_dec);   % Stringa dei secondi nel formato XX.XXXXX

end
