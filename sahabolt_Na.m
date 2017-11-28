function [ sahabolt ] = sahabolt_Na(temp, eldens, ionstage, level)

sahabolt = saha_Na(temp, eldens, ionstage).*Boltz_Na(temp, ionstage, level);

end

