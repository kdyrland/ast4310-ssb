function [ V ] = voigt( a,v )


z = v + 1i.*a;
V = fadf(z);


end

