function [ gammavdw ] = gammavdw_NaD(temp, pgas, s)

rsq_u = rsq_NaD(s-1);
rsq_l = rsq_NaD(2);

loggvdw = 6.33 + 0.4 .* log10(rsq_u - rsq_l) + log10(pgas) - 0.7 .* log10(temp);
gammavdw = 10.^loggvdw;


end
