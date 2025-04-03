param_file = 'params_male.txt';

numchck = 82;
lines = readlines(param_file);
if size(lines,1) ~= numchck
    error("incorrect number of parameters")
end

param_names = strings([1,numchck]);
param_vals = zeros(1,numchck);
for i=1:numchck
    param_names(i) = extractBefore(lines(i), "=");
    val_parse = extractAfter(lines(i), "=");
    pval = extractBefore(val_parse, ";");
    param_vals(i) = str2double(extractAfter(pval, " "));
end