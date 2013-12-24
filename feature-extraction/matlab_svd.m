function matlab_svd(file_loc, k)
load(file_loc);
[U, S, V] = svds(avgOP, str2num(k));
save(file_loc, 'U', 'S', 'V');
exit

