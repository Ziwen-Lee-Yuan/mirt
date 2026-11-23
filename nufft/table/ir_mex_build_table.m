% ir_mex_build_table
% run matlab's "mex" command to "compile" the table interpolation code
% into mex files.

dir_current = pwd;
% parent of this file's folder (nufft/table/), not path_find_dir('nufft'),
% which any other folder ending in "nufft" earlier on the path would shadow.
dir_nufft = fileparts(fileparts(mfilename('fullpath')));
dir_table = [dir_nufft filesep 'table'];
cd(dir_table)

% compile mex files using ir_mex_fun2 ("2" since there are two arguments)
fun = ir_mex_fun();
fun('interp1_table_adj_mex.c',	'interp1_table1_adj.c')
fun('interp1_table_mex.c',	'interp1_table1_for.c')
fun('interp2_table_adj_mex.c',	'interp2_table1_adj.c')
fun('interp2_table_mex.c',	'interp2_table1_for.c')
fun('interp3_table_adj_mex.c',	'interp3_table1_adj.c')
fun('interp3_table_mex.c',	'interp3_table1_for.c')


cd(dir_current)

% "test" mex files by running them (each should display usage)
try
	interp1_table_adj_mex
	interp1_table_mex
    interp2_table_adj_mex
    interp2_table_mex
    interp3_table_adj_mex
    interp3_table_mex
catch
    disp(['WARNING! Error compiling mex files in ', mfilename])
end

