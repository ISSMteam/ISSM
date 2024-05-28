%
%  least squares study for rosenbrock case
%  (see Users4.2.pdf, Sec. 2.4.1)
%
%  "Copyright 2009, by the California Institute of Technology.
%  ALL RIGHTS RESERVED. United States Government Sponsorship
%  acknowledged. Any commercial use must be negotiated with
%  the Office of Technology Transfer at the California Institute
%  of Technology.  (J. Schiermeier, NTR 47078)
%
%  This software may be subject to U.S. export control laws.
%  By accepting this  software, the user agrees to comply with
%  all applicable U.S. export laws and regulations. User has the
%  responsibility to obtain export licenses, or other export
%  authority as may be required before exporting such information
%  to foreign countries or providing access to foreign persons."
%

function [dout,ddat]=dakota_rosenbrock_ls()

%  define dakota variables as continuous design
dvar.x1=continuous_design('',-1.2,-2,2);
dvar.x2=continuous_design('', 1.0,-2,2);

%  define dakota response as least-squares terms
dresp.f1sq=least_squares_term('');
dresp.f2sq=least_squares_term('');

%  define dakota method and specify method-dependent parameters
dmeth=dakota_method('nl2sol');
dmeth=dmeth_params_set(dmeth,'max_iterations',100,...
                             'convergence_tolerance',1.e-4);

%  specify method-independent parameters
%  (dakota_in_params does not need to be called, but provides a template)
dparams=dakota_in_params([]);
dparams.direct=true;
dparams.analysis_driver='rosenbrock';
dparams.tabular_graphics_data=true;
dparams.tabular_graphics_file='dakota_rosenbrock_ls.dat';
dparams.analytic_gradients=true;

%  write out dakota input file
dakota_in_write(dmeth,dvar,dresp,dparams,'dakota_rosenbrock_ls.in')

%  execute dakota
!dakota -i dakota_rosenbrock_ls.in -o dakota_rosenbrock_ls.out

%  read dakota output and tabular data files
%  (output file for parameter studies has no interesting info)
[method,dout]=dakota_out_parse('dakota_rosenbrock_ls.out');
[~     ,ddat]=dakota_out_parse('dakota_rosenbrock_ls.dat');

%  perform any desired plotting
% plot_rvsv_surf(ddat,{'x1','x2'},ddat,{'f'})

end

