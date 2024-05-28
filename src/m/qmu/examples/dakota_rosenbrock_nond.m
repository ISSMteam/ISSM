%
%  monte carlo sampling for rosenbrock case
%  (see Users4.2.pdf, Sec. 2.4.9)
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

function [dout,ddat,scm,pcm,srcm,prcm]=dakota_rosenbrock_nond()

%  define dakota variables as uniform uncertain
%  (may use set 1 or set 2 of variables, but not both)
dvar(1).uuv(1)=uniform_uncertain('x1',-2,2);
dvar(1).uuv(2)=uniform_uncertain('x2',-2,2);
dvar(2).x1=uniform_uncertain('',-2,2);
dvar(2).x2=uniform_uncertain('',-2,2);

%  define dakota response as response function
%  (may use set 1 or set 2 of responses, but not both)
dresp(1).rf=response_function('f',[100]);
dresp(2).f=response_function('',[100]);

%  define dakota method and specify method-dependent parameters
dmeth=dakota_method('nond_samp');
dmeth=dmeth_params_set(dmeth,'samples',200,...
                             'seed',17,...
                             'sample_type','random');

%  specify method-independent parameters
%  (dakota_in_params does not need to be called, but provides a template)
dparams=dakota_in_params([]);
dparams.direct=true;
dparams.analysis_driver='rosenbrock';
dparams.tabular_graphics_data=true;
dparams.tabular_graphics_file='dakota_rosenbrock_nond.dat';

%  write out dakota input file
dakota_in_write(dmeth,dvar(2),dresp(2),dparams,'dakota_rosenbrock_nond.in')

%  execute dakota
!dakota -i dakota_rosenbrock_nond.in -o dakota_rosenbrock_nond.out

%  read dakota output and tabular data files
[method,dout,scm,pcm,srcm,prcm]=dakota_out_parse('dakota_rosenbrock_nond.out');
[~     ,ddat                  ]=dakota_out_parse('dakota_rosenbrock_nond.dat');

%  perform any desired plotting
plot_boxplot(ddat,{'f'})
plot_normplot(ddat,{'f'})
plot_sampdist_bars(ddat,{'f'})
plot_normdist_bars(ddat,{'f'})
plot_hist_norm(ddat,{'f'})
plot_hist_norm(ddat,{'f'},'hmin',0,'hmax',200)
plot_hist_norm_ci(ddat,{'f'})
plot_hist_norm_ci(ddat,{'f'},'ciplt','line','cdfplt','off','ymin1',0,'ymax1',0.15)
plot_prob_bars(dout,{'f'})
plot_rlev_bars_ci(ddat,{'f'},'xtlrot',90)

end

