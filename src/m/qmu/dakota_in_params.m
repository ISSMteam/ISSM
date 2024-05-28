%
%  populate a Dakota parameter structure.
%
%  [params]=dakota_in_params(params)
%
%  where the optional input is:
%    params        (structure array, method-independent parameters)
%
%  and the output is the same.
%
%  this function takes a structure of method-independent dakota
%  parameters, which may be empty, and adds default parameters
%  for those parameters which do not exist.
%
%  the field names of the structure are identical to the dakota
%  parameter names (and are in fact used to write them to the
%  files).  logical values are used for parameters which have
%  no associated data and are determined only by their presence
%  or absence.
%
%  note that the method-dependent parameters are contained in
%  the dakota_method class object.
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
function [params]=dakota_in_params(params)

if ~nargin
    help dakota_in_params
    return
end

%%  process the input parameters

if ~exist('params','var')
    params=struct();
end

%%  strategy section

if ~isfield(params,'graphics')
    params.graphics=false;
end
if ~isfield(params,'tabular_graphics_data')
    params.tabular_graphics_data=false;
end
% could use unique file name rather than 'dakota_tabular.dat'
if ~isfield(params,'tabular_graphics_file')
    params.tabular_graphics_file=false;
end

%%  method section

%  nearly all method parameters are in the dakota_method class
%  or result from the response level lists

if ~isfield(params,'compute')
    params.compute='probabilities';
end
if ~isfield(params,'distribution')
    params.distribution='cumulative';
end

%%  model section

%%  interface section

if ~isfield(params,'system')
    params.system=false;
end
if ~isfield(params,'fork')
    params.fork=false;
end
if ~isfield(params,'direct')
    params.direct=false;
end

%  interface parallelism controls

if ~isfield(params,'asynchronous')
    params.asynchronous=true;
end
if ~isfield(params,'evaluation_concurrency')
    params.evaluation_concurrency=false;
end
if ~isfield(params,'analysis_concurrency')
    params.analysis_concurrency=false;
end
if ~isfield(params,'evaluation_servers')
    params.evaluation_servers=false;
end
if ~isfield(params,'evaluation_self_scheduling')
    params.evaluation_self_scheduling=false;
end
if ~isfield(params,'evaluation_static_scheduling')
	params.evaluation_static_scheduling=true;
end
if ~isfield(params,'evaluation_scheduling')
	params.evaluation_scheduling=false;
end
if ~isfield(params,'processors_per_evaluation')
	params.processors_per_evaluation=false;
end
if ~isfield(params,'analysis_servers')
    params.analysis_servers=false;
end
if ~isfield(params,'analysis_self_scheduling')
    params.analysis_self_scheduling=false;
end
if ~isfield(params,'analysis_static_scheduling')
    params.analysis_static_scheduling=false;
end

%  algebraic mappings

if ~isfield(params,'algebraic_mappings')
    params.algebraic_mappings=false;
end

%  simulation interface controls

if ~isfield(params,'analysis_driver')
    params.analysis_driver='';
end
if ~isfield(params,'analysis_components')
    params.analysis_components='';
end
if ~isfield(params,'input_filter')
    params.input_filter='';
end
if ~isfield(params,'output_filter')
    params.output_filter='';
end

if ~isfield(params,'failure_capture')
    params.failure_capture='abort';
end
if ~isfield(params,'deactivate')
    params.deactivate='evaluation_cache restart_file';
end

%  system call or fork interface

if ~isfield(params,'parameters_file')
    params.parameters_file='params.in';
end
if ~isfield(params,'results_file')
    params.results_file='results.out';
end
if ~isfield(params,'verbatim')
    params.verbatim=false;
end
if ~isfield(params,'aprepro')
    params.aprepro=false;
end
if ~isfield(params,'file_tag')
    params.file_tag=true;
end
if ~isfield(params,'file_save')
    params.file_save=true;
end

%  direct function interface

if ~isfield(params,'processors_per_analysis')
    params.processors_per_analysis=false;
end

%%  responses section

if ~isfield(params,'numerical_gradients')
    params.numerical_gradients=false;
end
if ~isfield(params,'method_source')
    params.method_source='dakota';
end
if ~isfield(params,'interval_type')
    params.interval_type='forward';
end
if ~isfield(params,'fd_gradient_step_size')
    params.fd_gradient_step_size=0.001;
end
if ~isfield(params,'analytic_gradients')
    params.analytic_gradients=false;
end
%  mixed_gradients not fully implemented
if ~isfield(params,'mixed_gradients')
    params.mixed_gradients=false;
end
if ~isfield(params,'id_analytic_gradients')
    params.id_analytic_gradients=false;
end
if ~isfield(params,'id_numerical_gradients')
    params.id_numerical_gradients=false;
end
%  hessians not fully implemented
if ~isfield(params,'numerical_hessians')
    params.numerical_hessians=true;
end
if ~isfield(params,'hessian_gradient_step_size')
    params.hessian_gradient_step_size=0.001;
end

end
