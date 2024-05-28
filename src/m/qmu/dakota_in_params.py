from dakota_in_data import *
#move this later:
from helpers import *


def dakota_in_params(params):
    '''
  populate a Dakota parameter structure.

  params = dakota_in_params(params)

  where the optional input is:
    params        (structure array, method - independent parameters)

  and the output is the same.

  this function takes a structure of method - independent dakota
  parameters, which may be empty, and adds default parameters
  for those parameters which do not exist.

  the field names of the structure are identical to the dakota
  parameter names (and are in fact used to write them to the
  files).  logical values are used for parameters which have
  no associated data and are determined only by their presence
  or absence.

  note that the method - dependent parameters are contained in
  the dakota_method class object.
'''
    if params is None:
        help(dakota_in_params)
        return

    #  process the input parameters
    if len(fieldnames(params)) == 0:
        params = struct()

    #  strategy section
    if not isfield(params, 'graphics'):
        params.graphics = False

    if not isfield(params, 'tabular_graphics_data'):
        params.tabular_graphics_data = False

    # could use unique file name rather than 'dakota_tabular.dat'
    if not isfield(params, 'tabular_graphics_file'):
        params.tabular_graphics_file = False

    #  method section
    #  nearly all method parameters are in the dakota_method class
    #  or result from the response level lists
    if not isfield(params, 'compute'):
        params.compute = 'probabilities'

    if not isfield(params, 'distribution'):
        params.distribution = 'cumulative'

    #  model section

    #  interface section
    if not isfield(params, 'system'):
        params.system = False

    if not isfield(params, 'fork'):
        params.fork = False

    if not isfield(params, 'direct'):
        params.direct = False

    #  interface parallelism controls
    if not isfield(params, 'asynchronous'):
        params.asynchronous = True

    if not isfield(params, 'evaluation_concurrency'):
        params.evaluation_concurrency = False

    if not isfield(params, 'analysis_concurrency'):
        params.analysis_concurrency = False

    if not isfield(params, 'evaluation_servers'):
        params.evaluation_servers = False

    if not isfield(params, 'evaluation_self_scheduling'):
        params.evaluation_self_scheduling = False

    if not isfield(params, 'evaluation_static_scheduling'):
        params.evaluation_static_scheduling = True

    if not isfield(params, 'evaluation_scheduling'):
        params.evaluation_scheduling = False

    if not isfield(params, 'processors_per_evaluation'):
        params.processors_per_evaluation = False

    if not isfield(params, 'analysis_servers'):
        params.analysis_servers = False

    if not isfield(params, 'analysis_self_scheduling'):
        params.analysis_self_scheduling = False

    if not isfield(params, 'analysis_static_scheduling'):
        params.analysis_static_scheduling = False

    #  algebraic mappings
    if not isfield(params, 'algebraic_mappings'):
        params.algebraic_mappings = False

    #  simulation interface controls
    if not isfield(params, 'analysis_driver'):
        params.analysis_driver = ''

    if not isfield(params, 'analysis_components'):
        params.analysis_components = ''

    if not isfield(params, 'input_filter'):
        params.input_filter = ''

    if not isfield(params, 'output_filter'):
        params.output_filter = ''

    if not isfield(params, 'failure_capture'):
        params.failure_capture = 'abort'

    if not isfield(params, 'deactivate'):
        params.deactivate = 'evaluation_cache restart_file'

    #  system call or fork interface
    if not isfield(params, 'parameters_file'):
        params.parameters_file = 'params.in'

    if not isfield(params, 'results_file'):
        params.results_file = 'results.out'

    if not isfield(params, 'verbatim'):
        params.verbatim = False

    if not isfield(params, 'aprepro'):
        params.aprepro = False

    if not isfield(params, 'file_tag'):
        params.file_tag = True

    if not isfield(params, 'file_save'):
        params.file_save = True

    #  direct function interface
    if not isfield(params, 'processors_per_analysis'):
        params.processors_per_analysis = False

    #  responses section
    if not isfield(params, 'numerical_gradients'):
        params.numerical_gradients = False

    if not isfield(params, 'method_source'):
        params.method_source = 'dakota'

    if not isfield(params, 'interval_type'):
        params.interval_type = 'forward'

    if not isfield(params, 'fd_gradient_step_size'):
        params.fd_gradient_step_size = 0.001

    if not isfield(params, 'analytic_gradients'):
        params.analytic_gradients = False

    #  mixed_gradients not fully implemented
    if not isfield(params, 'mixed_gradients'):
        params.mixed_gradients = False

    if not isfield(params, 'id_analytic_gradients'):
        params.id_analytic_gradients = False

    if not isfield(params, 'id_numerical_gradients'):
        params.id_numerical_gradients = False

    #  hessians not fully implemented
    if not isfield(params, 'numerical_hessians'):
        params.numerical_hessians = True

    if not isfield(params, 'hessian_gradient_step_size'):
        params.hessian_gradient_step_size = 0.001

    return params
