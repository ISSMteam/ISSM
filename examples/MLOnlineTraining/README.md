# Native neural basal friction training

Use an ISSM build with CoDiPack and a compatible linear solver. 

This exmaple requires an ISSM configured with CoDiPack!

```matlab
runme('stressbalance')
runme('transient')
```

An optional second argument selects the MPI rank count. A third argument
sets the random seed (default 1). The driver writes
the objective history, fitted weights, normalization, and model outputs to
`Results/nn_<stage>_np<ranks>.mat`.

The network is a native 2-32-32-1 ReLU MLP with one shared row of 1,185
weights and biases. Inputs are the square root of the squared friction
coefficient and the square root of basal speed after smoothing near zero.
The denormalized network output is squared to produce the drag prefactor.
The driver uses seeded random weights and biases. It uses the reference
drag to set an output scale and starts near a constant drag. The fixed
speed targets come from Weertman with m=3. This synthetic case uses a uniform friction coefficient, so it
tests learning the speed dependence at that coefficient.

Weertman is used for the reference solves. The network supplies drag during
each optimization solve, and M1QN3 changes the shared weights between
objective evaluations. The stressbalance case trains on one stress-balance solve.
The transient case sums surface-speed misfits from two steps with mass
transport enabled. Both cases compare a CoDiPack derivative with a centered
finite difference through the full solve before fitting, then reload the
fitted weights for an ordinary forward run.
