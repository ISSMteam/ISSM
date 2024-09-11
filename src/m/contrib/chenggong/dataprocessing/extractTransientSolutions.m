function [transientSolutions]=extractTransientSolutions(md)
%extractTransientSolutions - take the transient solutions out from model
%                            and put each of them into an individual array.
transientSolutions.time = cell2mat({md.results.TransientSolution(:).time});
transientSolutions.vx = cell2mat({md.results.TransientSolution(:).Vx});
transientSolutions.vy = cell2mat({md.results.TransientSolution(:).Vy});
transientSolutions.vel = cell2mat({md.results.TransientSolution(:).Vel});
transientSolutions.volume = cell2mat({md.results.TransientSolution(:).IceVolume});
transientSolutions.thickness = cell2mat({md.results.TransientSolution(:).Thickness});
if (isfield(md.results.TransientSolution, 'SigmaVM'))
	transientSolutions.SigmaVM = cell2mat({md.results.TransientSolution(:).SigmaVM});
end
if (isfield(md.results.TransientSolution, 'SmbMassBalance'))
	transientSolutions.smb = cell2mat({md.results.TransientSolution(:).SmbMassBalance});
end
transientSolutions.ice_levelset = cell2mat({md.results.TransientSolution(:).MaskIceLevelset});
if (isfield(md.results.TransientSolution, 'CalvingCalvingrate'))
	transientSolutions.calvingRate = cell2mat({md.results.TransientSolution(:).CalvingCalvingrate});
end
if (isfield(md.results.TransientSolution, 'CalvingMeltingrate'))
	transientSolutions.meltingRate = cell2mat({md.results.TransientSolution(:).CalvingMeltingrate});
end
