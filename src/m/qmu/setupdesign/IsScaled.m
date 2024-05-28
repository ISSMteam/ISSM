function status=IsScaled(variablename)
%ISSCALED decide whether a variable should be scaled or not. 

switch variablename,
case {'MaterialsRhoIce','MaterialsRhoSeawater','MaterialsHeatCapacity','MaterialsThermalConductivity','Gravity','MaxVel'},

	status=0;

case {'GeometryThickness','GeometrySurface','GeometryBed','FrictionCoefficient','MaterialsRheologyB','MaterialsRheologyBbar'},

	status=1;

case {'RiftsFriction'},

	status=2; %special treatment

otherwise
	error(['IsScaled error  message: could not find ' variablename]);
end
