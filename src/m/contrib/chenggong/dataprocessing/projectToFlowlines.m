% Function to project values on the mesh to a flowline
%
%	md			-	ISSM model with mesh
%	pValue	-	data on the mesh
%	fx			-	x coordinates of the flowline
%	fy			-	y coordinates of the flowline
%

function valueC = projectToFlowlines(md, pValue, fx, fy)
    temp = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y,...
        pValue, fx, fy);
    temp = diag(temp); % not very efficient, but works
    valueC = temp(:)';
end
