function response=qmuresponse(models,results,processedresults,descriptor)
%QMURESPONSE - compute response function from model results.

if strcmpi(descriptor,'max_vel'),
	response=max(processedresults.vel);
elseif strcmpi(descriptor,'min_vel'),
	response=min(processedresults.vel);
elseif strcmpi(descriptor,'max_vx'),
	response=max(processedresults.vx);
elseif strcmpi(descriptor,'max_abs_vx'),
	response=max(abs(processedresults.vx));
elseif strcmpi(descriptor,'min_vx'),
	response=min(processedresults.vx);
elseif strcmpi(descriptor,'max_vy'),
	response=max(processedresults.vy);
elseif strcmpi(descriptor,'max_abs_vy'),
	response=max(abs(processedresults.vy));
elseif strcmpi(descriptor,'min_vy'),
	response=min(processedresults.vy);
elseif strncmpi(descriptor,'mass_flux',9),
	indx=str2int(descriptor(10:end));
	if isempty(indx) || ~indx
		indx=1;
	end

	%call mass flux module.
	m_dh=models.dh;
	m_dhu=models.dhu;
	m_ds=models.ds;
	isSIA=m_dhu.parameters.isSIA;
	isSSA=m_dh.parameters.isSSA;
	isHO=m_dh.parameters.isHO;
	isFS=m_ds.parameters.isFS;
	if isSIA,

% for now, separate all segments from double array for parallel to make cells
		if (length(m_dhu.parameters.qmu_mass_flux_num_segments) > 1)
			segments=m_dhu.parameters.qmu_mass_flux_segments;
			m_dhu.parameters.qmu_mass_flux_segments=cell(size(m_dhu.parameters.qmu_mass_flux_num_segments));
			ipt=1;
			for i=1:length(m_dhu.parameters.qmu_mass_flux_num_segments)
				if m_dhu.parameters.qmu_mass_flux_num_segments(i)
					m_dhu.parameters.qmu_mass_flux_segments{i}=segments(ipt:ipt+m_dhu.parameters.qmu_mass_flux_num_segments(i)-1,:);
					ipt=ipt+m_dhu.parameters.qmu_mass_flux_num_segments(i);
				end
			end
			clear segments
		end

		if isnumeric(m_dhu.parameters.qmu_mass_flux_segments)
			response=MassFlux(m_dhu.elements,m_dhu.nodes,m_dhu.vertices,m_dhu.loads,m_dhu.materials,m_dhu.parameters,results.u_g);
		else
			save=m_dhu.parameters.qmu_mass_flux_segments;
			m_dhu.parameters.qmu_mass_flux_segments=m_dhu.parameters.qmu_mass_flux_segments{indx};
			response=MassFlux(m_dhu.elements,m_dhu.nodes,m_dhu.vertices,m_dhu.loads,m_dhu.materials,m_dhu.parameters,results.u_g);
			m_dhu.parameters.qmu_mass_flux_segments=save;
			clear save
		end

	elseif isSSA || isHO,

% for now, separate all segments from double array for parallel to make cells
		if (length(m_dh.parameters.qmu_mass_flux_num_segments) > 1)
			segments=m_dh.parameters.qmu_mass_flux_segments;
			m_dh.parameters.qmu_mass_flux_segments=cell(size(m_dh.parameters.qmu_mass_flux_num_segments));
			ipt=1;
			for i=1:length(m_dh.parameters.qmu_mass_flux_num_segments)
				if m_dh.parameters.qmu_mass_flux_num_segments(i)
					m_dh.parameters.qmu_mass_flux_segments{i}=segments(ipt:ipt+m_dh.parameters.qmu_mass_flux_num_segments(i)-1,:);
					ipt=ipt+m_dh.parameters.qmu_mass_flux_num_segments(i);
				end
			end
			clear segments
		end

		if isnumeric(m_dh.parameters.qmu_mass_flux_segments)
			response=MassFlux(m_dh.elements,m_dh.nodes,m_dh.vertices,m_dh.loads,m_dh.materials,m_dh.parameters,results.u_g);
		else
			save=m_dh.parameters.qmu_mass_flux_segments;
			m_dh.parameters.qmu_mass_flux_segments=m_dh.parameters.qmu_mass_flux_segments{indx};
			response=MassFlux(m_dh.elements,m_dh.nodes,m_dh.vertices,m_dh.loads,m_dh.materials,m_dh.parameters,results.u_g);
			m_dh.parameters.qmu_mass_flux_segments=save;
			clear save
		end

	elseif isFS,

% for now, separate all segments from double array for parallel to make cells
		if (length(m_ds.parameters.qmu_mass_flux_num_segments) > 1)
			segments=m_ds.parameters.qmu_mass_flux_segments;
			m_ds.parameters.qmu_mass_flux_segments=cell(size(m_ds.parameters.qmu_mass_flux_num_segments));
			ipt=1;
			for i=1:length(m_ds.parameters.qmu_mass_flux_num_segments)
				if m_ds.parameters.qmu_mass_flux_num_segments(i)
					m_ds.parameters.qmu_mass_flux_segments{i}=segments(ipt:ipt+m_ds.parameters.qmu_mass_flux_num_segments(i)-1,:);
					ipt=ipt+m_ds.parameters.qmu_mass_flux_num_segments(i);
				end
			end
			clear segments
		end

		if isnumeric(m_ds.parameters.qmu_mass_flux_segments)
			response=MassFlux(m_ds.elements,m_ds.nodes,m_ds.vertices,m_ds.loads,m_ds.materials,m_ds.parameters,results.u_g);
		else
			save=m_ds.parameters.qmu_mass_flux_segments;
			m_ds.parameters.qmu_mass_flux_segments=m_ds.parameters.qmu_mass_flux_segments{indx};
			response=MassFlux(m_ds.elements,m_ds.nodes,m_ds.vertices,m_ds.loads,m_ds.materials,m_ds.parameters,results.u_g);
			m_ds.parameters.qmu_mass_flux_segments=save;
			clear save
		end
	else
		error('qmuresponse error message: unsupported analysis type for mass_flux computation!');
	end
else
	error(['qmuresponse error message: unknown descriptor ' descriptor]);
end
