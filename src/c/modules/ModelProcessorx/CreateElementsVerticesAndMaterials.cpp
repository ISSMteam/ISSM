/*
 * CreateElementsNodesAndMaterialsStressbalanceHoriz.c:
 */

#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

#define MAXCONNECTIVITY 15

bool IsVertexInRank(int* vertices_ranks,int* vertices_proc_count,int vid,int rank){/*{{{*/

	/*See if node is already in partition*/
	for(int k=0;k<vertices_proc_count[vid];k++){
		if(vertices_ranks[MAXCONNECTIVITY*vid+k] == rank) return true;
	}

	return false;
}/*}}}*/
void AddVertexToRank(int* vertices_ranks,int* vertices_proc_count,int vid,int rank){/*{{{*/

	/*See if node is already in partition, return if this is the case*/
	if(IsVertexInRank(vertices_ranks,vertices_proc_count,vid,rank)) return;

	/*This rank has not been marked for this node just yet so go ahead and add it*/
	if(vertices_proc_count[vid]==MAXCONNECTIVITY) _error_("This vertex is connected to more than "<<MAXCONNECTIVITY<<" partition. Either reduce the number of processors, or increase MAXCONNECTIVITY");
	vertices_ranks[MAXCONNECTIVITY*vid+vertices_proc_count[vid]] = rank;
	vertices_proc_count[vid]++;
}/*}}}*/

void CreateElements(Elements* elements,IoModel* iomodel,const int nummodels){/*{{{*/

	/*Intermediary*/
	bool control_analysis;
	bool adolc_analysis;

	/*Fetch parameters: */
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	iomodel->FindConstant(&adolc_analysis,"md.autodiff.isautodiff");

	/*Did we already create the elements? : */
	_assert_(elements->Size()==0);

	/*Create elements*/
	if(control_analysis && !adolc_analysis)iomodel->FetchData(2,"md.inversion.min_parameters","md.inversion.max_parameters");
	if(iomodel->domaintype==Domain2DverticalEnum || iomodel->domaindim==3)  iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");

	int count = 0;
	switch(iomodel->meshelementtype){
		case TriaEnum:
			for(int i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					elements->AddObject(new Tria(i+1,i,count,iomodel,nummodels));
					count++;
				}
			}
			break;
		case TetraEnum:
			for(int i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					elements->AddObject(new Tetra(i+1,i,count,iomodel,nummodels));
					count++;
				}
			}
			break;
		case PentaEnum:
			iomodel->FetchData(2,"md.mesh.upperelements","md.mesh.lowerelements");
			for(int i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					elements->AddObject(new Penta(i+1,i,count,iomodel,nummodels));
					count++;
				}
			}
			break;
		default:
			_error_("Mesh not supported yet");
	}

	/*Free data: */
	iomodel->DeleteData(6,"md.mesh.upperelements","md.mesh.lowerelements","md.inversion.min_parameters","md.inversion.max_parameters","md.mesh.vertexonbase","md.mesh.vertexonsurface");

}/*}}}*/
void CreateMaterials(Elements* elements,Inputs* inputs,Materials* materials,IoModel* iomodel,const int nummodels){/*{{{*/

	/*Intermediary*/
	int  i;
	int  nnat,dummy;
	int* nature=NULL;

	/*Fetch parameters: */
	int materials_type;
	iomodel->FindConstant(&materials_type,"md.materials.type");

	/*Did we already create the materials? : */
	_assert_(materials->Size()==0);

	/*Create materials*/
	switch(materials_type){
		case MaticeEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			for (i=0;i<iomodel->numberofelements;i++) if(iomodel->my_elements[i]) materials->AddObject(new Matice(i+1,i,iomodel));
			switch(iomodel->domaindim){
				case 2:
					inputs->DuplicateInput(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
					break;
				case 3:
					break;
				default:
					_error_("Mesh not supported yet");
			}
			break;
		case MatenhancediceEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_E",MaterialsRheologyEEnum);
			for (i=0;i<iomodel->numberofelements;i++) if(iomodel->my_elements[i]) materials->AddObject(new Matice(i+1,i,iomodel));
			switch(iomodel->domaindim){
				case 2:
					inputs->DuplicateInput(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
					inputs->DuplicateInput(MaterialsRheologyEEnum,MaterialsRheologyEbarEnum);
					break;
				case 3:
					break;
				default:
					_error_("Mesh not supported yet");
			}
			break;
		case MatdamageiceEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.damage.D",DamageDEnum);
			for (i=0;i<iomodel->numberofelements;i++) if(iomodel->my_elements[i]) materials->AddObject(new Matice(i+1,i,iomodel));
			switch(iomodel->domaindim){
				case 2:
					inputs->DuplicateInput(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
					inputs->DuplicateInput(DamageDEnum,DamageDbarEnum);
					break;
				case 3:
					break;
				default:
					_error_("Mesh not supported yet");
			}
			break;
		case MatestarEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_Ec",MaterialsRheologyEcEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_Es",MaterialsRheologyEsEnum);
			for(i=0;i<iomodel->numberofelements;i++) if(iomodel->my_elements[i]) materials->AddObject(new Matestar(i+1,i,iomodel));
			switch(iomodel->domaindim){
				case 2:
					inputs->DuplicateInput(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
					inputs->DuplicateInput(MaterialsRheologyEcEnum,MaterialsRheologyEcbarEnum);
					inputs->DuplicateInput(MaterialsRheologyEsEnum,MaterialsRheologyEsbarEnum);
					break;
				case 3:
					break;
				default:
					_error_("Mesh not supported yet");
			}
			break;
		case MaterialsEnum: 

			//we have several types of materials. Retrieve this info first: 
			iomodel->FetchData(&nature,&nnat,&dummy,"md.materials.nature");

			//make sure materials that are not tied to elements come last:  for now, only Matlitho qualifies.
			for(int i=0;i<nnat;i++){ 
				if (IoCodeToEnumNature(nature[i])==MatlithoEnum){
					int temp=nature[nnat-1];
					nature[nnat-1]=nature[i];
					nature[i]=temp;
				}
			}

			//go through list of materials, and create them: 
			for(int nat=0;nat<nnat;nat++){
				switch(IoCodeToEnumNature(nature[nat])){ 
					case MaticeEnum:{ /*{{{*/
							iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
							iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
							for (int k=0;k<iomodel->numberofelements;k++) if(iomodel->my_elements[k]) materials->AddObject(new Matice(k+1,k,MaticeEnum));
							switch(iomodel->domaindim){
								case 2:
									inputs->DuplicateInput(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
									break;
								case 3:
									break;
								default:
									_error_("Mesh not supported yet");
							}
						} /*}}}*/
						break;
					case MatlithoEnum: { /*{{{*/
							bool* issolid = NULL;
							int*  rheologymodel = NULL;
							iomodel->FetchData(&issolid, NULL, NULL, "md.materials.issolid");
							iomodel->FetchData(&rheologymodel, NULL, NULL, "md.materials.rheologymodel");
							iomodel->FetchData(11,"md.materials.radius","md.materials.viscosity","md.materials.lame_lambda","md.materials.lame_mu","md.materials.burgers_viscosity","md.materials.burgers_mu","md.materials.ebm_alpha","md.materials.ebm_delta","md.materials.ebm_taul","md.materials.ebm_tauh","md.materials.density");
							materials->AddObject(new Matlitho(iomodel->numberofelements+1, iomodel, issolid, rheologymodel));
							iomodel->DeleteData(11,"md.materials.radius","md.materials.viscosity","md.materials.lame_lambda","md.materials.lame_mu","md.materials.burgers_viscosity","md.materials.burgers_mu","md.materials.ebm_alpha","md.materials.ebm_delta","md.materials.ebm_taul","md.materials.ebm_tauh","md.materials.density");
							xDelete<bool>(issolid);
							xDelete<int>(rheologymodel);
						}
						/*}}}*/
						break;
					case MathydroEnum: {/*{{{*/
							/*If we don't have any materials pointed to by elements (meaning, if we are running only litho or hydro), 
							 * then we need to zero out the hmaterial pointers inside the elements dataset so that it won't error out 
							 * during configuration: */
							bool isice=false;
							for (int j=0;j<nnat;j++){
								if((IoCodeToEnumNature(nature[j])==MaticeEnum)||
										(IoCodeToEnumNature(nature[j])==MatenhancediceEnum)||
										(IoCodeToEnumNature(nature[j])==MatestarEnum)||
										(IoCodeToEnumNature(nature[j])==MatdamageiceEnum)){
									isice=true; break; }
							}
							if (!isice){
								/*go through elements, and zero the hmaterials pointers: */
								for(Object* & object : elements->objects){
									Element* element=xDynamicCast<Element*>(object);
									switch(element->ObjectEnum()){
										case TriaEnum: 
											{
												Tria* tria= xDynamicCast<Tria*>(element);
												if(tria->hmaterial)delete tria->hmaterial; 
												tria->hmaterial=NULL;
											}
											break;
										default: 
											_error_("Not implemented yet");
									}
								}
							}
						} /*}}}*/
						break;
					case MatenhancediceEnum: /*{{{*/
						iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
						iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
						iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_E",MaterialsRheologyEEnum);
						for(int j=0;j<iomodel->numberofelements;j++) if(iomodel->my_elements[j]) materials->AddObject(new Matice(j+1,j,MatenhancediceEnum));
						switch(iomodel->domaindim){
							case 2:
								inputs->DuplicateInput(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
								inputs->DuplicateInput(MaterialsRheologyEEnum,MaterialsRheologyEbarEnum);
								break;
							case 3:
								break;
							default:
								_error_("Mesh not supported yet");
						}
						/*}}}*/
						break;
					case MatdamageiceEnum: /*{{{*/
						iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
						iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
						iomodel->FetchDataToInput(inputs,elements,"md.damage.D",DamageDEnum);
						for (int j=0;j<iomodel->numberofelements;j++) if(iomodel->my_elements[j]) materials->AddObject(new Matice(j+1,j,MatdamageiceEnum));
						switch(iomodel->domaindim){
							case 2:
								inputs->DuplicateInput(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
								inputs->DuplicateInput(DamageDEnum,DamageDbarEnum);
								break;
							case 3:
								break;
							default:
								_error_("Mesh not supported yet");
						}
						/*}}}*/
						break;
					case MatestarEnum: /*{{{*/
						iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
						iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_Ec",MaterialsRheologyEcEnum);
						iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_Es",MaterialsRheologyEsEnum);
						for(int j=0;j<iomodel->numberofelements;j++) if(iomodel->my_elements[j]) materials->AddObject(new Matestar(j+1,j,iomodel));
						switch(iomodel->domaindim){
							case 2:
								inputs->DuplicateInput(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
								inputs->DuplicateInput(MaterialsRheologyEcEnum,MaterialsRheologyEcbarEnum);
								inputs->DuplicateInput(MaterialsRheologyEsEnum,MaterialsRheologyEsbarEnum);
								break;
							case 3:
								break;
							default:
								_error_("Mesh not supported yet");
						} 
						/*}}}*/
						break; 
					default:
						_error_("Materials nature type "<<EnumToStringx(IoCodeToEnumNature(nature[nat]))<<" not supported");
						break;
				}
			}
			//Free resources:
			xDelete<int>(nature);
			break;

		default:
			_error_("Materials "<<EnumToStringx(materials_type)<<" not supported");
	}

	/*Free data: */
	iomodel->DeleteData(3,"md.material.rheology_B","md.material.rheology_n","md.damage.D");
}/*}}}*/
void CreateVertices(Elements* elements,Vertices* vertices,IoModel* iomodel,int solution_type,bool isamr){/*{{{*/

	/*Get element partitionning*/
	int* epart = iomodel->epart;

	/*Determine element width*/
	int  elements_width;
	switch(iomodel->meshelementtype){
		case TriaEnum:  elements_width=3; break;
		case TetraEnum: elements_width=4; break;
		case PentaEnum: elements_width=6; break;
		default: _error_("mesh elements "<< EnumToStringx(iomodel->meshelementtype) <<" not supported yet");
	}

	/*Get my_rank:*/
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*create matrix that keeps track of all ranks that have vertex i, and initialize as -1 (Common to all CPUs)*/
	int* vertices_ranks = xNew<int>(MAXCONNECTIVITY*iomodel->numberofvertices);
	for(int i=0;i<MAXCONNECTIVITY*iomodel->numberofvertices;i++) vertices_ranks[i] = -1;

	/*For all vertices, count how many cpus hold vertex i (initialize with 0)*/
	int* vertices_proc_count = xNewZeroInit<int>(iomodel->numberofvertices);

	/*Go through all elements and mark all vertices for all partitions*/
	for(int i=0;i<iomodel->numberofelements;i++){
		for(int j=0;j<elements_width;j++){
			/*Get current vertex sid*/
			int vid = iomodel->elements[elements_width*i+j]-1;
			AddVertexToRank(vertices_ranks,vertices_proc_count,vid,epart[i]);
		}
	}

	/*Take care of penalties (only in non-AMR for now)*/
	if(!isamr){
		int numvertex_pairing;
		int *vertex_pairing = NULL;
		iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.stressbalance.vertex_pairing");
		for(int i=0;i<numvertex_pairing;i++){
			int id1 = vertex_pairing[2*i+0]-1;
			int id2 = vertex_pairing[2*i+1]-1;
			for(int e=0;e<num_procs;e++){
				if(IsVertexInRank(vertices_ranks,vertices_proc_count,id1,e)){
					AddVertexToRank(vertices_ranks,vertices_proc_count,id2,e);
				}
			}
		}
		xDelete<int>(vertex_pairing);
		iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.masstransport.vertex_pairing");
		for(int i=0;i<numvertex_pairing;i++){
			int id1 = vertex_pairing[2*i+0]-1;
			int id2 = vertex_pairing[2*i+1]-1;
			for(int e=0;e<num_procs;e++){
				if(IsVertexInRank(vertices_ranks,vertices_proc_count,id1,e)){
					AddVertexToRank(vertices_ranks,vertices_proc_count,id2,e);
				}
			}
		}
		xDelete<int>(vertex_pairing);
	}

	/*Create vector of size total numnodes, initialized with -1, that will keep track of local ids*/
	int  offset = 0;
	int* vertices_offsets  = xNew<int>(iomodel->numberofvertices);
	for(int i=0;i<iomodel->numberofvertices;i++){
		if(IsVertexInRank(vertices_ranks,vertices_proc_count,i,my_rank)){
			vertices_offsets[i] = offset++;
		}
		else{
			vertices_offsets[i] = -1;
		}
	}

	/*Now, Count how many clones we have with other partitions*/
	int*  common_send = xNew<int>(num_procs);
	int*  common_recv = xNew<int>(num_procs);
	int** common_send_ids = xNew<int*>(num_procs);
	int** common_recv_ids = xNew<int*>(num_procs);

	/*First step: allocate, Step 2: populate*/
	for(int step=0;step<2;step++){

		if(step==1){
			/*Allocate send and receive arrays of ids now*/
			for(int i=0;i<num_procs;i++){
				_assert_(common_send[i]>=0 && common_recv[i]>=0);
				common_send_ids[i] = xNew<int>(common_send[i]);
				common_recv_ids[i] = xNew<int>(common_recv[i]);
			}
		}

		/*Re/Initialize counters to 0*/
		for(int i=0;i<num_procs;i++){
			common_recv[i]=0;
			common_send[i]=0;
		}

		/*Go through table and find clones/masters etc*/
		for(int i=0;i<iomodel->numberofvertices;i++){

			/*If we did not find this vertex in our current partition, go to next vertex*/
			if(vertices_offsets[i] == -1) continue;

			/*Find in what column this rank belongs*/
			int col = -1;
			for(int j=0;j<MAXCONNECTIVITY;j++){
				if(vertices_ranks[MAXCONNECTIVITY*i+j] == my_rank){
					col = j;
					break;
				}
			}
			_assert_(col!=-1);

			/*If col==0, it is either not on boundary, or a master*/
			if(col==0){
				/*1. is this vertex on the boundary? Skip if not*/
				if(vertices_ranks[MAXCONNECTIVITY*i+col+1]==-1){
					continue;
				}
				else{
					for(int j=1;j<vertices_proc_count[i];j++){
						_assert_(vertices_ranks[MAXCONNECTIVITY*i+j]>=0);
						int rank = vertices_ranks[MAXCONNECTIVITY*i+j];
						if(step==1){
							common_send_ids[rank][common_send[rank]] = vertices_offsets[i];
						}
						common_send[rank]++;
					}
				}
			}
			else{
				/*3. It is a slave, record that we need to receive for this cpu*/
				int rank = vertices_ranks[MAXCONNECTIVITY*i+0];
				if(step==1){
					common_recv_ids[rank][common_recv[rank]] = vertices_offsets[i];
				}
				common_recv[rank]++;
			}
		}
	}

	/*Create Vertices, depending on the constructor type: */
	if(solution_type!=LoveSolutionEnum) CreateNumberNodeToElementConnectivity(iomodel);
	if(!isamr){
		int isoceancoupling;
		iomodel->FindConstant(&isoceancoupling,"md.transient.isoceancoupling");

		//iomodel->FetchData(6,"md.mesh.x","md.mesh.y","md.mesh.z","md.geometry.base","md.geometry.thickness","md.mask.ice_levelset");
		iomodel->FetchData(5,"md.mesh.x","md.mesh.y","md.mesh.z","md.geometry.base","md.geometry.thickness");
		if (iomodel->domaintype == Domain3DsurfaceEnum) iomodel->FetchData(3,"md.mesh.lat","md.mesh.long","md.mesh.r");
		if (isoceancoupling) iomodel->FetchData(2,"md.mesh.lat","md.mesh.long");

		for(int i=0;i<iomodel->numberofvertices;i++){
			if(vertices_offsets[i]!=-1){
				bool isclone = (vertices_ranks[MAXCONNECTIVITY*i+0]!=my_rank);
				vertices->AddObject(new Vertex(i+1,i,isclone,iomodel,isamr));
			}
		}

		/*Free data: */
		//iomodel->DeleteData(6,"md.mesh.x","md.mesh.y","md.mesh.z","md.geometry.base","md.geometry.thickness","md.mask.ice_levelset");
		iomodel->DeleteData(5,"md.mesh.x","md.mesh.y","md.mesh.z","md.geometry.base","md.geometry.thickness");
		if (iomodel->domaintype == Domain3DsurfaceEnum) iomodel->DeleteData(3,"md.mesh.lat","md.mesh.long","md.mesh.r");
		if (isoceancoupling) iomodel->DeleteData(2,"md.mesh.lat","md.mesh.long");
	}
	else{
		for(int i=0;i<iomodel->numberofvertices;i++){
			if(vertices_offsets[i]!=-1){
				bool isclone = (vertices_ranks[MAXCONNECTIVITY*i+0]!=my_rank);
				vertices->AddObject(new Vertex(i+1,i,isclone,iomodel,isamr));
			}
		}
	}
	xDelete<int>(vertices_offsets);

	/*Final step, create my_vertices*/
	_assert_(!iomodel->my_vertices);
	iomodel->my_vertices = xNew<bool>(iomodel->numberofvertices);
	for(int i=0;i<iomodel->numberofvertices;i++){
		if(IsVertexInRank(vertices_ranks,vertices_proc_count,i,my_rank)){
			iomodel->my_vertices[i] = true;
		}
		else{
			iomodel->my_vertices[i] = false;
		}
	}

	/*Free data: */
	xDelete<int>(vertices_ranks);
	xDelete<int>(vertices_proc_count);

	/*Assign communicators*/
	vertices->common_send=common_send;
	vertices->common_recv=common_recv;
	vertices->common_send_ids=common_send_ids;
	vertices->common_recv_ids=common_recv_ids;

	/*Finalize Initialization*/
	vertices->Finalize(iomodel);
}/*}}}*/
