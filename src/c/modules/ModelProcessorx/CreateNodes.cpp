/*
 * CreateNodes.c:
 */

#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

#define MAXCONNECTIVITY 15

bool IsNodeInRank(int* nodes_ranks,int* nodes_proc_count,int nid,int rank){/*{{{*/

	/*See if node is already in partition*/
	for(int k=0;k<nodes_proc_count[nid];k++){
		if(nodes_ranks[MAXCONNECTIVITY*nid+k] == rank) return true;
	}

	return false;
}/*}}}*/
void AddNodeToRank(int* nodes_ranks,int* nodes_proc_count,int nid,int rank){/*{{{*/

	/*See if node is already in partition, return if this is the case*/
	if(IsNodeInRank(nodes_ranks,nodes_proc_count,nid,rank)) return;

	/*This rank has not been marked for this node just yet so go ahead and add it*/
	if(nodes_proc_count[nid]==MAXCONNECTIVITY) _error_("This node is connected to more than "<<MAXCONNECTIVITY<<" partition. Either reduce the number of processors, or increase MAXCONNECTIVITY");
	nodes_ranks[MAXCONNECTIVITY*nid+nodes_proc_count[nid]] = rank;
	nodes_proc_count[nid]++;
}/*}}}*/

void CreateNodes(Nodes* nodes, IoModel* iomodel,int analysis,int finite_element,bool isamr,int approximation,int* approximations){

	/*Intermediaries*/
	int        numnodes;
	int        element_numnodes;
	int        element_node_ids[40] = {0};

	/*Get partitioning variables*/
	int  my_rank   = IssmComm::GetRank();
	int  num_procs = IssmComm::GetSize();
	int*     epart = iomodel->epart;

	/*Determine how many nodes we have in total (which depends on the type of finite element)*/
	/*{{{*/
	if(iomodel->meshelementtype==TriaEnum){
		switch(finite_element){
			case P0DGEnum:
				numnodes = iomodel->numberofelements;
				break;
			case P1Enum:
				numnodes = iomodel->numberofvertices;
				break;
			case P1DGEnum:
				numnodes = 3*iomodel->numberofelements;
				break;
			case P1bubbleEnum: case P1bubblecondensedEnum:
				numnodes = iomodel->numberofvertices+iomodel->numberofelements;
				break;
			case P2Enum:
				EdgesPartitioning(iomodel);
				numnodes = iomodel->numberofvertices+iomodel->numberofedges;
				break;
			case P1P1Enum: case P1P1GLSEnum: 
				/*P1 velocity + P1 pressure element with GLS stabilization*/
				numnodes = (iomodel->numberofvertices) + iomodel->numberofvertices;
				break;
			case MINIEnum: case MINIcondensedEnum:
				/*P1+ velocity (bubble statically condensed), P1 pressure*/
				numnodes = (iomodel->numberofvertices+iomodel->numberofelements) + (iomodel->numberofvertices);
				break;
			case TaylorHoodEnum: case XTaylorHoodEnum:
				/*P2 velocity, P1 pressure*/
				EdgesPartitioning(iomodel);
				numnodes = (iomodel->numberofvertices+iomodel->numberofedges) + (iomodel->numberofvertices);
				break;
			case LATaylorHoodEnum:
				/*P2 velocity*/
				EdgesPartitioning(iomodel);
				numnodes = (iomodel->numberofvertices+iomodel->numberofedges);
				break;
			case CrouzeixRaviartEnum:
				/*P2b velocity, P1 DG pressure*/
				EdgesPartitioning(iomodel);
				numnodes = (iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements) + (3*iomodel->numberofelements);
				break;
			case LACrouzeixRaviartEnum:
				/*P2b velocity*/
				EdgesPartitioning(iomodel);
				numnodes = (iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements);
				break;
			default:
				_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
		}
	}
	else if(iomodel->meshelementtype==PentaEnum){
		switch(finite_element){
			case P1Enum:
				numnodes = iomodel->numberofvertices;
				break;
			case P1bubbleEnum: case P1bubblecondensedEnum:
				numnodes = iomodel->numberofvertices+iomodel->numberofelements;
				break;
			case P1xP2Enum:
				EdgesPartitioning(iomodel);
				numnodes = iomodel->numberofvertices+iomodel->numberofverticaledges;
				break;
			case P1xP3Enum:
				EdgesPartitioning(iomodel);
				numnodes = iomodel->numberofvertices+2*iomodel->numberofverticaledges;
				break;
			case P1xP4Enum:
				EdgesPartitioning(iomodel);
				numnodes = iomodel->numberofvertices+3*iomodel->numberofverticaledges;
				break;
			case P2xP1Enum:
				EdgesPartitioning(iomodel);
				numnodes = iomodel->numberofvertices+iomodel->numberofhorizontaledges;
				break;
			case P2Enum:
				EdgesPartitioning(iomodel);
				FacesPartitioning(iomodel);
				numnodes = iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces;
				break;
			case P2xP4Enum:
				EdgesPartitioning(iomodel);
				FacesPartitioning(iomodel);
				numnodes = iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces;
				break;
			case P1P1Enum: case P1P1GLSEnum: 
				/*P1 velocity + P1 pressure element with GLS stabilization*/
				numnodes = (iomodel->numberofvertices) + iomodel->numberofvertices;
				break;
			case MINIEnum: case MINIcondensedEnum:
				/*P1+ velocity (bubble statically condensed), P1 pressure*/
				numnodes = (iomodel->numberofvertices+iomodel->numberofelements) + (iomodel->numberofvertices);
				break;
			case TaylorHoodEnum: case XTaylorHoodEnum:
				/*P2 velocity, P1 pressure*/
				EdgesPartitioning(iomodel);
				FacesPartitioning(iomodel);
				numnodes = (iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces) + (iomodel->numberofvertices);
				break;
			case OneLayerP4zEnum:
				/*P2xP4 velocity, P1 pressure*/
				EdgesPartitioning(iomodel);
				FacesPartitioning(iomodel);
				numnodes = (iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces) + (iomodel->numberofvertices);
				break;
			default:
				_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
		}
	}
	else{
		_error_("mesh elements "<< EnumToStringx(iomodel->meshelementtype) <<" not supported yet");
	}
	/*}}}*/

	/*create matrix that keeps track of all ranks that have node i, and initialize as -1 (Common to all CPUs)*/
	int* nodes_ranks = xNew<int>(MAXCONNECTIVITY*numnodes);
	for(int i=0;i<MAXCONNECTIVITY*numnodes;i++) nodes_ranks[i] = -1;

	/*For all nodes, count how many cpus have node i (initialize with 0)*/
	int* nodes_proc_count = xNewZeroInit<int>(numnodes);

	/*Create vector of approximation per node (used for FS: vel or pressure)*/
	int* nodes_approx = xNew<int>(numnodes);
	if(approximations){
		for(int i=0;i<numnodes;i++) nodes_approx[i] = approximations[i];
	}
	else{
		for(int i=0;i<numnodes;i++) nodes_approx[i] = approximation;
	}

	/*Go through all elements and mark all vertices for all partitions*/
	/*{{{*/
	for(int i=0;i<iomodel->numberofelements;i++){

		/*Define nodes sids for each element*/
		if(iomodel->meshelementtype==TriaEnum){
			switch(finite_element){
				case P0DGEnum:
					element_numnodes=1;
					element_node_ids[0]=i;
					break;
				case P1Enum:
					element_numnodes=3;
					element_node_ids[0]=iomodel->elements[3*i+0]-1;
					element_node_ids[1]=iomodel->elements[3*i+1]-1;
					element_node_ids[2]=iomodel->elements[3*i+2]-1;
					break;
				case P1bubbleEnum: case P1bubblecondensedEnum:
					element_numnodes=4;
					element_node_ids[0]=iomodel->elements[3*i+0]-1;
					element_node_ids[1]=iomodel->elements[3*i+1]-1;
					element_node_ids[2]=iomodel->elements[3*i+2]-1;
					element_node_ids[3]=iomodel->numberofvertices+i;
					break;
				case P1DGEnum:
					element_numnodes=3;
					element_node_ids[0]=3*i+0;
					element_node_ids[1]=3*i+1;
					element_node_ids[2]=3*i+2;
					break;
				case P2Enum:
					element_numnodes = 6;
					element_node_ids[0]=iomodel->elements[3*i+0]-1;
					element_node_ids[1]=iomodel->elements[3*i+1]-1;
					element_node_ids[2]=iomodel->elements[3*i+2]-1;
					element_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+0];
					element_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+1];
					element_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+2];
					break;
				case P1P1Enum: case P1P1GLSEnum:
					element_numnodes = 3+3;
					element_node_ids[0]=iomodel->elements[3*i+0]-1;
					element_node_ids[1]=iomodel->elements[3*i+1]-1;
					element_node_ids[2]=iomodel->elements[3*i+2]-1;
					for(int n=0;n<3;n++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					element_node_ids[3]=iomodel->numberofvertices+iomodel->elements[3*i+0]-1;
					element_node_ids[4]=iomodel->numberofvertices+iomodel->elements[3*i+1]-1;
					element_node_ids[5]=iomodel->numberofvertices+iomodel->elements[3*i+2]-1;
					for(int n=3;n<6;n++) nodes_approx[element_node_ids[n]] = FSpressureEnum;
					break;
				case MINIEnum: case MINIcondensedEnum:
					element_numnodes = 4+3;
					element_node_ids[0]=iomodel->elements[3*i+0]-1;
					element_node_ids[1]=iomodel->elements[3*i+1]-1;
					element_node_ids[2]=iomodel->elements[3*i+2]-1;
					element_node_ids[3]=iomodel->numberofvertices+i;
					for(int n=0;n<4;n++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					element_node_ids[4]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[3*i+0]-1;
					element_node_ids[5]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[3*i+1]-1;
					element_node_ids[6]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[3*i+2]-1;
					for(int n=4;n<7;n++) nodes_approx[element_node_ids[n]] = FSpressureEnum;
					break;
				case TaylorHoodEnum: case XTaylorHoodEnum:
					element_numnodes = 6+3;
					element_node_ids[0]=iomodel->elements[3*i+0]-1;
					element_node_ids[1]=iomodel->elements[3*i+1]-1;
					element_node_ids[2]=iomodel->elements[3*i+2]-1;
					element_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+0];
					element_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+1];
					element_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+2];
					for(int n=0;n<6;n++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					element_node_ids[6]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[3*i+0]-1;
					element_node_ids[7]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[3*i+1]-1;
					element_node_ids[8]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[3*i+2]-1;
					for(int n=6;n<9;n++) nodes_approx[element_node_ids[n]] = FSpressureEnum;
					break;
				case LATaylorHoodEnum:
					element_numnodes = 6;
					element_node_ids[0]=iomodel->elements[3*i+0]-1;
					element_node_ids[1]=iomodel->elements[3*i+1]-1;
					element_node_ids[2]=iomodel->elements[3*i+2]-1;
					element_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+0];
					element_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+1];
					element_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+2];
					for(int n=0;n<6;n++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					break;
				case CrouzeixRaviartEnum:
					element_numnodes = 7+3;
					element_node_ids[0]=iomodel->elements[3*i+0]-1;
					element_node_ids[1]=iomodel->elements[3*i+1]-1;
					element_node_ids[2]=iomodel->elements[3*i+2]-1;
					element_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+0];
					element_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+1];
					element_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+2];
					element_node_ids[6]=iomodel->numberofvertices+iomodel->numberofedges+i;
					for(int n=0;n<7;n++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					element_node_ids[7]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements+3*i+0;
					element_node_ids[8]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements+3*i+1;
					element_node_ids[9]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements+3*i+2;
					for(int n=7;n<10;n++) nodes_approx[element_node_ids[n]] = FSpressureEnum;
					break;
				case LACrouzeixRaviartEnum:
					element_numnodes = 7;
					element_node_ids[0]=iomodel->elements[3*i+0]-1;
					element_node_ids[1]=iomodel->elements[3*i+1]-1;
					element_node_ids[2]=iomodel->elements[3*i+2]-1;
					element_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+0];
					element_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+1];
					element_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*i+2];
					element_node_ids[6]=iomodel->numberofvertices+iomodel->numberofedges+i;
					for(int n=0;n<7;n++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					break;
				default:
					_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
			}
		}
		else if(iomodel->meshelementtype==PentaEnum){
			switch(finite_element){
				case P1Enum:
					element_numnodes=6;
					element_node_ids[0]=iomodel->elements[6*i+0]-1;
					element_node_ids[1]=iomodel->elements[6*i+1]-1;
					element_node_ids[2]=iomodel->elements[6*i+2]-1;
					element_node_ids[3]=iomodel->elements[6*i+3]-1;
					element_node_ids[4]=iomodel->elements[6*i+4]-1;
					element_node_ids[5]=iomodel->elements[6*i+5]-1;
					break;
				case P1bubbleEnum: case P1bubblecondensedEnum:
					element_numnodes=7;
					element_node_ids[0]=iomodel->elements[6*i+0]-1;
					element_node_ids[1]=iomodel->elements[6*i+1]-1;
					element_node_ids[2]=iomodel->elements[6*i+2]-1;
					element_node_ids[3]=iomodel->elements[6*i+3]-1;
					element_node_ids[4]=iomodel->elements[6*i+4]-1;
					element_node_ids[5]=iomodel->elements[6*i+5]-1;
					element_node_ids[6]=iomodel->numberofvertices+i;
					break;
				case P1xP2Enum:
					element_numnodes=9;
					element_node_ids[0]=iomodel->elements[6*i+0]-1;
					element_node_ids[1]=iomodel->elements[6*i+1]-1;
					element_node_ids[2]=iomodel->elements[6*i+2]-1;
					element_node_ids[3]=iomodel->elements[6*i+3]-1;
					element_node_ids[4]=iomodel->elements[6*i+4]-1;
					element_node_ids[5]=iomodel->elements[6*i+5]-1;
					element_node_ids[6]=iomodel->numberofvertices+iomodel->elementtoverticaledgeconnectivity[3*i+0];
					element_node_ids[7]=iomodel->numberofvertices+iomodel->elementtoverticaledgeconnectivity[3*i+1];
					element_node_ids[8]=iomodel->numberofvertices+iomodel->elementtoverticaledgeconnectivity[3*i+2];
					break;
				case P1xP3Enum:
					element_numnodes=12;
					element_node_ids[ 0]=iomodel->elements[6*i+0]-1;
					element_node_ids[ 1]=iomodel->elements[6*i+1]-1;
					element_node_ids[ 2]=iomodel->elements[6*i+2]-1;
					element_node_ids[ 3]=iomodel->elements[6*i+3]-1;
					element_node_ids[ 4]=iomodel->elements[6*i+4]-1;
					element_node_ids[ 5]=iomodel->elements[6*i+5]-1;
					element_node_ids[ 6]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*i+0]+0;
					element_node_ids[ 7]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*i+1]+0;
					element_node_ids[ 8]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*i+2]+0;
					element_node_ids[ 9]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*i+0]+1;
					element_node_ids[10]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*i+1]+1;
					element_node_ids[11]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*i+2]+1;
					break;
				case P1xP4Enum:
					element_numnodes=15;
					element_node_ids[ 0]=iomodel->elements[6*i+0]-1; /*Vertex 1*/
					element_node_ids[ 1]=iomodel->elements[6*i+1]-1; /*Vertex 2*/
					element_node_ids[ 2]=iomodel->elements[6*i+2]-1; /*Vertex 3*/
					element_node_ids[ 3]=iomodel->elements[6*i+3]-1; /*Vertex 4*/
					element_node_ids[ 4]=iomodel->elements[6*i+4]-1; /*Vertex 5*/
					element_node_ids[ 5]=iomodel->elements[6*i+5]-1; /*Vertex 6*/
					element_node_ids[ 6]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*i+0]+0; /*mid vertical edge 1*/
					element_node_ids[ 7]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*i+1]+0; /*mid vertical edge 2*/
					element_node_ids[ 8]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*i+2]+0; /*mid vertical edge 3*/
					element_node_ids[ 9]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*i+0]+1; /* 1/4 vertical edge 1*/
					element_node_ids[10]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*i+1]+1; /* 1/4 vertical edge 2*/
					element_node_ids[11]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*i+2]+1; /* 1/4 vertical edge 3*/
					element_node_ids[12]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*i+0]+2; /* 3/4 vertical edge 1*/
					element_node_ids[13]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*i+1]+2; /* 3/4 vertical edge 2*/
					element_node_ids[14]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*i+2]+2; /* 3/4 vertical edge 3*/
					break;
				case P2xP1Enum:
					element_numnodes=12;
					element_node_ids[ 0]=iomodel->elements[6*i+0]-1;
					element_node_ids[ 1]=iomodel->elements[6*i+1]-1;
					element_node_ids[ 2]=iomodel->elements[6*i+2]-1;
					element_node_ids[ 3]=iomodel->elements[6*i+3]-1;
					element_node_ids[ 4]=iomodel->elements[6*i+4]-1;
					element_node_ids[ 5]=iomodel->elements[6*i+5]-1;
					element_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*i+0];
					element_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*i+1];
					element_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*i+2];
					element_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*i+3];
					element_node_ids[10]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*i+4];
					element_node_ids[11]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*i+5];
					break;
				case P2Enum:
					element_numnodes = 18;
					element_node_ids[ 0]=iomodel->elements[6*i+0]-1;
					element_node_ids[ 1]=iomodel->elements[6*i+1]-1;
					element_node_ids[ 2]=iomodel->elements[6*i+2]-1;
					element_node_ids[ 3]=iomodel->elements[6*i+3]-1;
					element_node_ids[ 4]=iomodel->elements[6*i+4]-1;
					element_node_ids[ 5]=iomodel->elements[6*i+5]-1;
					element_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+0];
					element_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+1];
					element_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+2];
					element_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+3];
					element_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+4];
					element_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+5];
					element_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+6];
					element_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+7];
					element_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+8];
					element_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*i+0];
					element_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*i+1];
					element_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*i+2];
					break;
				case P2xP4Enum:
					element_numnodes = 30;
					element_node_ids[ 0]=iomodel->elements[6*i+0]-1; /*Vertex 1*/
					element_node_ids[ 1]=iomodel->elements[6*i+1]-1; /*Vertex 2*/
					element_node_ids[ 2]=iomodel->elements[6*i+2]-1; /*Vertex 3*/
					element_node_ids[ 3]=iomodel->elements[6*i+3]-1; /*Vertex 4*/
					element_node_ids[ 4]=iomodel->elements[6*i+4]-1; /*Vertex 5*/
					element_node_ids[ 5]=iomodel->elements[6*i+5]-1; /*Vertex 6*/
					element_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+0]; /*mid vertical edge 1*/
					element_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+1]; /*mid vertical edge 2*/
					element_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+2]; /*mid vertical edge 3*/
					element_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+3]; /*mid basal edge 1*/
					element_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+4]; /*mid basal edge 2*/
					element_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+5]; /*mid basal edge 3*/
					element_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+6]; /*mid top edge 1*/
					element_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+7]; /*mid top edge 2*/
					element_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+8]; /*mid top edge 3*/
					element_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+0]; /* 1/4 vertical edge 1*/
					element_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+1]; /* 1/4 vertical edge 2*/
					element_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+2]; /* 1/4 vertical edge 3*/
					element_node_ids[18]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+0]+1; /* 3/4 vertical edge 1*/
					element_node_ids[19]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+1]+1; /* 3/4 vertical edge 2*/
					element_node_ids[20]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+2]+1; /* 3/4 vertical edge 3*/
					element_node_ids[21]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+0]+0; /* 1/4 vertical face 1*/
					element_node_ids[22]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+1]+0; /* 1/4 vertical face 2*/
					element_node_ids[23]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+2]+0; /* 1/4 vertical face 3*/
					element_node_ids[24]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+0]+1; /* 2/4 vertical face 1*/
					element_node_ids[25]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+1]+1; /* 2/4 vertical face 2*/
					element_node_ids[26]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+2]+1; /* 2/4 vertical face 3*/
					element_node_ids[27]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+0]+2; /* 3/4 vertical face 1*/
					element_node_ids[28]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+1]+2; /* 3/4 vertical face 2*/
					element_node_ids[29]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+2]+2; /* 3/4 vertical face 3*/
					break;
				case P1P1Enum: case P1P1GLSEnum:
					/*P1-P1 very similar to MINI, but no DOF on the element
					TODO: add if(!approximations) or not?*/
					element_numnodes = 6+6;
					element_node_ids[ 0]=iomodel->elements[6*i+0]-1;
					element_node_ids[ 1]=iomodel->elements[6*i+1]-1;
					element_node_ids[ 2]=iomodel->elements[6*i+2]-1;
					element_node_ids[ 3]=iomodel->elements[6*i+3]-1;
					element_node_ids[ 4]=iomodel->elements[6*i+4]-1;
					element_node_ids[ 5]=iomodel->elements[6*i+5]-1;
					if(!approximations) for(int n=0;n<6;n ++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					element_node_ids[ 6]=iomodel->numberofvertices+iomodel->elements[6*i+0]-1;
					element_node_ids[ 7]=iomodel->numberofvertices+iomodel->elements[6*i+1]-1;
					element_node_ids[ 8]=iomodel->numberofvertices+iomodel->elements[6*i+2]-1;
					element_node_ids[ 9]=iomodel->numberofvertices+iomodel->elements[6*i+3]-1;
					element_node_ids[10]=iomodel->numberofvertices+iomodel->elements[6*i+4]-1;
					element_node_ids[11]=iomodel->numberofvertices+iomodel->elements[6*i+5]-1;
					if(!approximations) for(int n=6;n<12;n++) nodes_approx[element_node_ids[n]] = FSpressureEnum;
					break;
				case MINIEnum: case MINIcondensedEnum:
					element_numnodes = 7+6;
					element_node_ids[ 0]=iomodel->elements[6*i+0]-1;
					element_node_ids[ 1]=iomodel->elements[6*i+1]-1;
					element_node_ids[ 2]=iomodel->elements[6*i+2]-1;
					element_node_ids[ 3]=iomodel->elements[6*i+3]-1;
					element_node_ids[ 4]=iomodel->elements[6*i+4]-1;
					element_node_ids[ 5]=iomodel->elements[6*i+5]-1;
					element_node_ids[ 6]=iomodel->numberofvertices+i;
					if(!approximations) for(int n=0;n<7;n++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					element_node_ids[ 7]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*i+0]-1;
					element_node_ids[ 8]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*i+1]-1;
					element_node_ids[ 9]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*i+2]-1;
					element_node_ids[10]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*i+3]-1;
					element_node_ids[11]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*i+4]-1;
					element_node_ids[12]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*i+5]-1;
					if(!approximations) for(int n=7;n<13;n++) nodes_approx[element_node_ids[n]] = FSpressureEnum;
					break;
				case TaylorHoodEnum: case XTaylorHoodEnum:
					element_numnodes = 18+6;
					element_node_ids[ 0]=iomodel->elements[6*i+0]-1;
					element_node_ids[ 1]=iomodel->elements[6*i+1]-1;
					element_node_ids[ 2]=iomodel->elements[6*i+2]-1;
					element_node_ids[ 3]=iomodel->elements[6*i+3]-1;
					element_node_ids[ 4]=iomodel->elements[6*i+4]-1;
					element_node_ids[ 5]=iomodel->elements[6*i+5]-1;
					element_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+0];
					element_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+1];
					element_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+2];
					element_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+3];
					element_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+4];
					element_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+5];
					element_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+6];
					element_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+7];
					element_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+8];
					element_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*i+0];
					element_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*i+1];
					element_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*i+2];
					for(int n=0;n<18;n++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					element_node_ids[18]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*i+0]-1;
					element_node_ids[19]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*i+1]-1;
					element_node_ids[20]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*i+2]-1;
					element_node_ids[21]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*i+3]-1;
					element_node_ids[22]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*i+4]-1;
					element_node_ids[23]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*i+5]-1;
					for(int n=18;n<24;n++) nodes_approx[element_node_ids[n]] = FSpressureEnum;
					break;
				case OneLayerP4zEnum:
					element_numnodes = 30+6;
					element_node_ids[ 0]=iomodel->elements[6*i+0]-1; /*Vertex 1*/
					element_node_ids[ 1]=iomodel->elements[6*i+1]-1; /*Vertex 2*/
					element_node_ids[ 2]=iomodel->elements[6*i+2]-1; /*Vertex 3*/
					element_node_ids[ 3]=iomodel->elements[6*i+3]-1; /*Vertex 4*/
					element_node_ids[ 4]=iomodel->elements[6*i+4]-1; /*Vertex 5*/
					element_node_ids[ 5]=iomodel->elements[6*i+5]-1; /*Vertex 6*/
					element_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+0]; /*mid vertical edge 1*/
					element_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+1]; /*mid vertical edge 2*/
					element_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+2]; /*mid vertical edge 3*/
					element_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+3]; /*mid basal edge 1*/
					element_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+4]; /*mid basal edge 2*/
					element_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+5]; /*mid basal edge 3*/
					element_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+6]; /*mid top edge 1*/
					element_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+7]; /*mid top edge 2*/
					element_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*i+8]; /*mid top edge 3*/
					element_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+0]; /* 1/4 vertical edge 1*/
					element_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+1]; /* 1/4 vertical edge 2*/
					element_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+2]; /* 1/4 vertical edge 3*/
					element_node_ids[18]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+0]+1; /* 3/4 vertical edge 1*/
					element_node_ids[19]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+1]+1; /* 3/4 vertical edge 2*/
					element_node_ids[20]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*i+2]+1; /* 3/4 vertical edge 3*/
					element_node_ids[21]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+0]+0; /* 1/4 vertical face 1*/
					element_node_ids[22]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+1]+0; /* 1/4 vertical face 2*/
					element_node_ids[23]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+2]+0; /* 1/4 vertical face 3*/
					element_node_ids[24]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+0]+1; /* 2/4 vertical face 1*/
					element_node_ids[25]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+1]+1; /* 2/4 vertical face 2*/
					element_node_ids[26]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+2]+1; /* 2/4 vertical face 3*/
					element_node_ids[27]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+0]+2; /* 3/4 vertical face 1*/
					element_node_ids[28]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+1]+2; /* 3/4 vertical face 2*/
					element_node_ids[29]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*i+2]+2; /* 3/4 vertical face 3*/
					for(int n=0;n<30;n++) nodes_approx[element_node_ids[n]] = FSvelocityEnum;
					element_node_ids[30]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*i+0]-1;
					element_node_ids[31]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*i+1]-1;
					element_node_ids[32]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*i+2]-1;
					element_node_ids[33]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*i+3]-1;
					element_node_ids[34]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*i+4]-1;
					element_node_ids[35]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*i+5]-1;
					for(int n=30;n<36;n++) nodes_approx[element_node_ids[n]] = FSpressureEnum;
					break;
				default:
					_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
			}
		}
		else{
			_error_("mesh elements "<< EnumToStringx(iomodel->meshelementtype) <<" not supported yet");
		}

		/*Add rank epart[i] for all nodes belonging to this element*/
		for(int j=0;j<element_numnodes;j++){
			int nid = element_node_ids[j]; _assert_(nid<numnodes);
			AddNodeToRank(nodes_ranks,nodes_proc_count,nid,epart[i]);
		}
	}
	/*}}}*/
	if((finite_element==P0DGEnum || finite_element==P1DGEnum) && analysis!=UzawaPressureAnalysisEnum){/*Special case for DG...{{{*/
		if(finite_element==P1DGEnum){
			int node_list[4];
			if(iomodel->domaintype!=Domain2DhorizontalEnum) _error_("not implemented yet");
			CreateEdges(iomodel);
			CreateFaces(iomodel);
			for(int i=0;i<iomodel->numberoffaces;i++){
				int e1=iomodel->faces[4*i+2]-1; //faces are [node1 node2 elem1 elem2]
				int e2=iomodel->faces[4*i+3]-1; //faces are [node1 node2 elem1 elem2]
				if(e2!=-2){
					if(epart[e1]!=epart[e2]){
						int i1=iomodel->faces[4*i+0];
						int i2=iomodel->faces[4*i+1];
						int pos=-1;
						for(int j=0;j<3;j++) if(iomodel->elements[3*e2+j]==i1) pos=j;
						if(     pos==0){ node_list[0] = e2*3+0; node_list[1] = e2*3+2;}
						else if(pos==1){ node_list[0] = e2*3+1; node_list[1] = e2*3+0;}
						else if(pos==2){ node_list[0] = e2*3+2; node_list[1] = e2*3+1;}
						else _error_("not supposed to happen");
						pos=-1;
						for(int j=0;j<3;j++) if(iomodel->elements[3*e1+j]==i1) pos=j;
						if(     pos==0){ node_list[2] = e1*3+0; node_list[3] = e1*3+1;}
						else if(pos==1){ node_list[2] = e1*3+1; node_list[3] = e1*3+2;}
						else if(pos==2){ node_list[2] = e1*3+2; node_list[3] = e1*3+0;}
						else _error_("not supposed to happen");
						for(int j=0;j<4;j++){
							int  nid = node_list[j];
							AddNodeToRank(nodes_ranks,nodes_proc_count,nid,epart[e1]);
							AddNodeToRank(nodes_ranks,nodes_proc_count,nid,epart[e2]);
						}
					}
				}
			}
		}
		else if(finite_element==P0DGEnum){
			int node_list[2];
			if(iomodel->domaintype!=Domain2DhorizontalEnum) _error_("not implemented yet");
			CreateEdges(iomodel);
			CreateFaces(iomodel);
			for(int i=0;i<iomodel->numberoffaces;i++){
				int e1=iomodel->faces[4*i+2]-1; //faces are [node1 node2 elem1 elem2]
				int e2=iomodel->faces[4*i+3]-1; //faces are [node1 node2 elem1 elem2]
				if(e2!=-2){
					if(epart[e1]!=epart[e2]){
						AddNodeToRank(nodes_ranks,nodes_proc_count,e2,epart[e1]);
						//AddNodeToRank(nodes_ranks,nodes_proc_count,e1,epart[e2]);
					}
				}
			}
		}
		else{
			_error_("not supported");
		}
	}/*}}}*/
	/*Vertex pairing for stressbalance{{{*/
	if(!isamr && (analysis==StressbalanceAnalysisEnum || analysis==StressbalanceVerticalAnalysisEnum)){
		int *vertex_pairing = NULL;
		int  numvertex_pairing;
		iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.stressbalance.vertex_pairing");
		if(numvertex_pairing>0){
			if(finite_element!=P1Enum) _error_("Periodic boundary conditions only implemented for linear elements");
		}
		for(int i=0;i<numvertex_pairing;i++){
			int nid1 = vertex_pairing[2*i+0]-1;
			int nid2 = vertex_pairing[2*i+1]-1;
			for(int j=0;j<nodes_proc_count[nid1];j++) AddNodeToRank(nodes_ranks,nodes_proc_count,nid2,nodes_ranks[MAXCONNECTIVITY*nid1+j]);
			for(int j=0;j<nodes_proc_count[nid2];j++) AddNodeToRank(nodes_ranks,nodes_proc_count,nid1,nodes_ranks[MAXCONNECTIVITY*nid2+j]);
		}
		xDelete<int>(vertex_pairing);
	}
	if(!isamr && (analysis==MasstransportAnalysisEnum
					|| analysis==FreeSurfaceBaseAnalysisEnum
					|| analysis==FreeSurfaceTopAnalysisEnum
					|| analysis==DebrisAnalysisEnum
					|| analysis==ThermalAnalysisEnum
					|| analysis==EnthalpyAnalysisEnum
					)){
		int *vertex_pairing = NULL;
		int  numvertex_pairing;
		iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.masstransport.vertex_pairing");
		_assert_(numvertex_pairing==0 || finite_element==P1Enum);
		for(int i=0;i<numvertex_pairing;i++){
			int nid1 = vertex_pairing[2*i+0]-1;
			int nid2 = vertex_pairing[2*i+1]-1;
			for(int j=0;j<nodes_proc_count[nid1];j++) AddNodeToRank(nodes_ranks,nodes_proc_count,nid2,nodes_ranks[MAXCONNECTIVITY*nid1+j]);
			for(int j=0;j<nodes_proc_count[nid2];j++) AddNodeToRank(nodes_ranks,nodes_proc_count,nid1,nodes_ranks[MAXCONNECTIVITY*nid2+j]);
		}
		xDelete<int>(vertex_pairing);
	}
	/*}}}*/

	/*Create vector of size total numnodes, initialized with -1, that will keep track of local ids*/
	int  offset = 0;
	int* nodes_offsets  = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++){
		if(IsNodeInRank(nodes_ranks,nodes_proc_count,i,my_rank)){
			nodes_offsets[i] = offset++;
		}
		else{
			nodes_offsets[i] = -1;
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
		for(int i=0;i<numnodes;i++){
			/*If we did not find this vertex in our current partition, go to next vertex*/
			if(nodes_offsets[i] == -1) continue;
			/*Find in what column this rank belongs*/
			int col = -1;
			for(int j=0;j<MAXCONNECTIVITY;j++){
				if(nodes_ranks[MAXCONNECTIVITY*i+j] == my_rank){
					col = j;
					break;
				}
			}
			_assert_(col!=-1);

			/*If col==0, it is either not on boundary, or a master*/
			if(col==0){
				/*1. is this vertex on the boundary? Skip if not*/
				if(nodes_ranks[MAXCONNECTIVITY*i+col+1]==-1){
					continue;
				}
				else{
					for(int j=1;j<nodes_proc_count[i];j++){
						_assert_(nodes_ranks[MAXCONNECTIVITY*i+j]>=0);
						int rank = nodes_ranks[MAXCONNECTIVITY*i+j];
						if(step==1){
							common_send_ids[rank][common_send[rank]] = nodes_offsets[i];
						}
						common_send[rank]++;
					}
				}
			}
			else{
				/*3. It is a slave, record that we need to receive for this cpu*/
				int rank = nodes_ranks[MAXCONNECTIVITY*i+0];
				if(step==1){
					common_recv_ids[rank][common_recv[rank]] = nodes_offsets[i];
				}
				common_recv[rank]++;
			}
		}
	}
	xDelete<int>(nodes_proc_count);

	/*Go ahead and create vertices now that we have all we need*/
	for(int i=0;i<numnodes;i++){
		if(nodes_offsets[i]!=-1){
			bool isclone = (nodes_ranks[MAXCONNECTIVITY*i+0]!=my_rank);
			int io_index = 0;
			if(i<iomodel->numberofvertices) io_index = i;
			Node* node=new Node(i+1,i,io_index,isclone,iomodel,analysis,nodes_approx[i],isamr);
			if(finite_element==MINIcondensedEnum || finite_element==P1bubblecondensedEnum){
				/*Bubble function is collapsed, needs to constrain it, maybe this is not the best place to do this, but that's life!*/
				if(i>=iomodel->numberofvertices && i<iomodel->numberofvertices+iomodel->numberofelements){
					node->HardDeactivate();
				}
			}
			nodes->AddObject(node);
		}
	}

	/*Free data: */
	xDelete<int>(nodes_approx);
	xDelete<int>(nodes_ranks);
	xDelete<int>(nodes_offsets);

	/*Assign communicators*/
	nodes->common_send=common_send;
	nodes->common_recv=common_recv;
	nodes->common_send_ids=common_send_ids;
	nodes->common_recv_ids=common_recv_ids;

	/*Finalize Initialization*/
	nodes->Finalize();
}
