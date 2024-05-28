/*!\file NodesDofx
 * \brief: establish degrees of freedom for all nodes
 */

#include "./NodesDofx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void NodesDofx(Nodes* nodes, Parameters* parameters){

	/*Do we have any nodes for this analysis type? :*/
	if(!nodes->NumberOfNodes()) return;

	/*Do we really need to update dof indexings*/
	if(!nodes->RequiresDofReindexing()) return;

	if(VerboseModule()) _printf0_("   Renumbering degrees of freedom\n");

	/*Go through all nodes, and build degree of freedom lists. Each node gets a fixed number of dofs. When 
	 *a  node has already been distributed dofs on one cpu, all other cpus with the same node cannot distribute it 
	 *anymore. Use clone field to be sure of that: */
	nodes->DistributeDofs(GsetEnum);
	nodes->DistributeDofs(FsetEnum);
	nodes->DistributeDofs(SsetEnum);
}
