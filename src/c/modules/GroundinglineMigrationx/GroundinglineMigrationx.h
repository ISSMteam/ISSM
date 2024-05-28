/*!\file:  GroundinglineMigrationx.h
 * \brief header file for Grounding Line Migration
 */ 

#ifndef _GROUNDINGLINEMIGRATIONX_H
#define _GROUNDINGLINEMIGRATIONX_H

class Elements;
class Vertices;
class Nodes;
class Parameters;

/* local prototypes: */
void         GroundinglineMigrationx(Elements* elements,Nodes* nodes, Vertices* vertices,Loads* loads,Materials* materials, Parameters* parameters);
IssmDouble*  ContactFSLevelset(Elements* elements,Vertices* vertices);
IssmDouble*  PotentialUngrounding(Elements* elements,Vertices* vertices,Parameters* parameters);
IssmDouble*  PropagateFloatingiceToGroundedNeighbors(Elements* elements,Nodes* nodes,Vertices* vertices,Parameters* parameters,IssmDouble* vertices_potentially_ungrounding);
#endif  /* _GROUNDINGLINEMIGRATIONX_H */
