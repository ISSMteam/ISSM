/*!\file BarystaticContributions.h
 * \brief: header file for barystatic contribution object
 */

#ifndef _BARYSTATICCONTRIBUTIONS_H_
#define _BARYSTATICCONTRIBUTIONS_H_

/*Headers:*/
class IoModel;
class Parameters;
class Results;
template <class doubletype> class Vector;
#include "../shared/shared.h"

class BarystaticContributions {

	public: 

		Vector<IssmDouble>* ice;  //contributions to every ice partition (size nice x 1)
		Vector<IssmDouble>* cumice;  //cumulated contributions to every ice partition
		int                 nice; //number of ice partitions 
		IssmDouble*         pice; //ice partition (nel)

		Vector<IssmDouble>* hydro;  //contributions to every hydro partition (size nhydro x 1)
		Vector<IssmDouble>* cumhydro;  //cumulated contributions to every hydro partition
		int                 nhydro; //number of hydro partitions 
		IssmDouble*         phydro; //hydro partition (nel)

		Vector<IssmDouble>* ocean;  //contributions to every ocean partition (size nocean x 1)
		Vector<IssmDouble>* cumocean;  //cumulated contributions to every ocean partition
		int                 nocean; //number of ocean partitions 
		IssmDouble*         pocean; //ocean partition (nel)

		/*BarystaticContributions constructors, destructors :*/
		BarystaticContributions(IoModel* iomodel );
		~BarystaticContributions();

		/*routines:*/
		IssmDouble Total();
		IssmDouble CumTotal();
		void Cumulate(Parameters* parameters);
		void Save(Results* results, Parameters* parameters, IssmDouble oceanarea);
		void Set(int eid, IssmDouble icevalue, IssmDouble hydrovalue, IssmDouble oceanvalue);
		void Reset();

};
#endif  /* _BARYSTATICCONTRIBUTIONS_H_ */
