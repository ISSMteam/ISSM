/*!\file Radar.h
 * \brief: header file for Radar object
 */

#ifndef _RADAR_H_
#define _RADAR_H_

/*Headers:*/
#include "./Definition.h"
#include "./FemModel.h"

class Radar: public Object, public Definition{

	public: 
		char* name;
		int	definitionenum;

		/*Radar constructors, destructors :*/
		Radar();
		Radar(char* in_name, int in_definitionenum);
		~Radar();

		/*Object virtual function resolutoin: */
		Object* copy();
		void DeepEcho(void);
		void Echo(void);
		int  Id(void);
		void Marshall(MarshallHandle* marshallhandle);
		int ObjectEnum(void);

		/*Definition virtual function resolutoin: */
		int DefinitionEnum();
		char* Name();
		IssmDouble Response(FemModel* femmodel);
		void ComputeRadarAttenuation(Element* element);
		void ComputeRadarPower(Element* element);
};
#endif  /* _RADAR_H_ */
