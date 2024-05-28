/*!\file Segment.h
 * \brief: header file for node object
 */

#ifndef _SEGMENT_H_
#define _SEGMENT_H_

/*Headers:*/
/*{{{*/
#include "../datastructures/datastructures.h"
#include "../shared/Numerics/constants.h"
/*}}}*/

template <class doubletype> 
class Segment: public Object{

	public:
		int        eid;
		doubletype x1;
		doubletype y1;
		doubletype x2;
		doubletype y2;

		/*Segment constructors, destructors :*/
		Segment(){/*{{{*/
			this->eid = UNDEF;
			this->x1  = UNDEF;
			this->y1  = UNDEF;
			this->x2  = UNDEF;
			this->y2  = UNDEF;
		}
		/*}}}*/
		Segment(int segment_eid, doubletype segment_x1,doubletype segment_y1,doubletype segment_x2, doubletype segment_y2){/*{{{*/

			this->eid = segment_eid;
			this->x1  = segment_x1;
			this->y1  = segment_y1;
			this->x2  = segment_x2;
			this->y2  = segment_y2;

		}
		/*}}}*/
		~Segment(){/*{{{*/
		}
		/*}}}*/

		/*Object virtual functions definitions:*/
		Object* copy() {/*{{{*/
			return new Segment(this->eid,this->x1,this->y1,this->x2,this->y2);
		}
		/*}}}*/
		void DeepEcho(void){/*{{{*/
			this->Echo();
		}
		/*}}}*/
		void Echo(void){/*{{{*/

			_printf_("Segment:\n");
			_printf_("   eid: " << eid << "\n");
			_printf_("   node 1: " << this->x1 << "|" << this->y1 << "\n");
			_printf_("   node 2: " << this->x2 << "|" << this->y2 << "\n");

		}
		/*}}}*/
		int    Id(void){ return eid; }/*{{{*/
		/*}}}*/
		void Marshall(MarshallHandle* marshallhandle){/*{{{*/
			_error_("not implemented yet!"); 
		} 
		/*}}}*/
		int ObjectEnum(void){/*{{{*/

			return SegmentEnum;

		}
		/*}}}*/

};

#endif  /* _SEGMENT_H_ */
