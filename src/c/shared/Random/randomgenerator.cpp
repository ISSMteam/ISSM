/*!\file: randomgenerator
 * \brief random number generating functions
 */

#include <iostream>
#include "./randomgenerator.h"
#include <cmath>
#include <chrono>

#undef M_PI
#define M_PI 3.141592653589793238462643

namespace rnd{

  /* Linear congruential engine */

  linear_congruential_engine::linear_congruential_engine(){/*{{{*/

      pseed = new unsigned int;
		*pseed = std::chrono::steady_clock::now().time_since_epoch()/std::chrono::milliseconds(1);
      a = 1103515245;		       // BSD Formula
      c = 12345;					     // BSD Formula
      m = 2147483648;			     // BSD Formula
			return;
	}/*}}}*/
  linear_congruential_engine::linear_congruential_engine(unsigned int _a, unsigned int _b, unsigned int _m){/*{{{*/

      pseed = new unsigned int;
		*pseed = std::chrono::steady_clock::now().time_since_epoch()/std::chrono::milliseconds(1);
      a = _a;
      c = _b;
      m = _m;
			return;
	}/*}}}*/
	linear_congruential_engine::~linear_congruential_engine(){}

  unsigned int linear_congruential_engine::get_m() { return m; }

  void linear_congruential_engine::seed( int s ) {
      if(s<0)
        *pseed = std::chrono::steady_clock::now().time_since_epoch()/std::chrono::milliseconds(1);
      else
        *pseed = (unsigned) s;
  }

  unsigned int linear_congruential_engine::generator(){
    *pseed = ( a * *pseed + c ) % m ;
    return *pseed;
  }

  void linear_congruential_engine::free_resources(){
	  delete pseed;
	  return;
  }

  /* Uniform distribution */

	uniform_distribution::uniform_distribution(){/*{{{*/
			a = 0.0;
			b = 1.0;
			return;
	}
	/*}}}*/

	uniform_distribution::uniform_distribution(double _a,double _b){/*{{{*/
			a = _a;
			b = _b;
			return;
	}
	/*}}}*/

	uniform_distribution::~uniform_distribution(){}

	double uniform_distribution::generator(rnd::linear_congruential_engine random_engine) {
			return (b-a)*(double) random_engine.generator()/ random_engine.get_m() + a;
	}

  /* Normal distribution */

	normal_distribution::normal_distribution(){/*{{{*/
			mean   = 0;
			sdev  = 1.0;
			return;
	}
	/*}}}*/

	normal_distribution::normal_distribution(double m,double s){/*{{{*/
			mean   = m;
			sdev  = s;
			return;
	}
	/*}}}*/

	normal_distribution::~normal_distribution(){}

	double normal_distribution::generator(rnd::linear_congruential_engine random_engine){/*{{{*/

			rnd::uniform_distribution	distribution;

			double u1 = distribution.generator(random_engine);
			double u2 = distribution.generator(random_engine);

			double R = sqrt(-2*log(u1));
			double theta = 2*M_PI*u2;

			return mean + sdev * (R*cos(theta));

	}
	/*}}}*/

  /* Log-Normal distribution */

  lognormal_distribution::lognormal_distribution(){/*{{{*/
      logmean   = 0;
      logsdev  = 1.0;
      return;
  }
  /*}}}*/

  lognormal_distribution::lognormal_distribution(double m,double s){/*{{{*/
      logmean   = m;
      logsdev  = s;
      return;
  }
  /*}}}*/

  lognormal_distribution::~lognormal_distribution(){}

  double lognormal_distribution::generator(rnd::linear_congruential_engine random_engine){/*{{{*/

      rnd::normal_distribution	distribution(logmean,logsdev);

      return exp(distribution.generator(random_engine));

  }
  /*}}}*/

  /* Chi-squared distribution */

  chi_squared_distribution::chi_squared_distribution(){/*{{{*/
      k  = 1;
      return;
  }
  /*}}}*/

  chi_squared_distribution::chi_squared_distribution(unsigned int dof){/*{{{*/
      k   = dof;
      return;
  }
  /*}}}*/

  chi_squared_distribution::~chi_squared_distribution(){}

  double chi_squared_distribution::generator(rnd::linear_congruential_engine random_engine){/*{{{*/

      rnd::normal_distribution	distribution;

      double rand = 0;

      for(int i=0;i<k;i++)
        rand = rand + pow(distribution.generator(random_engine),2);

      return rand;

  }
  /*}}}*/


  /* Exponential distribution */

  exponential_distribution::exponential_distribution(){/*{{{*/
      lambda  = 1.0;
      return;
  }
  /*}}}*/

  exponential_distribution::exponential_distribution(double scale){/*{{{*/
      lambda   = scale;
      return;
  }
  /*}}}*/

  exponential_distribution::~exponential_distribution(){}

  double exponential_distribution::generator(rnd::linear_congruential_engine random_engine){/*{{{*/

      rnd::uniform_distribution	distribution;

      return -1.0/lambda*log(1.0-distribution.generator(random_engine));

  }
  /*}}}*/

  /* Student-t distribution */

  student_t_distribution::student_t_distribution(){/*{{{*/
      k  = 1;
      return;
  }
  /*}}}*/

  student_t_distribution::student_t_distribution(unsigned int dof){/*{{{*/
      k   = dof;
      return;
  }
  /*}}}*/

  student_t_distribution::~student_t_distribution(){}

  double student_t_distribution::generator(rnd::linear_congruential_engine random_engine){/*{{{*/

      rnd::normal_distribution	normal_distribution;
      rnd::chi_squared_distribution	chi_squared_distribution(k);

      return normal_distribution.generator(random_engine)/sqrt(chi_squared_distribution.generator(random_engine)/k);

  }
  /*}}}*/

}
