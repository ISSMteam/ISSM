/*!\file: randomgenerator.h
 * \brief prototypes for randomgenerator.h
 */

#ifndef _RANDOMGENERATOR_H_
#define _RANDOMGENERATOR_H_

#undef M_PI
#define M_PI 3.141592653589793238462643

namespace rnd{

  class linear_congruential_engine
  {
    private:
      unsigned int a;
      unsigned int c;
      unsigned int m;
      unsigned int *pseed;

    public:

      /*constructors, destructors: */
      linear_congruential_engine();
      linear_congruential_engine(unsigned int _a, unsigned int _b, unsigned int _m);
      ~linear_congruential_engine();

      unsigned int get_m();
      void seed( int s );
      unsigned int generator();
		void free_resources();
  };

  class uniform_distribution
  {

    private:
      double a;  // lower bound of range
      double b;  // upper bound of range

    public:

      /*constructors, destructors: */
      uniform_distribution();
      uniform_distribution(double _a,double _b);
      ~uniform_distribution();

      double generator(rnd::linear_congruential_engine random_engine);

  };

  class normal_distribution
  {

    private:
      double mean;
      double sdev;

    public:

      /*constructors, destructors: */
      normal_distribution();
      normal_distribution(double m,double s);
      ~normal_distribution();

      double generator(rnd::linear_congruential_engine random_engine);

  };

  class lognormal_distribution
  {

    private:
      double logmean;
      double logsdev;

    public:

      /*constructors, destructors: */
      lognormal_distribution();
      lognormal_distribution(double m,double s);
      ~lognormal_distribution();

      double generator(rnd::linear_congruential_engine random_engine);

  };

  class chi_squared_distribution
  {

    private:
      unsigned int k;

    public:

      /*constructors, destructors: */
      chi_squared_distribution();
      chi_squared_distribution(unsigned int k);
      ~chi_squared_distribution();

      double generator(rnd::linear_congruential_engine random_engine);

  };

  class exponential_distribution
  {

    private:
      double lambda;

    public:

      /*constructors, destructors: */
      exponential_distribution();
      exponential_distribution(double scale);
      ~exponential_distribution();

      double generator(rnd::linear_congruential_engine random_engine);

  };

  class student_t_distribution
  {

    private:
      unsigned int k;

    public:

      /*constructors, destructors: */
      student_t_distribution();
      student_t_distribution(unsigned int k);
      ~student_t_distribution();

      double generator(rnd::linear_congruential_engine random_engine);

  };


}

#endif //* _RANDOMGENERATOR_H_ */
