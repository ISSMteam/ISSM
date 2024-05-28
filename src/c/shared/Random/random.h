/*!\file: random.h
 * \brief prototypes for random.h
 */ 

#ifndef _RANDOM_H_
#define _RANDOM_H_

void univariateNormal(IssmDouble* prand, IssmDouble mean, IssmDouble sdev, int seed=-1);
void multivariateNormal(IssmDouble** prand, int dim, IssmDouble mean, IssmDouble* covariancematrix, int seed=-1);
void multivariateNormal(IssmDouble** prand, int dim, IssmDouble* mean, IssmDouble* covariancematrix, int seed=-1);

#endif //ifndef _RANDOM_H_
