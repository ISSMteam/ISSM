/* file: Interpolation .cpp
	Contains functions to interpolate from grid field to a given coordinate
 */
#include "./numerics.h"
#include "../Exceptions/exceptions.h"

/*triangleinterp{{{*/
double triangleinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y){
	/*split the rectangle in 2 triangle and
	 * use Lagrange P1 interpolation
	 *
	 *   +3----+2,3' Q12----Q22
	 *   |    /|     |    /|
	 *   |   / |     |   / |
	 *   |  /  |     |  /  |
	 *   | /   |     | /   |
	 *   |/    |     |/    |
	 *   1-----2'    Q11---Q21        */

	/*Intermediaries*/
	double area,area_1,area_2,area_3;

	/*Checks*/
	_assert_(x2>x1 && y2>y1);
	_assert_(x<=x2 && x>=x1 && y<=y2 && y>=y1);

	/*area of the rectangle*/
	area=(x2-x1)*(y2-y1);

	/*is it the upper left triangle?*/
	if ((x-x1)/(x2-x1)<(y-y1)/(y2-y1)){

		area_1=((y2-y)*(x2-x1))/area;
		area_2=((x-x1)*(y2-y1))/area;
		area_3=1-area_1-area_2;

		return area_1*Q11+area_2*Q22+area_3*Q12;
	}
	else {

		area_1=((y-y1)*(x2-x1))/area;
		area_2=((x2-x)*(y2-y1))/area;
		area_3=1-area_1-area_2;

		return area_1*Q22+area_2*Q11+area_3*Q21;
	}
}/*}}}*/
/*bilinearinterp{{{*/
IssmDouble bilinearinterp(IssmDouble* x_grid,IssmDouble* y_grid,IssmDouble* data,IssmDouble x,IssmDouble y,int m,int n,int Nx){
	IssmDouble x1, x2, y1, y2;
	IssmDouble Q11, Q12, Q21, Q22;
	x1=x_grid[n]; x2=x_grid[n+1];
	y1=y_grid[m]; y2=y_grid[m+1];

	Q11=data[m*Nx+n];
	Q12=data[(m+1)*Nx+n];
	Q21=data[m*Nx+n+1];
	Q22=data[(m+1)*Nx+n+1];
	return bilinearinterp<IssmDouble>(x1, x2, y1, y2, Q11, Q12, Q21, Q22, x, y);
}
/*}}}*/
/*nearestinterp{{{*/
double nearestinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y){
	/*Nearest neighbor interpolation*/

	/*    Q12             Q22
	 * y2 x--------x---------x
	 *    |        |         |
	 *    |        |  xP     |
	 * ym |--------+---------|
	 *    |        |         |
	 *    |        |         |
	 * y1 x--------x---------x Q21
	 *    x1       xm        x2
	 *
	 */
	/*Checks*/
	_assert_(x2>x1 && y2>y1);
	_assert_(x<=x2 && x>=x1 && y<=y2 && y>=y1);

	double xm=(x2-x1)/2;
	double ym=(y2-y1)/2;

	if (x<xm && y<ym) return Q11;
	if (x<xm && y>ym) return Q12;
	if (x>xm && y<ym) return Q21;
	else return Q22;
}
/*}}}*/
