//	werner-ng/src/config.h
//	This file is a part of 'Werner NG' project.
/** \file config.h
**	\brief Declares very common type names and macros.
** 
** \author	Igor Krivenko (igor@shg.ru)
** \author	Alexey Rubtsov (alex@shg.ru)
** \author	Andrey Antipov (antipov@shg.ru)
*/

#ifndef _CONFIG_H
#define _CONFIG_H

#include <vector>

/** Real floating point type. */
typedef double n_type;
/** Complex type. */
typedef std::complex<n_type> ComplexType;


/** Dense complex matrix. */
//typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign> MatrixType;
/** Dense real matrix. */
//typedef Eigen::Matrix<n_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign> RealMatrixType;

/** Dense complex vector. */
typedef std::vector<ComplexType> VectorType;
/** Dense real vector. */
typedef std::vector<n_type> RealVectorType;
/** Dense vector of integers. */
typedef std::vector<int> IntVectorType;

/** A short name for imaginary unit. */
static const ComplexType I = ComplexType(0.0,1.0);	// 'static' to prevent linking problems
static n_type Pi=4.*atan(n_type(1.0));

/** Generalized 'square' function. */
template<typename T> inline T sqr(T x) { return x*x; }

//@{
/** Do-It-Once environment from A. Rubtsov
**
** When you want a piece of code to run exactly once, just write:
** \verbatim
do_once
	... your code goes here...
end_do_once
\endverbatim
**/
#define do_once { static bool done_once=false; if (!done_once) {done_once=true;
#define end_do_once }; };
//@}

// Maybe temporary
using std::cout;
using std::endl;
using std::flush;
using std::complex;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::vector;
using std::string;
using std::stringstream;

#endif /* #ifndef _CONFIG_H */
