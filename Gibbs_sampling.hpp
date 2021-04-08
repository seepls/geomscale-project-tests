
#include <Eigen/Core>
#include <Eigen/Dense>
#include <armadillo>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>

#include "Sampling_functions.h"


using namespace std;
using namespace Eigen;
using namespace arma;
using namespace boost;


//boost::random::mt19937 gen(0);
boost::random::mt19937 gen(time(0));
//distributions
double runif(double lower, double higher)
{
	boost::random::uniform_real_distribution<> dist(lower, higher);
	return dist(gen);
}
