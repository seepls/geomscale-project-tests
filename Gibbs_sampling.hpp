#ifndef GIBBS_MONTE_CARLO_WALK_HPP
#define GIBBS_MONTE_CARLO_WALK_HPP


#include "generators/boost_random_number_generator.hpp"
#include "random_walks/gaussian_helpers.hpp"


struct GibbsMonteCarloWalk {

    template
    <
    	typename Polytope,
    	typename RandomNumberGenerator,
  	>

 	struct Walk 
  	{

    	typedef typename Polytope::PointType Point;
    	typedef typename Point::FT NT;
    	int size =0;



    	template <typename GenericPolytope>
    	Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator& rng)
    	{    	Point x_in = GetDirection<Point>::apply(P.dimension(), rng);
        		Point v_in = GetDirection<Point>::apply(P.dimension(), rng);
    	}
    	template <typename Point>
        inline void update(GenericPolytope const& P, Point const& p1, Point const& p2, RandomNumberGenerator& rng,int i, std::list<Point> &points,Point const & interiorPoint)
        {   dim = p.dimension(); // a point will have as many marginals as dimension .
    		NT p = interiorPoint.getCoefficients();
            const std::vector<vec_t> cov {{8.0f, -1.0f}, {-1.0f, 2.0f}};
            int j = 1- i ;
        	p_hat = gaussian_rand(p1-(p2-p)*cov[i][j] / cov[i][i], std::sqrt(1.0f / cov[i][i]));
        	interiorPoint = p1;

        	points.push_back(Point(p_hat));

        }

        template <typename Point>
        inline void apply(BallPolytope const& P,
                          Point& p1,   // a point to start
                          Point& p2,
                          unsigned int const& walk_length,
                          RandomNumberGenerator& rng,
                          std::list<Point> &points,
                          Point const & interiorPoint)
        {   
            for (auto j=0u; j<walk_length; ++j) 
            {
                if(rng_uniform(0,1)){// 
                	update(P,p1,p2,rng,0,interiorPoint); size++ ; // update p1
                }else{
                	update(P,p2,p1,rng,1, interiorPoint); size++;  // update p2 
                }
            }
        }



        // Discard the burn-in period
        const int n_burnin = points.size() / 10;
        std::rotate(points.begin(), points.begin() + n_burnin, points.end());
        points.resize(points.size() - n_burnin );

    };
};


#endif


