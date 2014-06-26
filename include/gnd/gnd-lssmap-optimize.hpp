/*
 * gnd-lssmap-optimize.hpp (GND's Laser Scan Statistics MAP OPTIMIZE)
 *
 *  Created on: 2014/09/09
 *      Author: Taichi Yamada
 */


#ifndef GND_LSSMAP_OPTIMIZE_HPP_
#define GND_LSSMAP_OPTIMIZE_HPP_

#include "gnd-lssmap-base.hpp"

#include "gnd-linalg.hpp"
#include "gnd-optimize.hpp"
#include "gnd-random.hpp"
#include "gnd-matrix-coordinate.hpp"


// ---> class declaration
namespace gnd  {
	namespace lssmap {
		class optimize_basic;
		class optimize_newton;
		typedef optimize_newton newton;
		class optimize_monte_calro_method;
		typedef optimize_monte_calro_method mcl;
		class optimize_quasi_monte_calro;
		typedef optimize_quasi_monte_calro qmc;
		class optimize_hybrid_qmc2newton;
		typedef optimize_hybrid_qmc2newton hybrid_q2n;
	}
} // <--- class declaration



//------------------------------------


// ---> function definition
namespace gnd {
	namespace lssmap {
		int gradient( lssmap_t *m, double x, double y, matrix::fixed<4,4> *c, matrix::fixed<3,1> *g, double *l );


		/**
		 * @brief optimization iterate
		 * @param[in]  x : laser scanner reflection point x
		 * @param[in]  y : laser scanner reflection point y
		 * @param[in]  c : coordinate convert matrix
		 * @param[in]  m : map
		 * @param[out] g : gradient   ( derivation of f(x,y) because of minimization problem)
		 * @param[out] l : likelihood ( f(x,y) )
		 * @return    0 :
		 * @details this function compute following value for newton's method:
		 * likelihood, gradient, hessian
		 */
		inline
		int gradient( lssmap_t *m, double x, double y, matrix::fixed<4,4> *c, matrix::fixed<3,1> *g, double *l ) {
			matrix::fixed<4,1> X;				// reflection point on global coordinate
			matrix::fixed<PosDim,3> J;			// jacobi
			matrix::fixed<PosDim,1> q_quad_33;

			{ // ---> initialize
				matrix::set_zero(g);
				if(l) *l = 0;
			} // <--- initialize


			{ // ---> compute reflection points on global coordinate
				matrix::fixed<4,1> xx;				// reflection point on global coordinate
				xx[0][0] = x;
				xx[1][0] = y;
				xx[2][0] = 0;
				xx[3][0] = 1;
				matrix::prod(c, &xx, &X);
			} // <--- compute reflection points on global coordinate

			{ // ---> jacobi matrix
				set_zero(&J);
				J[PosX][0] = 1;
				J[PosX][2] = - x * (*c)[1][0] - y * (*c)[0][0];
				J[PosY][1] = 1;
				J[PosY][2] =   x * (*c)[0][0] - y * (*c)[1][0];
			} // <--- jacobi matrix

			{ // ---> q quad_33
				q_quad_33[PosX][0] = - x * (*c)[0][0] + y * (*c)[1][0];
				q_quad_33[PosY][0] = - x * (*c)[1][0] - y * (*c)[0][0];
			} // <--- q quad_33;

			// ---> scanning loop of map plane
			for( size_t mi = 0; mi < PlaneNum; mi++){
				lssmap_pixel_t *px;
				matrix::fixed<PosDim,1> q;
				matrix::fixed<1,PosDim> qT_invS;		// q^t * Sigma^-1
				matrix::fixed<1,3> qT_invS_J;			// q^t * Sigma^-1 * J
				double chi_sq_2;						// q^T * Sigma^-1 * q / 2 / 2value
				double lkh = 0;
				matrix::fixed<1,3> gg;					// gradient
				matrix::fixed<1,1> ws1x1;				// work space

				{ // ---> get pixel data and compute q(difference from mean)
					matrix::fixed<2,1> ws2x1;
					long pr, pc;
					// get index of pixel
					if( m->plane[mi].pindex( X[PosX][0], X[PosY][0], &pr, &pc ) < 0)
						continue; // no data
					px = m->plane[mi].pointer( pr, pc );
					// zero weight
					if(px->K <= 0.0)	continue;
					// get pixel core position
					m->plane[mi].pget_pos_core(pr, pc, pointer(&ws2x1, PosX, 0), pointer(&ws2x1, PosY, 0));
					// compute sensor reading position on a focus pixel
					sub(&X, &ws2x1, &ws2x1);
					// compute difference from mean
					sub(&ws2x1, &px->mean, &q);
				} // <--- get pixel data and compute q(difference from mean)


				{ // ---> likelihood
					// q^T * Sigma^-1
					prod_transpose1(&q, &px->inv_cov, &qT_invS);
					// chi-square / 2.0
					prod(&qT_invS, &q, &ws1x1);
					chi_sq_2 = ws1x1[0][0] / 2.0;
					// likelihood exp( -q^T * Sigma^-1 * q / 2 )
					lkh = ::exp( -chi_sq_2 ) * px->K;
					if( l ) (*l) += lkh;
				} // <--- likelihood

				{ // ---> gradient
					// q^T * Sigma^-1 * q' * exp( - q^T * Sigma^-1 * q / 2)
					prod(&qT_invS, &J, &qT_invS_J);
					// gradient
					matrix::scalar_prod(&qT_invS_J, lkh, &gg);
					matrix::sub_transpose2(g, &gg, g);
				} // <--- gradient
			}
			return 0;
		}
	}
} // <--- function definition


// ---> class definition
namespace gnd {
	namespace lssmap {



		/**
		 * @ingroup GNDPSM
		 * @brief probabilistic scan matching optimizer
		 */
		class optimize_basic {
		public:
			/// @brief position type
			typedef matrix::fixed<3,1> 			pos_t;

			// ---> declaration
		public:
			/// map type
			typedef lssmap::lssmap_t					map_t;
			/// map type pointer
			typedef map_t*						map_pt;
			/// pixel type
			typedef lssmap::lssmap_pixel_t 				pixel_t;
			/// pixel type pointer
			typedef pixel_t*					pixel_pt;
			/// reflection point
			typedef matrix::fixed<2,1> 			point_t;
			/// reflection points
			typedef queue< matrix::fixed<2,1> > points_t;
			// <--- declaration


			// ---> constructor, destructor
		public:
			optimize_basic();
			optimize_basic(map_pt m);
			~optimize_basic();
			// <--- constructor, destructor


			// ---> map
		protected:
			map_pt _map;	///< @brief map reference
		public:
			virtual int set_map(map_pt m);
			virtual int release_map();
			// <--- map

			// ---> reflection point
		protected:
			/// @brief laser scanner reflection point
			points_t _points;
		public:
			virtual int set_scan_point(double x, double y);
			virtual int nscan_point() const;
			// <--- reflection point


			// ---> starting vlaue of optimization
		public:
			/**
			 * @brief create initial paramete
			 */
			virtual int initial_parameter_create(void** p) = 0;
			/**
			 * @brief set starting value
			 */
			virtual int initial_parameter_set_position(void* p, double x, double y, double theta) = 0;
			/**
			 * @brief delete starting value
			 */
			virtual int initial_parameter_delete(void** p) = 0;
			// <--- starting value of optimization


			// ---> optimization
		public:
			/**
			 * @brief begin
			 * @note clear all (scan data, previous optimization result) and set initial parameter
			 * @param[in] v : initial parameter
			 * @return    0 :
			 */
			virtual int begin(void *v) = 0;
			/**
			 * @brief iterate
			 * @param[out] d : delta
			 * @param[out] p : pos
			 * @param[out] l : likelihood
			 * @return    0 :
			 */
			virtual int iterate(matrix::fixed<3,1> *d, matrix::fixed<3,1> *p, double *l) = 0;
			// <--- optimization


			// ---> converge criteria
		public:
			/// @brief default convergenct threshold (distance)
			static const double DefConvergeDist = gnd_m2dist(0.01);
			/// @brief default convergenct threshold (orient)
			static const double DefConvergeAngle = gnd_deg2ang(0.5);
			/// @brief variables ofconvergence test
			struct converge_var {
				double sqdist;	///< distance threshold
				double orient;	///< orient threshold
				pos_t delta;	///< delta
				converge_var();
			};
		protected:
			/// converge parameter
			struct converge_var _converge;
		public:
			int set_converge_threshold(double d, double a);
			int converge_test();
			static int converge_test(double sqd, double o, double sqdt, double ot);
			// <--- converge criteria
		};

		/**
		 * @brief constructor
		 */
		inline
		optimize_basic::optimize_basic()
		{
			_map = 0;
		}

		/**
		 * @brief constructor
		 * @param[in] m: map
		 */
		inline
		optimize_basic::optimize_basic(map_pt m) {
			set_map(m);
		}

		/**
		 * @brief destructor
		 */
		inline
		optimize_basic::~optimize_basic()
		{
			release_map();
		}


		/**
		 * @brief set map
		 * @param[in] m : map pointer
		 * @return  0 : map is null
		 * 		   >0 : map is exist
		 */
		inline
		int optimize_basic::set_map(map_pt m)
		{
			if(m)	_map = m;
			return (_map != 0);
		}

		/**
		 * @brief release map
		 */
		inline
		int optimize_basic::release_map()
		{
			_map = 0;
			return 0;
		}


		/**
		 * @brief set scan point
		 * @param[in]  x : reflect point x
		 * @param[in]  y : reflect point y
		 * @return    number of points
		 */
		inline
		int optimize_basic::set_scan_point(double x, double y) {
			point_t p;
			p[0][0] = x;
			p[1][0] = y;
			return _points.push_back(&p);
		}

		/**
		 * @brief get number of scan point
		 * @return    number of points
		 */
		inline
		int optimize_basic::nscan_point() const{
			return _points.size();
		}


		/**
		 * @brief constructor of optimizer_basic::converge_var
		 */
		optimize_basic::converge_var::converge_var() : sqdist( gnd_square( DefConvergeDist )), orient(DefConvergeAngle)
		{
		}

		/**
		 * @brief set convergence criteria threshold
		 * @param[in] d: distance
		 * @param[in] o: orient
		 */
		int optimize_basic::set_converge_threshold(double d, double o) {
			_converge.sqdist = gnd_square(d);
			_converge.orient = ::fabs(o);
			return 0;
		}

		/**
		 * @brief converge test
		 */
		int optimize_basic::converge_test() {
			return converge_test( gnd_square( _converge.delta[0][0] ) + gnd_square( _converge.delta[1][0] ) , _converge.delta[2][0],
					_converge.sqdist, _converge.orient);
		}


		/**
		 * @brief converge test
		 * @param[in]  sqd : square distance
		 * @param[in]    o : orient
		 * @param[in] sqdt : square distance threshold
		 * @param[in]   ot : orient threshold
		 * @return !0 : converge
		 * @return  0 : not converge
		 */
		int optimize_basic::converge_test(double sqd, double o, double sqdt, double ot) {
			return sqd < sqdt && ::fabs(o) < ot;
		}

	}
}
// optimizer_basic
// ---> class definition



// ---> class definition
// optimizer_newton
namespace gnd {
	namespace lssmap {

		/**
		 * @ingroup GNDPSM
		 * @brief probabilistic scan matching optimizer
		 */
		class optimize_newton
				: public optimize_basic {

			// ---> type declaration
				public:
			/// map type
			typedef lssmap::lssmap_t					map_t;
			/// map type pointer
			typedef map_t*						map_pt;
			/// pixel type
			typedef lssmap::lssmap_pixel_t 				pixel_t;
			/// pixel type pointer
			typedef pixel_t*					pixel_pt;

			/// input variance
			typedef matrix::fixed<3,1>			initial_parameter;
			/// coordinate matrix type
			typedef matrix::fixed<4,4>			coordm_t;

			/// @brief variables
			struct variables {
				matrix::fixed<3,1> pos;		///< position
				double likelihood;			///< likelihood
			};
			/// variables
			typedef struct variables variables;
			// <--- type declaration


			// ---> constructor, destructor
				public:
			optimize_newton();
			optimize_newton(map_pt m);
			~optimize_newton();
			// <--- constructor, destructor

			// ---> variance
				protected:
			/// @brief variables
			variables	_var;
			/// @brief coordinate convert matrix
			coordm_t	_coordm;
			// <--- variance


			// ---> starting value of optimization
				public:
			virtual int initial_parameter_create(void** p);
			virtual int initial_parameter_delete(void** p);
			virtual int initial_parameter_set_position(void* p, double x, double y, double theta);
			// <--- starting value of optimization



			// ---> optimization
				public:
			virtual int begin(void *v);
			virtual int iterate(matrix::fixed<3,1> *d, matrix::fixed<3,1> *p, double *l) ;
			// <--- optimization

			static int _newton_method_variables_( double x, double y, matrix::fixed<4,4> *c, map_pt m,  double *l, matrix::fixed<3,1> *g, matrix::fixed<3,3> *h );
		};


		/**
		 * @brief constructor
		 */
		inline
		optimize_newton::optimize_newton()
		{
		}

		/**
		 * @brief constructor
		 */
		inline
		optimize_newton::optimize_newton(map_pt m) : optimize_basic(m)
		{
		}

		/**
		 * @brief destructor
		 */
		inline
		optimize_newton::~optimize_newton()
		{
		}


		/**
		 * @brief create starting value
		 * @param[in,out] p : starting value buffer pointer
		 */
		inline
		int optimize_newton::initial_parameter_create(void** p) {
			*p = static_cast<void*>( new initial_parameter );
			return 0;
		}


		/**
		 * @brief delete starting value
		 * @param[in,out] p : starting value
		 */
		inline
		int optimize_newton::initial_parameter_delete(void** p) {
			initial_parameter* pp = static_cast<initial_parameter*>( *p );
			delete pp;
			*p = 0;
			return 0;
		}


		/**
		 * @brief set starting value
		 * @param[in]     p : starting value
		 * @param[in]     x : robot position x
		 * @param[in]     y : robot position y
		 * @param[in] theta : robot position theta
		 */
		inline
		int optimize_newton::initial_parameter_set_position(void* p, double x, double y, double theta) {
			initial_parameter *pp = static_cast<initial_parameter*>(p);
			(*pp)[0][0] = x;
			(*pp)[1][0] = y;
			(*pp)[2][0] = theta;
			return 0;
		}

		/**
		 * @brief input starting value
		 * @param[in] v : starting value
		 */
		inline
		int optimize_newton::begin(void *v) {
			gnd_assert(!v, -1, "invalid null argument");
			initial_parameter *p = static_cast<initial_parameter*>(v);

			// set position
			copy(&_var.pos, p);

			matrix::coordinate_converter(&_coordm,
					_var.pos[0][0], _var.pos[1][0], 0,
					::cos(_var.pos[2][0]), ::sin(_var.pos[2][0]), 0,
					 0, 0, 1);

			// set zero
			_var.likelihood = 0;
			_points.clear();
			return 0;
		}

		/**
		 * @brief optimization iterate
		 * @param[out] d : difference
		 * @param[out] p : pos
		 * @param[out] v : variance
		 * @param[out] l : likelihood
		 * @return    0 :
		 */
		inline
		int optimize_newton::iterate(matrix::fixed<3,1> *d, matrix::fixed<3,1> *p, double *l) {
			double likelihood = 0;
			uint64_t pi;
			matrix::fixed<3,1> grad;
			matrix::fixed<3,3> hess;

			// ---> scanning loop of reflection points
			for( pi = 0; pi < _points.size(); pi++ ){
				double l;
				matrix::fixed<3,1> g;
				matrix::fixed<3,3> h;

				// compute likelihood, gradient, hessian
				_newton_method_variables_( _points[pi][0][0], _points[pi][1][0], &_coordm, _map, &l, &g, &h );
				// summation
				likelihood += l;
				add(&grad, &g, &grad);
				add(&hess, &h, &hess);
			} // <--- scanning loop of reflection points


			{ // ---> optimization
				int ret;
				matrix::fixed<3,1> delta;

				// modified newton's method for unconstrained minimizer
				if( (ret = optimize::newtons_method_unconstrainted(&grad, &hess, &delta)) < 0)
					return ret;

				// store delta for convergence test
				copy(&_converge.delta, &delta);
				// position optimization
				add(&delta, &_var.pos, &_var.pos);

				// set coordinate convert matrix
				matrix::coordinate_converter(&_coordm,
						_var.pos[0][0], _var.pos[1][0], 0,
						::cos(_var.pos[2][0]), ::sin(_var.pos[2][0]), 0,
						 0, 0, 1);

				LogDebugf("     : newton - delta = (%lf, %lf, %lf):\n", delta[0][0], delta[1][0], delta[2][0]);
				// set output data
				if( d ) {
					copy(d, &delta);
				}
				if( p ) {
					copy(p, &_var.pos);
				}
				if( l ) {
					*l = _var.likelihood;
				}
			} // <--- set output

			{ // ---> set zero
				_var.likelihood = 0;
			} // <--- set zero

			return 0;
		}



		/**
		 * @brief optimization iterate
		 * @param[in]  x : laser scanner reflection point x
		 * @param[in]  y : laser scanner reflection point y
		 * @param[in]  c : coordinate convert matrix
		 * @param[in]  m : likelihood map
		 * @param[out] l : likelihood ( f(x,y) )
		 * @param[out] g : gradient   ( derivation of -f(x,y) because of minimization problem)
		 * @param[out] h : hessian    ( quadratic of -f(x,y) because of minimization problem)
		 * @return    0 :
		 * @details this function compute following value for newton's method:
		 * likelihood, gradient, hessian
		 */
		inline
		int optimize_newton::_newton_method_variables_( double x, double y, matrix::fixed<4,4> *c, map_pt m,  double *l, matrix::fixed<3,1> *g, matrix::fixed<3,3> *h ) {
			matrix::fixed<4,1> X;				// reflection point on global coordinate
			matrix::fixed<PosDim,3> J;			// jacobi
			matrix::fixed<PosDim,1> q_quad_33;

			{ // ---> compute reflection points on global coordinate
				matrix::fixed<4,1> xx;
				xx[0][0] = x;
				xx[1][0] = y;
				xx[2][0] = 0;
				xx[3][0] = 1;
				prod( c, &xx, &X );
			} // <--- compute reflection points on global coordinate

			{ // ---> jacobi matrix
				set_zero(&J);
				J[PosX][0] = 1;
				J[PosX][2] = - x * (*c)[1][0] - y * (*c)[0][0];
				J[PosY][1] = 1;
				J[PosY][2] =   x * (*c)[0][0] - y * (*c)[1][0];
			} // <--- jacobi matrix

			{ // ---> q quad_33
				q_quad_33[PosX][0] = - x * (*c)[0][0] + y * (*c)[1][0];
				q_quad_33[PosY][0] = - x * (*c)[1][0] - y * (*c)[0][0];
			} // <--- q quad_33;

			{
				*l = 0;
				matrix::set_zero(g);
				matrix::set_zero(h);
			}

			// ---> scanning loop of map plane
			for( size_t mi = 0; mi < PlaneNum; mi++){
				pixel_t *px;
				matrix::fixed<3,1> q;
				matrix::fixed<1,3> qT_invS;		// q^t * Sigma^-1
				matrix::fixed<1,3> qT_invS_J;			// q^t * Sigma^-1 * J
				double chi_sq_2;						// q^T * Sigma^-1 * q / 2 / 2value
				double lkh = 0;
				matrix::fixed<1,3> gg;					// gradient
				matrix::fixed<3,3> hh;					// gradient
				matrix::fixed<1,1> ws1x1;				// work space

				{ // ---> get pixel data and compute q(difference from mean)
					matrix::fixed<2,1> ws2x1;
					long pr, pc;
					// get index of pixel
					if( m->plane[mi].pindex( X[PosX][0], X[PosY][0], &pr, &pc ) < 0)
						continue; // no data
					px = m->plane[mi].pointer( pr, pc );
					// zero weight
					if(px->K <= 0.0)	continue;
					// get pixel core position
					m->plane[mi].pget_pos_core(pr, pc, pointer(&ws2x1, PosX, 0), pointer(&ws2x1, PosY, 0));
					// compute sensor reading position on a focus pixel
					sub(&X, &ws2x1, &ws2x1);
					// compute difference from mean
					sub(&ws2x1, &px->mean, &q);
				} // <--- get pixed data and compute q(difference from mean)


				{ // ---> likelihood
					// q^T * Sigma^-1
					prod_transpose1(&q, &px->inv_cov, &qT_invS);
					// chi-square / 2.0
					prod(&qT_invS, &q, &ws1x1);
					chi_sq_2 = ws1x1[0][0] / 2.0;
					// likelihood exp( -q^T * Sigma^-1 * q / 2 )
					lkh = ::exp( -chi_sq_2 ) * px->K;
					*l += lkh;
				} // <--- likelihood

				{ // ---> gradient
					// q^T * Sigma^-1 * q' * exp( - q^T * Sigma^-1 * q / 2)
					matrix::prod(&qT_invS, &J, &qT_invS_J);
					// gradient
					scalar_prod(&qT_invS_J, lkh, &gg);
					matrix::add(g, &gg, g);
				} // <--- gradient

				{ // ---> hessian
					matrix::fixed<3,3> ws3x3;
					matrix::fixed<3,2> ws3x2;

					// (-q*Sigma^-1*q')^T * (-q*Sigma^-1*q')
					prod_transpose1( &qT_invS_J, &qT_invS_J, &hh );

					// -q*Sigma^-1*q''
					matrix::prod( &qT_invS, &q_quad_33, &ws1x1 );
					hh[2][2] -= ws1x1[0][0];

					// -q'*Sigma^-1*q'
					matrix::prod_transpose1( &J, &px->inv_cov, &ws3x2 );
					matrix::prod( &ws3x2, &J, &ws3x3 );
					matrix::sub( &hh, &ws3x3, &hh );
					// hessian
					scalar_prod( &hh, -(lkh), &hh );

					optimize::_model_hessian_(&hh);
					matrix::add(h, &hh, h);
				} // <--- hessian
			}
			return 0;
		}

	}

};
// <--- class function definition











// ---> class definition
namespace gnd {
	namespace lssmap {

		/**
		 * @brief probabilistic scan matching optimizer (monte-calro)
		 */
		class optimize_monte_calro_method : public optimize_basic {
			// ---> type declaration
		public:
			/// map type
			typedef lssmap::lssmap_t		map_t;
			/// map type pointer
			typedef map_t*			map_pt;
			/// pixel type
			typedef lssmap::lssmap_pixel_t	pixel_t;
			/// pixel type pointer
			typedef pixel_t*		pixel_pt;
			// <--- type declaration

			// ---> constructor, destructor
		public:
			optimize_monte_calro_method();
			optimize_monte_calro_method(map_pt m);
			~optimize_monte_calro_method();
			// <--- constructor, destructor

			// ---> particle
		protected:
			/// @brief particle
			struct particle {
				vector::fixed_column<3> pos;	///< position
				matrix::fixed<4,4> coordm;		///< coordinate convert matrix
				double likelihood;				///< likelihood
			};
			/// @brief particles
			queue< struct particle > particles;
			/// @brief particles
			queue< struct particle > _ws_resmpl;
			// <--- particle

			// ---> starting value
		public:
			/// @brief starting value of optimization
			typedef struct initial_parameter {
				vector::fixed_column<3> pos;		///<! position (column vector)
				matrix::fixed<3,3>	var_ini;		///<! initial variance
				matrix::fixed<3,3>	var_rsmp;		///<! variance for resampling random
				uint32_t			n;				///<! number of particle
				double				alpha;			///<! parameter for seeking around max likelihood particle
				initial_parameter();
			} initial_parameter;

		protected:
			/// @brief optimize
			initial_parameter _v;
			// <--- starting value


			// ---> starting value of optimization
		public:
			virtual int initial_parameter_create(void** p);
			virtual int initial_parameter_delete(void** p);
			virtual int initial_parameter_set_position(void* p, double x, double y, double theta);
			// <--- starting of optimization


			// ---> optimization
		public:
			virtual int begin(void *v);
			virtual int iterate(matrix::fixed<3,1> *d, matrix::fixed<3,1> *p, double *l);
			// <--- optimization
		protected:
			int init_particle( initial_parameter *v );
		};


		/**
		 * @brief constructor
		 */
		inline
		optimize_monte_calro_method::optimize_monte_calro_method() {
			return;
		}

		/**
		 * @brief constructor
		 * @param[in] m : map
		 */
		inline
		optimize_monte_calro_method::optimize_monte_calro_method(map_pt m)
		: optimize_basic(m) {
			return;
		}

		/**
		 * @brief destructor
		 */
		inline
		optimize_monte_calro_method::~optimize_monte_calro_method(){
		}


		/**
		 * @brief create starting value
		 * @param[in,out] p : starting value buffer pointer
		 */
		inline
		int optimize_monte_calro_method::initial_parameter_create(void** p) {
			*p = static_cast<void*>( new initial_parameter );
			return 0;
		}


		/**
		 * @brief delete starting value
		 * @param[in,out] p : starting value
		 */
		inline
		int optimize_monte_calro_method::initial_parameter_delete(void** p) {
			initial_parameter* pp = static_cast<initial_parameter*>( *p );
			delete pp;
			*p = 0;
			return 0;
		}


		/**
		 * @brief set starting value
		 * @param[in]     p : starting value
		 * @param[in]     x : robot position x
		 * @param[in]     y : robot position y
		 * @param[in] theta : robot position theta
		 */
		inline
		int optimize_monte_calro_method::initial_parameter_set_position(void* p, double x, double y, double theta) {
			initial_parameter *pp = static_cast<initial_parameter*>(p);
			pp->pos[0] = x;
			pp->pos[1] = y;
			pp->pos[2] = theta;
			return 0;
		}



		/**
		 * @brief input starting value
		 * @param[in] v : input variance
		 * @return 0
		 */
		inline
		int optimize_monte_calro_method::begin(void *v) {
			LogDebug("Begin - mcl begin\n");
			LogIndent();

			::memcpy(&_v, v, sizeof(_v));

			// ckolesky
			linalg::cholesky_decomposition(&_v.var_rsmp, 3);

			// reflesh particle
			particles.clear();
			init_particle(&_v);

			// clear laser scanner reflection point
			_points.clear();

			LogUnindent();
			LogDebug("End   - mcl begin\n");

			return 0;
		}


		/**
		 * @brief create particle
		 * @param[in] v : input variance
		 * @return 0
		 */
		inline
		int optimize_monte_calro_method::init_particle( initial_parameter *v ) {
			matrix::fixed<3,1> rnd;
			matrix::fixed<3,3> L;
			particle p;

			// cholesky decomposition
			matrix::copy(&L, &v->var_ini);
			linalg::cholesky_decomposition(&L, 3);

			// set 0
			p.likelihood = 0;

			// set position
			matrix::copy(&p.pos, &v->pos);

			// get coordinate convert matrix
			matrix::coordinate_converter(&p.coordm,
					p.pos[0], p.pos[1], 0,
					::cos(p.pos[2]), ::sin(p.pos[2]), 0,
					 0, 0, 1);
			particles.push_back(&p);

			while ( particles.size() < v->n ){
				// create particle according to gaussian
				rnd[0][0] = random_gaussian(1.0);
				rnd[1][0] = random_gaussian(1.0);
				rnd[2][0] = random_gaussian(1.0);
				matrix::prod(&L, &rnd, &p.pos);
				// add average
				matrix::add(&p.pos, &v->pos, &p.pos);

				// get coordinate convert matrix
				matrix::coordinate_converter(&p.coordm,
						p.pos[0], p.pos[1], 0,
						::cos(p.pos[2]), ::sin(p.pos[2]), 0,
						 0, 0, 1);

				particles.push_back(&p);
			}

			return 0;
		}



		/**
		 * @brief optimization iterate
		 * @param[out] d : difference
		 * @param[out] p : pos
		 * @param[out] v : variance ()
		 * @param[out] l : likelihood
		 * @return    0 :
		 */
		inline
		int optimize_monte_calro_method::iterate(matrix::fixed<3,1> *d, matrix::fixed<3,1> *p, double *l){
			gnd_error(particles.size() == 0, -1, "no data" );
			gnd_error(_points.size() == 0, -1, "no scan data" );

			LogDebug("Begin - mcl iterate\n");
			LogIndent();

			{ // ---> operate
				uint32_t i, j;
				double lk = 0;
				double sum;
				double max;
				size_t imax;
				vector::fixed_column<3>	delta;	// delta
				vector::fixed_column<3>	ws3x1;	// workspace
				matrix::fixed<3,3>		ws3x3;	// workspace

				matrix::set_zero(&delta);
				sum = 0;
				for( j = 0; j < _points.size(); j++ ){
					matrix::fixed<4,1> p;

					p[0][0] = _points[j][0][0];
					p[1][0] = _points[j][1][0];
					p[2][0] = 0;
					p[3][0] = 1;

					for(size_t i = 0; i < (unsigned)particles.size(); i++){
						matrix::fixed<4,1> x;

						matrix::prod( &particles[i].coordm,  &p, &x );

						lssmap::likelihood(_map, x[0][0], x[1][0], &lk);
						particles[i].likelihood += ::log( (double)lk + DBL_EPSILON );
					}
				} // ---> loop for compute likelihood

				// ---> get max
				max = -DBL_MAX;;
				imax = 0;
				for( i = 0 ; i <  (unsigned)particles.size(); i++ ) {
					// get max
					if( max < particles[i].likelihood )	{
						max = particles[i].likelihood;
						imax = i;
					}
				}
				// <--- get max

				// ---> likelihood sum
				for( i = 0; i < (unsigned)particles.size(); i++){
					particles[i].likelihood = ::exp( particles[i].likelihood - max + DBL_EPSILON) + DBL_EPSILON;

					// sum
					sum += particles[i].likelihood;

					LogVerbosef("      : particle (%lf, %lf, %lf) lkh = %lf\n",
							particles[i].pos[0], particles[i].pos[1], particles[i].pos[2], particles[i].likelihood);
				} // ---> likelihood sum

				// if sum == zero, error
				gnd_error( sum == 0, -1, "all likelihood is 0" );

				{ // ---> save maximum likelihood particle
					matrix::sub(&particles[imax].pos, &_v.pos, &delta);

					matrix::copy(&_v.pos, &particles[imax].pos);
				} // <--- save maximum likelihood particle

				LogDebugf("      : delta = (%lf, %lf, %lf):\n", delta[0], delta[1], delta[2]);

				// store delta for converge test
				matrix::copy(&_converge.delta, &delta);

				{ // ---> output
					if( p ) { // ---> set position
						matrix::copy(p, &particles[imax].pos);
					} // <--- set position

					if( d ) { // ---> set delta
						matrix::copy(d, &delta);
					} // <--- set delta

					if( l ) {// ---> set likelihood
						*l = particles[imax].likelihood;
					} // <--- set likelihood
				} // <--- output (nearest neighbor)



				{ // ---> resampling
					particle p;

					// clear
					_ws_resmpl.clear();

					// save max likelihood particle
					_ws_resmpl.push_back(&particles[imax]);

					// save max particle
					// and seek around max one
					while(_ws_resmpl.size() < _v.n * _v.alpha )
						_ws_resmpl.push_back(&particles[imax]);

					// ---> select particle considering its likelihood
					while( _ws_resmpl.size() < _v.n ) {
						double lrand = sum * random_uniform();

						// select sample
						for( i = 0; i < particles.size() - 1 && lrand > 0; i++ )
							lrand -= particles[i].likelihood;

						// save
						_ws_resmpl.push_back(&particles[i]);
					} // <--- select particle considering its likelihood


					// cholesky decomposition
					LogVerbosef("L    : %.4lf, %.4lf, %.4lf\n", _v.var_rsmp[0][0], _v.var_rsmp[0][1], _v.var_rsmp[0][2] );
					LogVerbosef("     : %.4lf, %.4lf, %.4lf\n", _v.var_rsmp[1][0], _v.var_rsmp[1][1], _v.var_rsmp[1][2] );
					LogVerbosef("     : %.4lf, %.4lf, %.4lf\n", _v.var_rsmp[2][0], _v.var_rsmp[2][1], _v.var_rsmp[2][2] );

					// clear
					particles.clear();
					p.likelihood = 0;

					// set max
					particles.push_back(&_ws_resmpl[0]);

					// ---> add random value following error
					for( i = particles.size(); i < _ws_resmpl.size(); i++ ) {
						// random
						ws3x1[0] = random_gaussian(1.0);
						ws3x1[1] = random_gaussian(1.0);
						ws3x1[2] = random_gaussian(1.0);

						matrix::prod( &_v.var_rsmp, &ws3x1, &p.pos );
						LogVerbosef("rand : %.4lf, %.4lf, %.4lf\n", p.pos[0], p.pos[1], p.pos[2] );
						matrix::add( &p.pos, &_ws_resmpl[i].pos, &p.pos );

						// get coordinate convert matrix
						matrix::coordinate_converter(&p.coordm,
								p.pos[0], p.pos[1], 0,
								::cos(p.pos[2]), ::sin(p.pos[2]), 0,
								 0, 0, 1);

						// set
						particles.push_back(&p);
					} // ---> add random value following error

				} // <--- resampling
			} // <--- operate

			LogUnindent();
			LogDebug("End  - mcl iterate\n");
			return 0;
		}


		inline
		optimize_monte_calro_method::initial_parameter::initial_parameter(){
			matrix::set_zero(&pos);

			n = 250;
			var_ini[0][0] = gnd_square( gnd_m2dist( 0.1 ) );
			var_ini[1][1] = gnd_square( gnd_m2dist( 0.1 ) );
			var_ini[2][2] = gnd_square( gnd_deg2ang( 10.0 ) );

			var_rsmp[0][0] = gnd_square( gnd_m2dist( 0.005 ) );
			var_rsmp[1][1] = gnd_square( gnd_m2dist( 0.005 ) );
			var_rsmp[2][2] = gnd_square( gnd_deg2ang( 0.25 ) );

			alpha = 0.0;
		}
	}
};
// <--- class definition





// ---> class definition
namespace gnd {
	namespace lssmap {

		/**
		 * @brief probabilistic scan matching optimizer (quasi-monte-calro)
		 */
		class optimize_quasi_monte_calro : public optimize_basic {
			// ---> type declaration
		public:
			/// map type
			typedef lssmap::lssmap_t		map_t;
			/// map type pointer
			typedef map_t*			map_pt;
			/// pixel type
			typedef lssmap::lssmap_pixel_t	pixel_t;
			/// pixel type pointer
			typedef pixel_t*		pixel_pt;
			// <--- type declaration

			// ---> constructor, destructor
		public:
			optimize_quasi_monte_calro();
			optimize_quasi_monte_calro(map_pt m);
			~optimize_quasi_monte_calro();
			// <--- constructor, destructor

		private:
			/// @brief particle table
			queue< vector::fixed_column<3> > table;
			/// @brief create sphere-table
			int create_particletable(uint32_t n);

			// ---> particle
		protected:
			/// @brief particle
			struct particle {
				vector::fixed_column<3> pos;		///< position
				matrix::fixed<4,4> coordm;	///< coordinate convert matrix
				double likelihood;			///< likelihood
			};
			/// @brief particles
			queue< struct particle > particles;
			// <--- particle

			// ---> starting value
		public:
			/// @brief starting value of optimization
			typedef struct initial_parameter {
				vector::fixed_column<3>	pos;			///<! position
				matrix::fixed<3,3>	var;		///<! variance
				uint32_t			n;				///<! resolution
				initial_parameter();
			} initial_parameter;
		protected:
			/// @brief optimize starting value
			initial_parameter _v;
			// <--- starting value


			// ---> starting value of optimization
		public:
			virtual int initial_parameter_create(void** p);
			virtual int initial_parameter_delete(void** p);
			virtual int initial_parameter_set_position(void* p, double x, double y, double theta);
			// <--- starting of optimization


			// ---> optimization
		public:
			virtual int begin(void *v);
			virtual int iterate(matrix::fixed<3,1> *d, matrix::fixed<3,1> *p, double *l);
			// <--- optimization
		protected:
			int create_particle(initial_parameter *v);
		};


		/**
		 * @brief constructor
		 */
		inline
		optimize_quasi_monte_calro::optimize_quasi_monte_calro() {
			return;
		}

		/**
		 * @brief constructor
		 * @param[in] m : map
		 */
		inline
		optimize_quasi_monte_calro::optimize_quasi_monte_calro(map_pt m)
		: optimize_basic(m) {
			return;
		}

		/**
		 * @brief destructor
		 */
		inline
		optimize_quasi_monte_calro::~optimize_quasi_monte_calro(){
		}

		inline
		int optimize_quasi_monte_calro::create_particletable(uint32_t n) {
			uint32_t x, y, t;
			vector::fixed_column<3> tmp;

			LogVerbose("Begin : qmc create_particletable\n");
			LogIndent();

			table.clear();
			matrix::set_zero(&tmp);
			table.push_back(&tmp);

			for(x = 0; x < n * 2 + 3; x++) {
				tmp[0] = ((double)x / (n+1)) - 1;
				for(y = 0; y < n * 2 + 3; y++) {
					tmp[1] = ((double)y / (n+1)) - 1;
					for(t = 0; t < n * 2 + 3; t++) {
						tmp[2] = ((double)t / (n+1)) - 1;
						LogVerbosef("%.04lf %.04lf %.04lf\n", tmp[0], tmp[1], tmp[2]);
						table.push_back(&tmp);
					} // for (t)
				} // for (y)
			} // for (x)

			LogUnindent();
			LogVerbose(" End   : qmc create_particletable\n");
			return 0;
		}



		/**
		 * @brief create starting value
		 * @param[in,out] p : starting value buffer pointer
		 */
		inline
		int optimize_quasi_monte_calro::initial_parameter_create(void** p) {
			*p = static_cast<void*>( new initial_parameter );
			return 0;
		}


		/**
		 * @brief delete starting value
		 * @param[in,out] p : starting value
		 */
		inline
		int optimize_quasi_monte_calro::initial_parameter_delete(void** p) {
			initial_parameter* pp = static_cast<initial_parameter*>( *p );
			delete pp;
			p = 0;
			return 0;
		}


		/**
		 * @brief set starting value
		 * @param[in]     p : starting value
		 * @param[in]     x : robot position x
		 * @param[in]     y : robot position y
		 * @param[in] theta : robot position theta
		 */
		inline
		int optimize_quasi_monte_calro::initial_parameter_set_position(void* p, double x, double y, double theta) {
			initial_parameter *pp = static_cast<initial_parameter*>(p);
			pp->pos[0] = x;
			pp->pos[1] = y;
			pp->pos[2] = theta;
			return 0;
		}

		/**
		 * @brief input starting value
		 * @param[in] v : input variance
		 * @return 0
		 */
		inline
		int optimize_quasi_monte_calro::begin(void *v) {
			initial_parameter *varp = static_cast<initial_parameter*>(v);

			create_particletable(varp->n);

			// reflesh particle
			particles.clear();
			create_particle(varp);

			// clear laser scanner reflection point
			_points.clear();
			::memcpy( &_v, varp, sizeof(_v) );

			return 0;
		}


		/**
		 * @brief create particle
		 * @param[in] v : input variance
		 * @return 0
		 */
		inline
		int optimize_quasi_monte_calro::create_particle(initial_parameter *v) {
			matrix::fixed<3,3> L;
			particle p;
			uint32_t i;

			LogVerbose("Begin : qmc create_particle\n");
			LogIndent();

			// cholesky decomposition
			matrix::copy(&L, &v->var);
			linalg::cholesky_decomposition(&L, 3);
			LogVerbosef("L : %.04lf %.04lf %.04lf\n", L[0][0], L[0][1], L[0][2]);
			LogVerbosef("  : %.04lf %.04lf %.04lf\n", L[1][0], L[1][1], L[1][2]);
			LogVerbosef("  : %.04lf %.04lf %.04lf\n", L[2][0], L[2][1], L[2][2]);

			// set 0
			p.likelihood = 0;
			matrix::set_zero(&p.pos);
			matrix::coordinate_converter(&p.coordm,
					p.pos[0] + v->pos[0], p.pos[1] + v->pos[1], 0,
					::cos(p.pos[2] + v->pos[2]), ::sin(p.pos[2] + v->pos[2]), 0,
					 0, 0, 1);
			particles.push_back(&p);

			// ---> scanning loop for Table
			for( i = 0; i < table.size(); i++ ){
				prod(&L, table + i, &p.pos);
				matrix::coordinate_converter(&p.coordm,
						p.pos[0] + v->pos[0], p.pos[1] + v->pos[1], 0,
						::cos(p.pos[2] + v->pos[2]), ::sin(p.pos[2] + v->pos[2]), 0,
						 0, 0, 1);
				LogVerbosef("%.04lf %.04lf %.04lf\n", p.pos[0], p.pos[1], p.pos[2]);
				particles.push_back(&p);
			} // <--- scanning loop for SphereTable

			LogUnindent();
			LogVerbose(" End   : qmc create_particle\n");
			return 0;
		}



		/**
		 * @brief optimization iterate
		 * @param[out] d : difference
		 * @param[out] p : pos
		 * @param[out] v : variance
		 * @param[out] l : likelihood
		 * @return    0 :
		 */
		inline
		int optimize_quasi_monte_calro::iterate(matrix::fixed<3,1> *d, matrix::fixed<3,1> *p, double *l){
			gnd_error(particles.size() == 0, -1, "no data" );

			{ // ---> operate
				double lk = 0;
				double sum;
				matrix::fixed<3,1> delta;
				matrix::fixed<3,1> ws3x1;

				set_zero(&delta);
				sum = 0;
				// ---> loop for compute likelihood
				for(size_t i = 0; i < (unsigned)particles.size(); i++){
					particles[i].likelihood = 0;
				}
				for(size_t j = 0; j < (unsigned)_points.size(); j++ ){
					matrix::fixed<4,1> p;

					p[0][0] = _points[j][0][0];
					p[1][0] = _points[j][1][0];
					p[2][0] = 0;
					p[3][0] = 1;

					for(size_t i = 0; i < (unsigned)particles.size(); i++){
						matrix::fixed<4,1> x;

						matrix::prod( &particles[i].coordm,  &p, &x );

						lssmap::likelihood(_map, x[0][0], x[1][0], &lk);
						particles[i].likelihood += lk;
					}
				} // ---> loop for compute likelihood
				for(size_t i = 0; i < (unsigned)particles.size(); i++){
					sum += particles[i].likelihood;
					matrix::scalar_prod( &particles[i].pos, particles[i].likelihood, &ws3x1 );
					matrix::add(&delta, &ws3x1, &delta);
				}

				if(sum == 0)	return -1;
				// weighted average
				matrix::scalar_div(&delta, sum, &delta);

				LogDebugf("     : qmc - delta = (%lf, %lf, %lf):\n", delta[0][0], delta[1][0], delta[2][0]);

				// store delta
				matrix::copy(&_converge.delta, &delta);

				// output delta
				if( d )	matrix::copy(d, &delta);
				// compute
				matrix::add(&delta, &_v.pos, &_v.pos);

				// output position
				if( p ) matrix::copy(p, &_v.pos);


				{ // ---> compute distribution
					matrix::fixed<3,1> ws3x1;
					matrix::fixed<3,3> ws3x3;

					// ---> compute weighted summation
					set_zero(&_v.var);
					for(size_t i = 0; i < (unsigned)particles.size(); i++) {
						matrix::sub(&particles[i].pos, &delta, &ws3x1);
						matrix::prod_transpose2(&ws3x1, &ws3x1, &ws3x3);
						matrix::scalar_prod(&ws3x3, particles[i].likelihood, &ws3x3);
						matrix::add(&_v.var, &ws3x3, &_v.var);
					} // <--- compute weighted summation

					matrix::scalar_div(&_v.var, sum, &_v.var);
				} // <--- compute distribution

				if( l ) {// ---> set likelihood
					*l = particles[0].likelihood;
				} // <--- set likelihood


				{ // ---> refresh particle
					// store starting value
					particles.clear();
					create_particle(&_v);
				} /// <--- refresh particle
			} // <--- operate
			return 0;
		}


		inline
		optimize_quasi_monte_calro::initial_parameter::initial_parameter(){
			n = 1;
			var[0][0] = gnd_square( gnd_m2dist( 0.05 ) );
			var[1][1] = gnd_square( gnd_m2dist( 0.05 ) );
			var[2][2] = gnd_square( gnd_deg2ang( 3 ) );
		}
	}
};
// <--- class definition





// ---> class definition
namespace gnd {
	namespace lssmap {

		/**
		 * @brief probabilistic scan matching optimizer (quasi-monte-calro and newton's method)
		 */
		class optimize_hybrid_qmc2newton : public optimize_quasi_monte_calro {

			// ---> type declaration
		public:
			/// map type
			typedef lssmap::lssmap_t		map_t;
			/// map type pointer
			typedef map_t*			map_pt;
			/// pixel type
			typedef lssmap::lssmap_pixel_t	pixel_t;
			/// pixel type pointer
			typedef pixel_t*		pixel_pt;
			// <--- type declaration

			// ---> constructor, destructor
		public:
			optimize_hybrid_qmc2newton();
			optimize_hybrid_qmc2newton(map_pt m);
			~optimize_hybrid_qmc2newton();
			// <--- constructor, destructor

			// ---> optimization
		public:
			virtual int iterate(matrix::fixed<3,1> *d, matrix::fixed<3,1> *p, double *l);
			// <--- optimization
		};


		/**
		 * @brief constructor
		 */
		inline
		optimize_hybrid_qmc2newton::optimize_hybrid_qmc2newton() {
			return;
		}

		/**
		 * @brief constructor
		 * @param[in] m : map
		 */
		inline
		optimize_hybrid_qmc2newton::optimize_hybrid_qmc2newton(map_pt m)
		: optimize_quasi_monte_calro(m) {
			return;
		}

		/**
		 * @brief destructor
		 */
		inline
		optimize_hybrid_qmc2newton::~optimize_hybrid_qmc2newton(){
		}





		/**
		 * @brief optimization iterate
		 * @param[out] d : difference
		 * @param[out] p : pos
		 * @param[out] v : variance
		 * @param[out] l : likelihood
		 * @return    0 :
		 */
		inline
		int optimize_hybrid_qmc2newton::iterate(matrix::fixed<3,1> *d, matrix::fixed<3,1> *p, double *l){

			{ // ---> operate
				double lk, likelihood = 0;
				double sum;
				matrix::fixed<3,1> delta;
				matrix::fixed<3,1> ws3x1;

				// ---> quasi monte-calro method
				if( particles.size() > 0) {
					set_zero(&delta);
					sum = 0;

					LogVerbose("     : qmc2newton - quasi-monte-calro:\n");
					// ---> loop for compute likelihood
					for(size_t i = 0; i < (unsigned)particles.size(); i++){
						particles[i].likelihood = 0;
						for(size_t j = 0; j < (unsigned)_points.size(); j++ ){
							matrix::fixed<4,1> x;
							matrix::fixed<4,1> X;

							x[0][0] = _points[j][0][0];
							x[1][0] = _points[j][1][0];
							x[2][0] = 0;
							x[3][0] = 1;

							prod( &particles[i].coordm,  &x, &X );

							lssmap::likelihood(_map, X[0][0], X[1][0], &lk);
							particles[i].likelihood += lk;
						}
						sum += particles[i].likelihood;
						matrix::scalar_prod( &particles[i].pos, particles[i].likelihood, &ws3x1 );
						add(&delta, &ws3x1, &delta);
					} // ---> loop for compute likelihood

					likelihood = particles[0].likelihood;

					// exception likelihood = 0
					if(sum == 0)	return -1;
					// weighted average
					scalar_div( &delta, sum, &delta );

					{ // ---> compute distribution
						matrix::fixed<3,1> ws3x1;
						matrix::fixed<3,3> ws3x3;

						// ---> compute weighted summation
						matrix::set_zero(&_v.var);
						for(size_t i = 0; i < (unsigned)particles.size(); i++) {
							matrix::sub(&particles[i].pos, &delta, &ws3x1);
							matrix::prod_transpose2(&ws3x1, &ws3x1, &ws3x3);
							matrix::scalar_prod(&ws3x3, particles[i].likelihood, &ws3x3);
							matrix::add(&_v.var, &ws3x3, &_v.var);
						} // <--- compute weighted summation

						matrix::scalar_div(&_v.var, sum, &_v.var);
					} // <--- compute distribution

				} // <--- quasi monte-calro method
				// ---> newton's method
				else {
					int ret;
					uint64_t pi;
					matrix::fixed<3,1> grad;
					matrix::fixed<3,3> hess;
					matrix::fixed<4,4> coordm;

					LogVerbose("     : qmc2newton - newton:\n");
					// coordinate convert matrix
					matrix::coordinate_converter(&coordm,
							_v.pos[0], _v.pos[1], 0,
							::cos(_v.pos[2]), ::sin(_v.pos[2]), 0,
							 0, 0, 1);

					// ---> scanning loop of reflection points
					for( pi = 0; pi < _points.size(); pi++ ){
						matrix::fixed<3,1> g;
						matrix::fixed<3,3> h;

						// compute likelihood, gradient, hessian
						optimize_newton::_newton_method_variables_( _points[pi][0][0], _points[pi][1][0], &coordm, _map, &lk, &g, &h );

						// summation
						likelihood += lk;
						add( &grad, &g, &grad );
						add( &hess, &h, &hess );
					} // <--- scanning loop of reflection points

					// modified newton's method for unconstrained minimizer
					if( (ret = optimize::newtons_method_unconstrainted( &grad, &hess, &delta )) < 0)
						return ret;
				} // <--- newton's method

				LogDebugf("     : qmc2newton - delta = (%lf, %lf, %lf):\n", delta[0][0], delta[1][0], delta[2][0]);
				// merge delta
				matrix::add( &delta, &_v.pos, &_v.pos );
				// store delta
				matrix::copy( &_converge.delta, &delta );

				// output position
				if( p ) matrix::copy( p, &_v.pos );
				// set likelihood
				if( l ) *l = likelihood;
				// output delta
				if( d )	matrix::copy( d, &delta );

				// select qmc or neton's method
				if( !converge_test( gnd_square( _converge.delta[0][0] ) + gnd_square( _converge.delta[1][0] ) , _converge.delta[2][0],
						gnd_square( gnd_m2dist(0.001) ), gnd_deg2ang(1)) ) {
					// quasi monte calro method
					particles.clear();
					create_particle(&_v);
				}
				else {
					// newton's method
					particles.clear();
				}

			} // <--- operate
			return 0;
		}

	}
};



#endif /* GND_LSSMAP_OPTIMIZE_HPP_ */

