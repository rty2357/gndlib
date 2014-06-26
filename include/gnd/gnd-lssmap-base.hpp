/*
 * gnd-lssmap-base.hpp (GND's Laser Scan Statistics MAP BASE)
 *
 *  Created on: 2014/09/09
 *      Author: Taichi Yamada
 */

#ifndef GND_LSSMAP_BASE_HPP_
#define GND_LSSMAP_BASE_HPP_

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdint.h>

#include "gnd-gridmap.hpp"
#include "gnd-bmp.hpp"
#include "gnd-util.h"
#include "gnd-queue.hpp"
#include "gnd-matrix-base.hpp"
#include "gnd-vector-base.hpp"
#include "gnd-lib-error.h"



// include debug logging function
#define GND_DEBUG_LOG_NAMESPACE1 gnd
#define GND_DEBUG_LOG_NAMESPACE2 lssmap
#include "gnd-debug-log.hpp"
#undef GND_DEBUG_LOG_NAMESPACE2
#undef GND_DEBUG_LOG_NAMESPACE1
#undef GND_DEBUG_LOG

#include "gnd-debug-log-util-def.h"

/**
 * @ifnot GNDLSSMAP
 * @defgroup GNDLSSMAP laser scan statistics map
 * @endif
 */


// ---> constant and structure definition
namespace gnd {
	namespace lssmap {

		/**
		 * @brief number of map data plane
		 */
		static const size_t PlaneNum = 4;

		static const double SqrtDBL_EPSILON = ::sqrt(DBL_EPSILON);
		static const double ErrorMargin     = SqrtDBL_EPSILON;

		/**
		 * @brief default file out directory
		 */
		static const char CMapDirectoryDefault[] = ".";

		/**
		 * @brief counting map file name format
		 */
		static const char CMapFileNameFormat[] = "%s/%s.%02d.%s";

		/**
		 * @brief default counting map file name
		 */
		static const char CMapFileNameDefault[] = "lssmap";

		/**
		 * @brief default counting map file extension
		 */
		static const char CMapFileExtension[] = "cmap";

		/**
		 * @brief default map data cell size
		 */
		static const double DefaultMapCellSize = 5;

		/**
		 * @brief additional smoothing parameter
		 */
		static const uint32_t AddSmooth_SampleNum = 3;

		const double AddSmooth_SumCovXX = 2.0/3.0;
		const double AddSmooth_SumCovXY = 0.0;
		const double AddSmooth_SumCovYX = 0.0;
		const double AddSmooth_SumCovYY = 2.0/3.0;
		const double AddSmooth_Eta = 0.831220;



		/**
		 * @brief reflection point position index
		 */
		enum {
			PosX = 0,
			PosY = 1,
			PosDim = 2,
		};


		/**
		 * @brief statistics counting map
		 */
		struct counting_map_pixel_ {
			/// @brief sum of reflection point position value
			gnd::vector::fixed_column< PosDim > pos_sum;
			/// @brief sum of covariance
			matrix::fixed< PosDim,PosDim > cov_sum;
			/// @brief count of reflection point
			uint64_t cnt;
			/// @brief constructor
			counting_map_pixel_() : cnt(0) {}
		};

		/**
		 * @typedef counting_pixel_t
		 * @see counting_map_pixel
		 */
		typedef struct counting_map_pixel_ cmap_pixel_t;

		/**
		 * @ingroup GNDLSSMAP
		 * @brief counting map
		 */
		struct counting_map_ {
			/// @brief four planes
			gridmap::gridplane<cmap_pixel_t> plane[PlaneNum];
		};
		/**
		 * @ingroup GNDLSSMAP
		 * @typedef cmap_t
		 * @see counting_map
		 */
		typedef struct counting_map_ cmap_t;



		/**
		 * @brief pixel of laser scan statistics map
		 */
		struct lssmap_pixel_ {
			/// @brief mean
			gnd::vector::fixed_column< PosDim > mean;
			/// @brief inverse matrix of covariance
			gnd::matrix::fixed< PosDim,PosDim > inv_cov;
			/// @brief weight ( N / sqrt(|Sigma|)  ) * eta
			double K;
			/// @brief constructor
			lssmap_pixel_() : K(0.0) {}
		};
		/**
		 * @typedef lssmap_pixel_t
		 * @see map_pixel
		 */
		typedef struct lssmap_pixel_ lssmap_pixel_t;


		/**
		 * @ingroup GNDLSSMAP
		 * @brief map
		 */
		struct lssmap_ {
			/// @brief four planes
			gridmap::gridplane<lssmap_pixel_t> plane[PlaneNum];
		};
		/**
		 * @typedef lssmap_t
		 * @see lssmap
		 */
		typedef struct lssmap_ lssmap_t;

	}
};
// <--- constant and structure definition



// ---> function declaration
namespace gnd {
	namespace lssmap {

		int init_counting_map(cmap_t *m, double p, double u = DefaultMapCellSize);
		int clear_counting_map(cmap_t *m);
		int destroy_counting_map(cmap_t *m);

		int counting_map(cmap_t *m, double x, double y);

		int update_map(cmap_t *cnt, lssmap_t *map, double x, double y, double err = ErrorMargin );
		int update_ndt_map(cmap_t *c, lssmap_t *m, double x, double y);
		int build_map(lssmap_t *map, cmap_t *cnt, double sr = 20, uint32_t alpha = 100, double blur = ErrorMargin, uint64_t np = 10 );
		int build_ndt_map(lssmap_t *map, cmap_t *cnt, double err = ErrorMargin);
		int destroy_map(lssmap_t *m);

		int likelihood(lssmap_t *m, double x, double y, double *l);

		int read_counting_map(cmap_t *c,  const char* d = CMapDirectoryDefault, const char* f = CMapFileNameDefault, const char* e = CMapFileExtension);
		int write_counting_map(cmap_t *c,  const char* d = CMapDirectoryDefault, const char* f = CMapFileNameDefault, const char* e = CMapFileExtension);

		int build_bmp(bmp8_t *b, lssmap_t *m, double p = 0.1, double sr = 0.0, double cp = 4.0);
		int build_bmp(bmp32_t *b, lssmap_t *m, double p = 0.1, double sr = 0.0, double cp = 4.0);
		int build_bmp8(bmp8_t *b, lssmap_t *m, double p = 0.1, double sr = 0.0, double cp = 4.0);
		int build_bmp32(bmp32_t *b, lssmap_t *m, double p = 0.1, double sr = 0.0, double cp = 4.0);
		int likelihood(bmp8_t *m, double sx, double sy, double stheta, double r, double maxr, double *l);
		int likelihood(bmp32_t *m, double sx, double sy, double stheta, double r, double maxr, double *l);
	}
};
// <--- function declaration



// ---> function definition
namespace gnd {
	namespace lssmap {

		/**
		 * @ingroup GNDLSSMAP
		 * @brief initialize counting map
		 * @param[out] m : counting map
		 * @param[in]  p : pixel size
		 * @param[in]  u : grid map allocate unit size
		 */
		inline
		int init_counting_map(cmap_t *m, double p, double u)
		{
			gnd_assert(!m, -1, "invalid null pointer");
			gnd_assert(p <= 0, -1, "invalid argument. grid size must be greater than 0.");
			gnd_assert(u < p, -1, "invalid argument. cell size must be greater than 0.");

			LogDebugf("Begin - init_counting_map(%p, %lf, %lf)\n", m, p, u);
			LogIndent();

			// ---> plane scan loop
			for( size_t i = 0; i < PlaneNum; i++){
				// check error
				if( m->plane[i].is_allocate() ) {
					LogDebug("this map is buzy\n");
					LogUnindent();
					LogDebugf(" Fail - init_counting_map(%p, %lf, %lf)\n", m, p, u);
					gnd_exit(-1, "invalid argument, this map is buzy");
				}

				// allocate memory cell
				if( m->plane[i].pallocate(u, u, p, p) < 0) {
					LogDebug("fail to allocate\n");
					LogUnindent();
					LogDebugf("Fail  - init_counting_map(%p, %lf, %lf)\n", m, p, u);
					gnd_exit(-1, "fail to allocate");
				}
				else {
					LogVerbosef("map[%d] allocated\n", i);
				}
				// set origin
				m->plane[i].pset_core( (i % 2) * (p / 2.0), ((i / 2) % 2) * (p / 2.0) );
			} // <--- plane scan loop

			LogUnindent();
			LogDebugf("End   - init_counting_map(%p, %lf, %lf)\n", m, p, u);
			return 0;
		}


		/**
		 * @ingroup GNDLSSMAP
		 * @brief clear counting map
		 * @param[out] m :
		 */
		int clear_counting_map(cmap_t *m) {
			gnd_assert(!m, -1, "invalid null pointer");
			gnd_assert(!(m->plane+0) || !(m->plane + 1) || !(m->plane + 2) || !(m->plane + 3), -1, "invalid null pointer");

			{ // ---> operation
				gnd::lssmap::cmap_pixel_t ini;

				for( size_t i = 0; i < lssmap::PlaneNum; i++){
					m->plane[i].set_uniform(&ini);
				}
			} // <--- operation
			return 0;
		}

		/**
		 * @ingroup GNDLSSMAP
		 * @brief release counting map
		 * @param[out] m : counting map
		 */
		inline
		int destroy_counting_map(cmap_t *m)
		{
			gnd_assert(!m, -1, "invalid null pointer");

			LogDebugf("Begin - destroy_counting_map(%p)\n", m);
			LogIndent();

			for( size_t i = 0; i < PlaneNum; i++){
				m->plane[i].deallocate();
			}

			LogUnindent();
			LogDebugf("  End - destroy_counting_map(%p)\n", m);
			return 0;
		}

		/**
		 * @ingroup GNDLSSMAP
		 * @brief map counting function
		 * @param[out] m : counting map
		 * @param[in] x : laser scanner reflection point x
		 * @param[in] y : laser scanner reflection point y
		 */
		inline
		int counting_map(cmap_t *m, double x, double y)
		{
			gnd_assert(!m, -1, "invalid null argument");

			{ // ---> operate
				matrix::fixed<PosDim, 1> xx;

				set(&xx, PosX, 0, x);
				set(&xx, PosY, 0, y);

				// ---> plene scanning loop
				for(size_t i = 0; i < PlaneNum; i++){
					cmap_pixel_t *pp = 0;					// reference of pixel data
					gnd::vector::fixed_column<PosDim> 	core;	// pixel core position
					gnd::vector::fixed_column<PosDim> 	pos;	// reflection point on pixel
					matrix::fixed<PosDim, PosDim>	cov;	// reflection point covariance in pixel
					long r = 0, c = 0;						//

					// if memory is lacking, map data reallocate
					for( pp = m->plane[i].ppointer( x, y );
							pp == 0;
							pp = m->plane[i].ppointer( x, y ) ){
						m->plane[i].reallocate(x, y);
					}

					// get row and column number of reflection point pixel
					m->plane[i].pindex(x, y, &r, &c);
					// get core position of reflection point pixel
					m->plane[i].pget_pos_core(r, c, pointer(&core, PosX, 0), pointer(&core, PosY, 0));

					// compute reflection point on pixel
					matrix::sub(&xx, &core, &pos);

					// compute covariance
					matrix::prod_transpose2(&pos, &pos, &cov);

					// add date
					matrix::add(&pp->pos_sum, &pos, &pp->pos_sum);
					matrix::add(&pp->cov_sum, &cov, &pp->cov_sum);
					pp->cnt++;
				} // <--- plene scanning loop
			} // <--- operate
			return 0;
		}


		/**
		 * @ingroup GNDLSSMAP
		 * @brief map update
		 * @param[out] cnt : counting map
		 * @param[out]    map : map
		 * @param[in]       x : laser scanner reflection point x
		 * @param[in]       y : laser scanner reflection point y
		 * @param[in]     err : margin of some kinds of error ( such as rounding error sensor resolution  )
		 */
		inline
		int update_map(cmap_t *cnt, lssmap_t *map, double x, double y, double err )
		{
			gnd_assert(!map, -1, "invalid null argument");

			{ // ---> operation
				int ret;
				matrix::fixed<PosDim,PosDim> ws2x2;	// workspace 2x2 matrix
				matrix::fixed<PosDim,PosDim> cov;		// covariance matrix
				long r = 0, c = 0;

				if( (ret = counting_map(cnt, x, y)) < 0){
					return ret;
				}
				// ---> for each plane
				for(size_t i = 0; i < PlaneNum; i++){
					cmap_pixel_t *cpp;
					lssmap_pixel_t *pp;

					if( !map->plane[i].is_allocate() ) {
						map->plane[i].allocate( cnt->plane[i]._unit_row_(), cnt->plane[i]._unit_column_(),
								cnt->plane[i]._plane_row_(), cnt->plane[i]._plane_column_());
						map->plane[i].pset_origin( cnt->plane[i].xlower(), cnt->plane[i].ylower());
						map->plane[i].pset_rsl( cnt->plane[i].xrsl(), cnt->plane[i].yrsl());
					}

					// if memory is lacking, map data reallocate
					for( ret = map->plane[i].pindex(x, y, &r, &c);
							ret < 0;
							ret = map->plane[i].pindex(x, y, &r, &c) ){
						map->plane[i].reallocate(x, y);
						LogDebug("reallocate\n");
					}

					// get pointer
					if( !(pp = map->plane[i].ppointer(x,y)) ) continue;
					if( !(cpp = cnt->plane[i].ppointer(x,y)) ) continue;

//					{ // ---> obtain the number of points
//						pp->N = cpp->cnt;
//					} // <--- obtain the number of points

					// ---> obtain mean and inverse matrix of co-variance
					if(cpp->cnt <= 3){
						pp->K = 0;
						continue;
					}
					else {
						// compute mean
						matrix::scalar_div(&cpp->pos_sum, (double)cpp->cnt, &pp->mean );

						// compute covariance
						matrix::prod_transpose2(&pp->mean, &cpp->pos_sum, &ws2x2);
						matrix::sub(&cpp->cov_sum, &ws2x2, &ws2x2);
						matrix::scalar_div(&ws2x2, (double)cpp->cnt, &cov);

						// add minimal diagonal matrix
						matrix::set_unit(&ws2x2);
						matrix::scalar_prod(&ws2x2, gnd_square(err) , &ws2x2);
						matrix::add(&cov, &ws2x2, &cov);

						// compute inverse covariance
						matrix::inverse(&cov, &pp->inv_cov);
					}
					// --->obtain mean and inverse matrix of co-variance


					{ // ---> compute evaluation gain
						double det;

						// obtain determinant of co-variance matrix
						if( matrix::det(&cov, &det) < 0) {
							pp->K = 0;
							continue;
						}

						pp->K = (double) (cpp->cnt) / ( ::sqrt(det) );
					} // <--- compute evaluation gain

				} // <--- for each plane

			} // <--- operation

			return 0;
		}


		/**
		 * @privatesection
		 * @ingroup GNDLSSMAP
		 * @brief map counting function (ndt)
		 * @param[in,out] cnt : counting map
		 * @param[out]    map : map
		 * @param[in]       x : laser scanner reflection point x
		 * @param[in]       y : laser scanner reflection point y
		 * @param[in]     err : minimum error
		 */
		inline
		int update_ndt_map(cmap_t *cnt, lssmap_t *map, double x, double y, double err)
		{
			gnd_assert(!map, -1, "invalid null argument");

			{ // ---> operation
				int ret;
				matrix::fixed<PosDim,PosDim> ws2x2;	// workspace 2x2 matrix
				matrix::fixed<PosDim,PosDim> cov;		// covariance matrix
				long r, c;

				if( (ret = counting_map(cnt, x, y)) < 0) return ret;
				// ---> for each plane
				for(size_t i = 0; i < PlaneNum; i++){
					cmap_pixel_t *cpp;
					lssmap_pixel_t *pp;

					if( !map->plane[i].is_allocate() ) {
						map->plane[i].allocate( cnt->plane[i]._unit_row_(), cnt->plane[i]._unit_column_(),
								cnt->plane[i]._plane_row_(), cnt->plane[i]._plane_column_());
						map->plane[i].pset_origin( cnt->plane[i].xlower(), cnt->plane[i].ylower());
						map->plane[i].pset_rsl( cnt->plane[i].xrsl(), cnt->plane[i].yrsl());
					}

					// if memory is lacking, map data reallocate
					for( ret = map->plane[i].pindex(x, y, &r, &c);
							ret < 0;
							ret = map->plane[i].pindex(x, y, &r, &c) ){
						map->plane[i].reallocate(x, y);
					}

					// get pointer
					if( !(pp = map->plane[i].ppointer(x,y)) ) continue;
					if( !(cpp = cnt->plane[i].ppointer(x,y)) ) continue;

//					{ // ---> obtain the number of points
//						pp->N = 1;
//					} // <--- obtain the number of points

					// if 0
					if(cpp->cnt <= 3){
						pp->K = 0;
						continue;
					}
					else {
						// compute mean
						scalar_div(&cpp->pos_sum, (double)cpp->cnt, &pp->mean );

						// compute covariance
						prod_transpose2(&pp->mean, &cpp->pos_sum, &ws2x2);
						sub(&cpp->cov_sum, &ws2x2, &ws2x2);
						scalar_div(&ws2x2, (double)cpp->cnt, &cov);

						// add minimal diagonal matrix
						set_unit(&ws2x2);
						scalar_prod(&ws2x2, gnd_square(err) , &ws2x2);
						add(&cov, &ws2x2, &cov);

						// compute inverse covariance
						inverse(&cov, &pp->inv_cov);
					}

					{ // ---> compute evaluation gain
						pp->K = 1.0;
					} // <--- compute evaluation gain

				} // <--- for each plane

			} // <--- operation

			return 0;
		}


		/**
		 * @ingroup GNDLSSMAP
		 * @brief map building function
		 * @param[out]    map : builded map
		 * @param[in]     cnt : laser scanner reflection point counting data
		 * @param[in]      sr : sensor range
		 * @param[in]   alpha : additional smoothing parameter
		 * @param[in]    blur :
		 * @param[in]      np :
		 */
		inline
		int build_map(lssmap_t *map, cmap_t *cnt, double sr, uint32_t alpha, double blur, uint64_t np ) {
			gnd_assert(!cnt, -1, "invalid null pointer");
			gnd_assert(!map, -1, "map is null");

			LogDebugf("Begin - int build_map(%p, %p, %lf, %lf)\n", map, cnt, blur, np, alpha);
			LogIndent();

			{ // ---> operation
				matrix::fixed<PosDim,PosDim> ws2x2;	// workspace 2x2 matrix
				matrix::fixed<PosDim,PosDim> cov;	// covariance matrix

				// ---> for each plane
				for(size_t i = 0; i < PlaneNum; i++){
					LogVerbosef("plane %d ...\n", i);

					if( map->plane[i].is_allocate() ){
						if( map->plane[i].xrsl() != cnt->plane[i].xrsl() || map->plane[i].yrsl() != cnt->plane[i].yrsl() ){
							LogDebug("fail to memeory allocate\n");
							LogUnindent();
							LogDebugf("Fail  - int build_map(%p, %p, %lf, %lf)\n", map, cnt, blur, np, alpha );
							return -1;
						}
					}
					else {
						map->plane[i].allocate( cnt->plane[i]._unit_row_(), cnt->plane[i]._unit_column_(),
								cnt->plane[i]._plane_row_(), cnt->plane[i]._plane_column_());

						map->plane[i].pset_origin( cnt->plane[i].xlower(), cnt->plane[i].ylower());
						map->plane[i].pset_rsl( cnt->plane[i].xrsl(), cnt->plane[i].yrsl());
					}

					// ---> for each row
					for( unsigned long r = 0; r < cnt->plane[i].row(); r++){
						// ---> for each column
						for( unsigned long c = 0; c < cnt->plane[i].column(); c++){
							cmap_pixel_t cp;
							lssmap_pixel_t *pp;
							double x, y;
							uint64_t cnt_cell = 0;
							uint64_t sum;

							// get ndt data
							cnt->plane[i].get(r, c, &cp);
							cnt->plane[i].pget_pos_core(r, c, &x, &y);
							for( pp = map->plane[i].ppointer( x, y );
									pp == 0;
									pp = map->plane[i].ppointer( x, y ) ){
								map->plane[i].reallocate(x, y);
							}


							// ---> obtain mean and inverse matrix of co-variance
							{
								matrix::fixed<2, 2> ws2x2;


								ws2x2[0][0] = AddSmooth_SumCovXX * alpha / gnd_square(AddSmooth_SampleNum) * gnd_square(cnt->plane[i].xrsl());
								ws2x2[0][1] = AddSmooth_SumCovXY * alpha / gnd_square(AddSmooth_SampleNum) * cnt->plane[i].xrsl() * cnt->plane[i].yrsl();
								ws2x2[1][0] = AddSmooth_SumCovYX * alpha / gnd_square(AddSmooth_SampleNum) * cnt->plane[i].xrsl() * cnt->plane[i].yrsl();
								ws2x2[1][1] = AddSmooth_SumCovYY * alpha / gnd_square(AddSmooth_SampleNum) * gnd_square(cnt->plane[i].yrsl());

								cp.cnt += alpha;
								matrix::add(&cp.cov_sum, &ws2x2, &cp.cov_sum);

								// compute mean
								matrix::scalar_div(&cp.pos_sum, (double)cp.cnt, &pp->mean );

								// compute covariance
								matrix::prod_transpose2(&pp->mean, &cp.pos_sum, &ws2x2);
								matrix::sub(&cp.cov_sum, &ws2x2, &ws2x2);
								matrix::scalar_div(&ws2x2, (double)cp.cnt, &cov);

								// add minimal diagonal matrix
								matrix::set_unit(&ws2x2);
								matrix::scalar_prod(&ws2x2, gnd_square(blur) , &ws2x2);
								matrix::add(&cov, &ws2x2, &cov);

								// compute inverse covariance
								matrix::inverse(&cov, &pp->inv_cov);
							}
							// --->obtain mean and inverse matrix of co-variance


							{ // ---> compute evaluation gain
								double det = 0;
								double inv_normal = 0;


								// obtain determinant of co-variance matrix
								if( matrix::det(&cov, &det) < 0) {
									pp->K = 0;
									continue;
								}
								// ---> obtain the sum of laser points in the sensor range for each plane
								if( sr > 0 ) {
									size_t f = ::floor( sr / map->plane[0].xrsl());
									cmap_pixel_t *tmp_cpp;
									unsigned long lowerr;	// search range lower row
									unsigned long upperr;	// search range upper row
									unsigned long lowerc;	// search range lower column
									unsigned long upperc;	// search range upper column

									lowerr = r < f ? 0 : r - f;
									upperr = r + f >= map->plane[i].row() ? map->plane[i].row() : r + f;
									lowerc = c < f ? 0 : c - f;
									upperc = c + f >= map->plane[i].column() ? map->plane[i].column() : c + f;
									sum = 0;

									// compute average of local area number of observed point
									for( unsigned long rr = lowerr; rr < upperr; rr++){
										for( unsigned long cc = lowerc; cc < upperc; cc++){
											tmp_cpp = cnt->plane[i].pointer( rr, cc );
											sum += tmp_cpp->cnt;
											cnt_cell++;
										}
									}
								} // <--- obtain the sum of laser points in the sensor range for each plane
								else {
									sum = 1;
								}

								if( np > 0 ){ // --->  normalize factor
									// ---> in a case that the cell have enough data
									if(cp.cnt >= alpha + 3) {
										double k = 1.0 / (2.0 * M_PI * ( ::sqrt(det) ));
										uint32_t pcnt = 0;
										uint32_t ix, iy;
										gnd::vector::fixed_column<2> xx;
										gnd::vector::fixed_row<2> ws1x2;
										gnd::matrix::fixed<1,1> ws1x1;

										for(ix = 0; -(map->plane[i].xrsl() * ix / np) + pp->mean[0] >= - map->plane[i].xrsl() / 2; ix++){
											for(iy = 0; -(map->plane[i].yrsl() * iy / np) + pp->mean[1] >= - map->plane[i].yrsl() / 2; iy++){
												xx[0] = -(map->plane[i].xrsl() * ix / np) + pp->mean[0];
												xx[1] = -(map->plane[i].yrsl() * iy / np) + pp->mean[1];
												// difference from mean
												gnd::matrix::sub(&xx, &pp->mean, &xx);
												// compute likelihood
												gnd::matrix::prod_transpose1(&xx, &pp->inv_cov, &ws1x2);
												gnd::matrix::prod(&ws1x2, &xx, &ws1x1);
												ws1x1[0][0] = k * ::exp( - ws1x1[0][0] / 2.0) * (map->plane[i].xrsl() * map->plane[i].yrsl());
												inv_normal += ws1x1[0][0];
												pcnt++;
											}
										}

										for(ix = 1; -(map->plane[i].xrsl() * ix / np) + pp->mean[0] >= - map->plane[i].xrsl() / 2; ix++){
											for(iy = 1; (map->plane[i].yrsl() * iy / np) + pp->mean[1] < map->plane[i].yrsl() / 2; iy++){
												xx[0] = -(map->plane[i].xrsl() * ix / np) + pp->mean[0];
												xx[1] = (map->plane[i].yrsl() * iy / np) + pp->mean[1];
												// difference from mean
												gnd::matrix::sub(&xx, &pp->mean, &xx);
												// compute likelihood
												gnd::matrix::prod_transpose1(&xx, &pp->inv_cov, &ws1x2);
												gnd::matrix::prod(&ws1x2, &xx, &ws1x1);
												ws1x1[0][0] = k * ::exp( - ws1x1[0][0] / 2.0) * (map->plane[i].xrsl() * map->plane[i].yrsl());
												inv_normal += ws1x1[0][0];
												pcnt++;
											}
										}

										for(ix = 1; (map->plane[i].xrsl() * ix / np) + pp->mean[0] < map->plane[i].xrsl() / 2; ix++){
											for(iy = 1; -(map->plane[i].yrsl() * iy / np) + pp->mean[1] >= - map->plane[i].yrsl() / 2; iy++){
												xx[0] = (map->plane[i].xrsl() * ix / np) + pp->mean[0];
												xx[1] = -(map->plane[i].yrsl() * iy / np) + pp->mean[1];
												// difference from mean
												gnd::matrix::sub(&xx, &pp->mean, &xx);
												// compute likelihood
												gnd::matrix::prod_transpose1(&xx, &pp->inv_cov, &ws1x2);
												gnd::matrix::prod(&ws1x2, &xx, &ws1x1);
												ws1x1[0][0] = k * ::exp( - ws1x1[0][0] / 2.0) * (map->plane[i].xrsl() * map->plane[i].yrsl());
												inv_normal += ws1x1[0][0];
												pcnt++;
											}
										}

										for(ix = 1; (map->plane[i].xrsl() * ix / np) + pp->mean[0] < map->plane[i].xrsl() / 2; ix++){
											for(iy = 1; (map->plane[i].yrsl() * iy / np) + pp->mean[1] < map->plane[i].yrsl() / 2; iy++){
												xx[0] = (map->plane[i].xrsl() * ix / np) + pp->mean[0];
												xx[1] = (map->plane[i].yrsl() * iy / np) + pp->mean[1];
												// difference from mean
												gnd::matrix::sub(&xx, &pp->mean, &xx);
												// compute likelihood
												gnd::matrix::prod_transpose1(&xx, &pp->inv_cov, &ws1x2);
												gnd::matrix::prod(&ws1x2, &xx, &ws1x1);
												ws1x1[0][0] = k * ::exp( - ws1x1[0][0] / 2.0) * (map->plane[i].xrsl() * map->plane[i].yrsl());
												inv_normal += ws1x1[0][0];
												pcnt++;
											}
										}
										if( pcnt > 0)	inv_normal /= pcnt;
										else 			return -1;
									}	// <--- in a case that the cell have enough data
									else if(alpha > 0){
										inv_normal = AddSmooth_Eta;
									}
									else {
										inv_normal = 0;
									}
								} // <---  normalize factor

								// obtain determinant of co-variance matrix
								if( inv_normal <= 0) {
									pp->K = 0;
									continue;
								}

								pp->K = ( ((double) (cp.cnt)) / ((double)sum + cnt_cell * alpha)) / ( ::sqrt(det) );

							} // <--- compute evaluation gain

						} // <-- for each column
					} // <--- for each row

				} // <--- for each plane
			}  // <--- operation

			LogUnindent();
			LogDebugf("End   - int build_map(%p, %p, %lf, %lf)\n", map, cnt, blur, np);
			return 0;
		}


		/**
		 * @ingroup GNDPSM
		 * @brief map building function
		 * @param[out] map : builded map
		 * @param[in]  cnt : laser scanner reflection point counting data
		 * @param[in]  err : laser scanner field range
		 */
		inline
		int build_ndt_map( lssmap_t *map, cmap_t *cnt, double err )
		{
			gnd_assert(!cnt, -1, "invalid null pointer");
			gnd_assert(!map, -1, "map is null");

			LogDebugf("Begin - build_ndt_map(%p, %p, %lf)\n", map, cnt, err);
			LogIndent();

			{ // ---> operation
				matrix::fixed<PosDim,PosDim> ws2x2;	// workspace 2x2 matrix
				matrix::fixed<PosDim,PosDim> cov;		// covariance matrix

				// ---> for each plane
				for(size_t i = 0; i < PlaneNum; i++){
					LogVerbosef("plane %d ...\n", i);

					if( map->plane[i].is_allocate() ){
						if( map->plane[i].xrsl() != cnt->plane[i].xrsl() || map->plane[i].yrsl() != cnt->plane[i].yrsl() ){
							LogDebug("fail to allocate\n");
							LogUnindent();
							LogDebugf("Fail - build_ndt_map(%p, %p, %lf)\n", map, cnt, err);
							return -1;
						}
					}
					else {
						map->plane[i].allocate( cnt->plane[i]._unit_row_(), cnt->plane[i]._unit_column_(),
								cnt->plane[i]._plane_row_(), cnt->plane[i]._plane_column_());

						map->plane[i].pset_origin( cnt->plane[i].xlower(), cnt->plane[i].ylower());
						map->plane[i].pset_rsl( cnt->plane[i].xrsl(), cnt->plane[i].yrsl());
					}

					// ---> for each row
					for( unsigned long r = 0; r < cnt->plane[i].row(); r++){
						// ---> for each column
						for( unsigned long c = 0; c < cnt->plane[i].column(); c++){
							cmap_pixel_t *cpp;
							lssmap_pixel_t *pp;
							double x, y;


							// get ndt data
							cpp = cnt->plane[i].pointer( r, c );
							cnt->plane[i].pget_pos_core(r, c, &x, &y);
							for( pp = map->plane[i].ppointer( x, y );
									pp == 0;
									pp = map->plane[i].ppointer( x, y ) ){
								map->plane[i].reallocate(x, y);
							}

							//					{ // ---> obtain the number of points
							//						pp->N = 1;
							//					} // <--- obtain the number of points

							// if 0
							if(cpp->cnt <= 3){
								pp->K = 0;
								continue;
							}
							else {
								// compute mean
								scalar_div(&cpp->pos_sum, (double)cpp->cnt, &pp->mean );

								// compute covariance
								prod_transpose2(&pp->mean, &cpp->pos_sum, &ws2x2);
								sub(&cpp->cov_sum, &ws2x2, &ws2x2);
								scalar_div(&ws2x2, (double)cpp->cnt, &cov);

								// add minimal diagonal matrix
								set_unit(&ws2x2);
								scalar_prod(&ws2x2, gnd_square(err) , &ws2x2);
								add(&cov, &ws2x2, &cov);

								// compute inverse covariance
								inverse(&cov, &pp->inv_cov);
							}

							pp->K = 1.0;

						} // <-- for each column
					} // <--- for each row

				} // <--- for each plane
			}  // <--- operation

			LogUnindent();
			LogDebugf("End - build_ndt_map(%p, %p, %lf)\n", map, cnt, err);
			return 0;
		}


		/**
		 * @ingroup GNDPSM
		 * @brief destory map
		 * @param[out] m :  map
		 */
		inline
		int destroy_map(lssmap_t *m)
		{
			gnd_assert(!m, -1, "invalid null pointer");

			LogDebugf("Begin - destroy_map(%p)\n", m);
			LogIndent();

			for( size_t i = 0; i < PlaneNum; i++)
				m->plane[i].deallocate();

			LogUnindent();
			LogDebugf("END   - destroy_map(%p)\n", m);
			return 0;
		}


		/**
		 * @ingroup GNDPSM
		 * @brief compute likelihood
		 * @param[in] map : map data
		 * @param[in]   x : laser scanner reflection point x
		 * @param[in]   y : laser scanner reflection point y
		 * @param[in]  pg : position gain
		 * @param[out]  l : likelihood
		 */
		inline
		int likelihood(lssmap_t *map, double x, double y, double *l){
			gnd_assert(!l, -1, "invalid null pointer");
			gnd_assert(!map, -1, "map is null");

			{ // ---> operate
				matrix::fixed<PosDim,1> pos;
				matrix::fixed<PosDim,1> q;
				matrix::fixed<1,PosDim> ws1x2;	// workspace 2x2 matrix
				matrix::fixed<1,1> ws1x1;		// workspace 1x1 matrix

				{ // ---> initialize
					set(&pos, PosX, 0, x);
					set(&pos, PosY, 0, y);
					*l = 0;
				} // <--- initialize


				// ---> for
				for( size_t i = 0; i < PlaneNum; i++){
					lssmap_pixel_t *pp = 0;
					int ret;
					long pr, pc;

					// get index of pixel
					ret = map->plane[i].pindex( value(&pos, PosX, 0), value(&pos, PosY, 0), &pr, &pc );
					// no data
					if( ret < 0 )			continue;

					// get pixel data
					pp = map->plane[i].pointer( pr, pc );
					// zero weight
					if( pp->K <= 0.0 )	continue;

					// get pixel core pos on ndt data pixel
					map->plane[i].pget_pos_core(pr, pc, pointer(&q, PosX, 0), pointer(&q, PosY, 0));
					matrix::sub(&pos, &q, &q);

					// difference from mean
					matrix::sub(&q, &pp->mean, &q);
					// compute likelihood
					matrix::prod_transpose1(&q, &pp->inv_cov, &ws1x2);
					matrix::prod(&ws1x2, &q, &ws1x1);

					*l += (pp->K) * ::exp( - value(&ws1x1, 0, 0) / 2.0);
				} // ---> for

				// normalization
				*l /= PlaneNum;
				return 0;
			} // <--- operate
		}

		/**
		 * @ingroup GNDPSM
		 * @brief compute likelihood
		 * @param[in] map : map data
		 * @param[in]     sx : sensor position x
		 * @param[in]     sy : sensor position y
		 * @param[in] stheta : sensor direction
		 * @param[in]      r : sensor reading (range)
		 * @param[in]   maxr : sensor range
		 * @param[out]     l : likelihood
		 */
		inline
		int likelihood(bmp32_t *m, double sx, double sy, double stheta, double r, double maxr, double *l) {
			gnd_assert(!l, -1, "invalid null pointer");
			gnd_assert(!m, -1, "map is null");

			{ // ---> operate
				double sinv, cosv;
				double ur = m->xrsl() < m->yrsl() ? m->xrsl() : m->yrsl();

				{ // ---> initialize
					sinv = ::sin(stheta);
					cosv = ::cos(stheta);

					*l = 0;
				} // <--- initialize

				{ // ---> compute likelihood
					double ieta = 0;
					double xx, yy;
					double rr;
					uint32_t v;

					if( m->pget(sx + r * cosv, sy + r * sinv, &v) < 0 ) {
						*l = 0;
						return -1;
					}

					// normalize factor
					for( rr = 0; rr < maxr; rr += ur ){
						uint32_t vv;
						xx = sx + rr * cosv;
						yy = sy + rr * sinv;
						if( m->pget(xx, yy, &vv) < 0 ) {
							vv = 0;
						}
						ieta += vv;
					}
					*l = (double) v / ieta;

				} // <--- compute likelihood
				return 0;
			} // <--- operate
		}

		/**
		 * @ingroup GNDPSM
		 * @brief compute likelihood
		 * @param[in] map : map data
		 * @param[in]     sx : sensor position x
		 * @param[in]     sy : sensor position y
		 * @param[in] stheta : sensor direction
		 * @param[in]      r : sensor reading (range)
		 * @param[in]   maxr : sensor range
		 * @param[out]     l : likelihood
		 */
		inline
		int likelihood(bmp8_t *m, double sx, double sy, double stheta, double r, double maxr, double *l) {
			gnd_assert(!l, -1, "invalid null pointer");
			gnd_assert(!m, -1, "map is null");

			{ // ---> operate
				double sinv, cosv;
				double ur = m->xrsl() / 2;

				{ // ---> initialize
					sinv = ::sin(stheta);
					cosv = ::cos(stheta);

					*l = 0;
				} // <--- initialize

				{ // ---> compute likelihood
					double ieta = 0;
					double xx, yy;
					double rr;
					uint8_t v;

					if( m->pget(sx + r * cosv, sy + r * sinv, &v) < 0 ) {
						*l = 0;
						return -1;
					}

					// normalize factor
					for( rr = 0; rr < maxr; rr += ur ){
						uint8_t vv;
						xx = sx + rr * cosv;
						yy = sy + rr * sinv;
						if( m->pget(xx, yy, &vv) < 0 ) {
							vv = 0;
						}
						ieta += vv;
					}

					*l = (double) v / ieta;

				} // <--- compute likelihood
				return 0;
			} // <--- operate
		}




		/**
		 * @ingroup GNDPSM
		 * @brief counting map file read
		 * @param[out] c : counting map
		 * @param[in] d : directory path
		 * @param[in] f : file name template
		 * @param[in] e : extention
		 */
		inline
		int read_counting_map(cmap_t *c,  const char* d, const char* f, const char* e)
		{
			gnd_assert(!c, -1, "invalid null argument");
			gnd_assert(!d, -1, "invalid null argument");
			gnd_assert(!f, -1, "invalid null argument");
			gnd_assert(!e, -1, "invalid null argument");

			LogDebugf("Begin - read_counting_map(%p, %p, %lf, %lf, %lf)\n", c, d, f, e);
			LogIndent();

			{ // ---> operation
				char path[1024];
				// ---> map plane data scanning loop
				for( size_t i = 0; i < PlaneNum; i++){
					::sprintf(path, CMapFileNameFormat, d, f, i, e);
					LogDebugf("file #%d path \"%s\"\n", i, path);
					if( c->plane[i].fread(path) < 0 ) {
						LogUnindent();
						LogDebugf("Fail  - read_counting_map(%p, %p, %lf, %lf, %lf)\n", c, d, f, e);
						gnd_exit(-1, "fail to file-open");
					}
				} // <--- map plane data scanning loop
			} // <--- operation
			LogUnindent();
			LogDebugf("End  - read_counting_map(%p, %p, %lf, %lf, %lf)\n", c, d, f, e);
			return 0;
		}

		/**
		 * @ingroup GNDPSM
		 * @brief counting map file out
		 * @param[in] c : counting map
		 * @param[in] d : directory path
		 * @param[in] f : file name template
		 * @param[in] e : extention
		 */
		inline
		int write_counting_map(cmap_t *c,  const char* d, const char* f, const char* e)
		{
			gnd_assert(!c, -1, "invalid null argument");
			gnd_assert(!d, -1, "invalid null argument");
			gnd_assert(!f, -1, "invalid null argument");
			gnd_assert(!e, -1, "invalid null argument");

			LogDebugf("Begin - write_counting_map(%p, %p, %lf, %lf, %lf)\n", c, d, f, e);
			LogIndent();

			{ // ---> operation
				char path[1024];
				// ---> map plane data scanning loop
				for( size_t i = 0; i < PlaneNum; i++){
					::sprintf(path, CMapFileNameFormat, d, f, i, e);
					LogVerbosef("file #%d path \"%s\"\n", i, path);
					if( c->plane[i].fwrite(path) < 0 ) {
						LogUnindent();
						LogDebugf("Fail - write_counting_map(%p, %p, %lf, %lf, %lf)\n", c, d, f, e);
						gnd_exit(-1, "fail to file-open");
					}
				} // <--- map plane data scanning loop
			} // <--- operation
			LogUnindent();
			LogDebugf("End  - write_counting_map(%p, %p, %lf, %lf, %lf)\n", c, d, f, e);
			return 0;
		}


		/**
		 * @brief build bitmap data (gray)
		 * @param[out] bmp : gray scale bit map
		 * @param[in] map : map data
		 * @param[in] ps : pixel size
		 * @param[in] sr : smoothing parameter ( sensor range[m])
		 * @param[in] cp : contranst parameter
		 */
		inline
		int build_bmp(bmp8_t *bmp, lssmap_t *map, double ps, double sr, double cp) {
			return build_bmp8(bmp, map, ps, sr, cp);
		}

		/**
		 * @brief build bitmap data (gray)
		 * @param[out] bmp : gray scale bit map
		 * @param[in] map : map data
		 * @param[in] ps : pixel size
		 * @param[in] sr : sensor range[m]
		 * @param[in] cp : contranst parameter
		 */
		inline
		int build_bmp(bmp32_t *bmp, lssmap_t *map, double ps, double sr, double cp) {
			return build_bmp32(bmp, map, ps, sr, cp);
		}


		/**
		 * @brief build bitmap data (gray)
		 * @param[out] bmp : gray scale bit map
		 * @param[in] map  : map data
		 * @param[in] ps   : pixel size
		 * @param[in] sr   : sensor range (for smoothing)
		 * @param[in] cp   : contrast parameter
		 */
		inline
		int build_bmp8(bmp8_t *bmp, lssmap_t *map, double ps, double sr, double cp)
		{
			gridmap::gridplane<double> ws;	// workspace

			gnd_assert(!bmp, -1, "invalid null argument");
			gnd_assert(!map, -1, "invalid null argument");
			gnd_assert(ps <= 0, -1, "invalid argument. pixel size must be greater than 0.");

			{ // ---> initialize
				// allocate bmp data
				if(bmp->is_allocate())	bmp->deallocate();
				bmp->pallocate(map->plane[0].xupper() - map->plane[3].xlower(), map->plane[0].yupper() - map->plane[3].ylower(), ps, ps);
				bmp->pset_origin(map->plane[3].xlower(), map->plane[3].ylower());

				ws.pallocate(map->plane[0].xupper() - map->plane[3].xlower(), map->plane[0].yupper() - map->plane[3].ylower(), ps, ps);
				ws.pset_origin(map->plane[3].xlower(), map->plane[3].ylower());
			} // <--- initialize

			{ // ---> operation
				double lkh = 0;		// likelihood
				double x, y;		// bitmap pixel core position
				double mu;
				double max = 0;
				double sigma;
				unsigned long cnt = 0;

				LogVerbose("compute likelihood");

				{ // ---> compute likelihood and average
					mu = 0;

					for( unsigned long r = 0; r < ws.row(); r++){
						for( unsigned long c = 0; c < ws.column(); c++){

							// get pixel core position
							ws.pget_pos_core(r, c, &x, &y);

							// compute likelihood
							likelihood(map, x, y, &lkh);
							lkh *= (bmp->xrsl() * bmp->yrsl());

							ws.set(r, c, &lkh);

							max = lkh > max ? lkh : max;
							if( lkh != 0 ) {
								mu += lkh;
								cnt++;
							}
						}
					}
					if( cnt == 0 ) return -1;
					mu /= cnt;
				} // ---> compute likelihood and average

				{ // ---> compute likelihood and average
					sigma = 0;
					for( unsigned long r = 0; r < ws.row(); r++){
						for( unsigned long c = 0; c < ws.column(); c++){
							ws.get(r, c, &lkh);
							if( lkh == 0 ) continue;
							sigma += gnd_square(lkh - mu);
						}
					}
					sigma /= cnt;
					sigma = ::sqrt(sigma);
				} // ---> compute likelihood and average

				{ // ---> set bitmap
					unsigned char bpv;	// bitmap pixel value

					for( unsigned long r = 0; r < bmp->row(); r++){
						// ---> grid map scanning loop column
						for( unsigned long c = 0; c < bmp->column(); c++){
							ws.get(r, c, &lkh);

							// normalize for bitmap
							if( lkh <= 0 ){
								bpv = 0;
							}
							else {
								lkh = (lkh) / (max);
								lkh *= 0xff;
								bpv = lkh > 0xff ? 0xff : (unsigned char)lkh;
							}

							bmp->set( r, c, &bpv);
						} // <--- grid map scanning loop column
					} // <--- grid map scanning loop
				} // <--- set bitmap

			} // <--- operation

			{ // ---> finalize
				ws.deallocate();
			} // <--- finalize

			return 0;
		}

		/**
		 * @brief build bitmap data (gray)
		 * @param[out] bmp : gray scale bit map
		 * @param[in] map  : map data
		 * @param[in] ps   : pixel size
		 * @param[in] sr   : sensor range (for smoothing)
		 * @param[in] cp   : contrast parameter
		 */
		inline
		int build_bmp32(bmp32_t *bmp, lssmap_t *map, double ps, double sr, double cp)
		{
			gridmap::gridplane<double> ws;	// workspace

			gnd_assert(!bmp, -1, "invalid null argument");
			gnd_assert(!map, -1, "invalid null argument");
			gnd_assert(ps <= 0, -1, "invalid argument. pixel size must be greater than 0.");


			{ // ---> initialize
				// allocate bmp data
				if(bmp->is_allocate())	bmp->deallocate();
				bmp->pallocate(map->plane[0].xupper() - map->plane[3].xlower(), map->plane[0].yupper() - map->plane[3].ylower(), ps, ps);
				bmp->pset_origin(map->plane[3].xlower(), map->plane[3].ylower());

				ws.pallocate(map->plane[0].xupper() - map->plane[3].xlower(), map->plane[0].yupper() - map->plane[3].ylower(), ps, ps);
				ws.pset_origin(map->plane[3].xlower(), map->plane[3].ylower());
			} // <--- initialize

			{ // ---> operation
				double lkh = 0;		// likelihood
				double x, y;		// bitmap pixel core position
				double mu;
				double sigma;
				double lkh_max;
				unsigned long cnt = 0;

				LogVerbose("compute likelihood");

				{ // ---> compute likelihood and average
					mu = 0;
					lkh_max = 0;

					for( unsigned long r = 0; r < ws.row(); r++){
						for( unsigned long c = 0; c < ws.column(); c++){

							// get pixel core position
							ws.pget_pos_core(r, c, &x, &y);

							// compute likelihood
							likelihood(map, x, y, &lkh);
							lkh *= (bmp->xrsl() * bmp->yrsl());

							ws.set(r, c, &lkh);

							if( lkh != 0 ) {
								mu += lkh;
								cnt++;
							}
							lkh_max = lkh_max < lkh ? lkh : lkh_max;
						}
					}
					if( cnt == 0 ) return -1;
					mu /= cnt;
				} // ---> compute likelihood and average

				{ // ---> compute likelihood and average
					sigma = 0;
					for( unsigned long r = 0; r < ws.row(); r++){
						for( unsigned long c = 0; c < ws.column(); c++){
							ws.get(r, c, &lkh);
							if( lkh == 0 ) continue;
							sigma += gnd_square(lkh - mu);
						}
					}
					sigma /= cnt;
					sigma = ::sqrt(sigma);
				} // ---> compute likelihood and average

				{ // ---> set bitmap
					uint32_t bpv;	// bitmap pixel value
					double max = lkh_max;


					for( unsigned long r = 0; r < bmp->row(); r++){
						// ---> grid map scanning loop column
						for( unsigned long c = 0; c < bmp->column(); c++){
							ws.get(r, c, &lkh);

							// normalize for bitmap
							if( lkh <= 0 ){
								bpv = 0;
							}
							else {
								lkh = (lkh) / (max);
								lkh *= 0xffffff;
								bpv = lkh > (0xffffff) ? (0xffffff) : (uint32_t)lkh;
							}

							bmp->set( r, c, &bpv);
						} // <--- grid map scanning loop column
					} // <--- grid map scanning loop
				} // <--- set bitmap

			} // <--- operation

			{ // ---> finalize
				ws.deallocate();
			} // <--- finalize

			return 0;
		}
	}
};
// <--- function definition




#endif /* GND_LSSMAP_BASE_HPP_ */
