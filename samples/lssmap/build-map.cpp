/*
 * sample.cpp
 *
 *  Created on: 2013/06/02
 *      Author: tyamada
 */

#include <stdio.h>
#include "gnd-lssmap-base.hpp"
#include "float.h"

int main(int argc, char* argv[]) {
	FILE						*fp;		// data file of laser scanner
	gnd::lssmap::cmap_t		cnt_map;	// laser scanner data collection map (it use for building environment map)


	gnd::lssmap::debug_set_log_level(2);
	gnd::gridmap::debug_set_log_level(2);

	{ // ---> initialize
		if( argc < 2 ) {
			::fprintf(stderr, "error: missing data file operand\n");
			::fprintf(stderr, "samples$ ./%s laser.dat\n", argv[0]);
			return 1;
		}
		else if( !(fp = ::fopen(argv[1], "r")) ){
			::fprintf(stderr, "fail to open %s\n", argv[1]);
			return 1;
		}

		gnd::lssmap::init_counting_map(&cnt_map, 1.0, 10 );
		//		gnd::lssmap::read_counting_map(&cnt_map, "lssmap-map1");
	} // <--- initialize


	{ // ---> operation

		{ // ---> build map
			double x, y;
			int cnt = 0;
			int ret;

			// ---> scanning loop (data fiile)
			while( !::feof(fp)  ) {
				// read data file
				ret = ::fscanf(fp, "%lf %lf\n", &x, &y);
                          fprintf(stdout, "%lf %lf\n", x, y);
				// when read all data, break
				if( ret <= 0 ) break;
				if( cnt % 1000 )
					::fprintf(stdout, ".");
				cnt++;
				// counting
				gnd::lssmap::counting_map(&cnt_map, x, y);
			} // ---> scanning loop (data fiile)
			::fprintf(stdout, "\n");


		} // <--- build map


	} // <--- operation



	{ // ---> finalize
		gnd::bmp8_t bmp;
		gnd::bmp32_t bmp32;
		gnd::lssmap::lssmap_t			lssmap;	// environment map for scan matching

		// save counting map
		// output in current directory
		gnd::lssmap::write_counting_map(&cnt_map, "./");
		::fprintf(stdout, "make counting map data observ-prob.**.cmap\n");


		{ // ---> build bmp image (to visualize for human)
			// build environmental map
			gnd::lssmap::build_map(&lssmap, &cnt_map, 20, 100);
			::fprintf(stdout, "make bmp image\n");
			// make bmp image: it show the likelihood field
			gnd::lssmap::build_bmp(&bmp, &lssmap, 1.0 / 20);
			// make bmp image: it show the likelihood field
			gnd::lssmap::build_bmp(&bmp32, &lssmap, 1.0 / 20);
			// file out
			gnd::bmp::write8("map-image8.bmp", &bmp);
			// file out
			gnd::bmp::write32("map-image32.bmp", &bmp32);

			gnd::lssmap::destroy_map(&lssmap);
			gnd::lssmap::destroy_counting_map(&cnt_map);
			{ // ---> origin
				char fname[512];
				FILE *fp = 0;
				double x, y;

				if( ::snprintf(fname, sizeof(fname), "%s", "origin.txt"  ) == sizeof(fname) ){
					::fprintf(stderr, "  ... \x1b[1m\x1b[31mError\x1b[39m\x1b[0m: fail to open. file name is too long\n");
				}
				else if( !(fp = fopen(fname, "w")) ) {
					::fprintf(stderr, "  ... \x1b[1m\x1b[31mError\x1b[39m\x1b[0m: fail to open \"\x1b[4m%s\x1b[0m\"\n", fname);
				}
				bmp.pget_origin(&x, &y);
				fprintf(fp, "%lf %lf\n", x, y);
				fclose(fp);
			} // --->  origin


			::fprintf(stdout, "make map image %s\n", "map-image.bmp");
			bmp.deallocate();
			bmp32.deallocate();
		} // <--- build bmp image (to visualize for human)


	} // <--- finalize

	return 0;
}




