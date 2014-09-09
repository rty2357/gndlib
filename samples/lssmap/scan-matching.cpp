/*
 * scan-matching.cpp
 *
 *  Created on: 2013/06/08
 *      Author: tyamada
 */


#include <stdio.h>
#include "gnd-lssmap.hpp"
#include "gnd-config-file.hpp"

#ifdef __linux__
#include "gnd-timer.hpp"

#define TIMER_INIT() 			gnd::stopwatch _sw_
#define TIMER_BEGIN() 			_sw_.begin( CLOCK_REALTIME )
#define TIMER_REC(t,l)										\
		do {													\
			_sw_.get( 0, (l) );									\
			*t += *l;											\
		} while(0)

#define TIMER_SHOW(t, l)									\
		do {													\
			::fprintf(stdout, "   lap time = %lf\n", *(l));		\
			::fprintf(stdout, " total time = %lf\n", *(t));		\
		} while(0)

#else

#define TIMER_INIT()
#define TIMER_BEGIN()
#define TIMER_REC(t, l)
#define TIMER_SHOW(t, l)

#endif

struct particle {
	double x;
	double y;
	double th;
	double lh;	// likelihood
} ;

struct position {
	double x;
	double y;
};

int main( int argc, char* argv[] ) {
	gnd::lssmap::cmap_t				cmap;					// counting map
	gnd::lssmap::lssmap_t			lssmap;					// scan matching map
	gnd::bmp32_t					bmp;					// bmp map (likelihood field)
	struct position					*scandata;				// scan data
	int								nscandata;
	struct particle 				*particles;
	int max;
	double max_lh;


	gnd::conf::parameter_array<char, 256> cmap_dir = {
			"cmap-dir", 							// item name
			"./", 									// default value
			"cmap directory path" 					// comment
	};
	gnd::conf::parameter_array<double, 3> pos_ini = {
			"init-pos", 							// item name
			{0, 0, 0}, 								// default value
			"initial position for scan matching (x[m], y[m], theta[deg])" 	// comment
	};
	gnd::conf::parameter<int> nparticle = {
			"number-of-particle", 							// item name
			1000, 					// default value
			"number of particles" 	// comment
	};
	gnd::conf::parameter_array<double, 3> error_dist = {
			"error-distribution", 							// item name
			{2.0, 2.0, 30 * M_PI / 180.0}, 					// default value
			"error distribution" 	// comment
	};
	gnd::conf::parameter_array<char, 256> scan_file = {
			"scan-data",							// item name
			"scan.txt",						// default value
			"scan data file path"					// comment
	};



	{ // ---> debug mode
		gnd::lssmap::debug_set_log_level(2);
		gnd::gridmap::debug_set_log_level(2);
	} // <--- debug mode

	{ // ---> initialize
		int ret;
		gnd::conf::file_stream fs;

		// set parameter item and default parameter
		if( argc < 2 ) {
			static const char fname[] = "_scan-matching.conf";
			::fprintf(stderr, " error: missing file operand\n");
			::fprintf(stderr, "        please input configuration file path\n");
			::fprintf(stderr, "        e.g.: %s scan-matching.conf\n", argv[0]);
			::fprintf(stderr, "        ... \n");
			// write configuration file
			gnd::conf::set_parameter(&fs, &cmap_dir);
			pos_ini.value[2] = pos_ini.value[2] * 180 / M_PI;
			gnd::conf::set_parameter(&fs, &pos_ini);
			error_dist.value[2] = error_dist.value[2] * 180 / M_PI;
			gnd::conf::set_parameter(&fs, &error_dist);
			gnd::conf::set_parameter(&fs, &nparticle);
			gnd::conf::set_parameter(&fs, &scan_file);
			fs.write(fname);
			::fprintf(stderr, "         => generate configuration file template \"%s\"\n", fname);

			return -1;
		}

		// configuration file read
		ret = fs.read(argv[1]);
		if( ret < 0 ) {
			::fprintf(stderr, " error: fail to read configuration file\n");
			return -1;
		}

		// get parameter value
		gnd::conf::get_parameter(&fs, &cmap_dir);
		if( gnd::conf::get_parameter(&fs, &pos_ini) )		pos_ini.value[2] = pos_ini.value[2] * M_PI / 180;
		if( gnd::conf::get_parameter(&fs, &error_dist) )	error_dist.value[2] = error_dist.value[2] * M_PI / 180;
		gnd::conf::get_parameter(&fs, &nparticle);
		gnd::conf::get_parameter(&fs, &scan_file);

		if( !cmap_dir.value[0] ) {
			::fprintf(stdout, " error: missing counting map\n");
			::fprintf(stdout, " please set counting map directory\n");
			return -1;
		}

		// read counting map
		ret = gnd::lssmap::read_counting_map( &cmap, cmap_dir.value );
		if( ret < 0 ) {
			::fprintf(stderr, " error: fail to read counting map form directory \"%s\"\n", cmap_dir.value);
			return 0;
		}

		{ // ---> read scan data file
			FILE	*fp;					// matched scan data (input)
			FILE	*fp_ini;				// scan data on initial position (file out)
			int		cnt = 0;
			double	cosv = cos(pos_ini.value[2]);
			double	sinv = sin(pos_ini.value[2]);


			// scan data file open
			if( !scan_file.value[0] ) {
				::fprintf(stdout, " error: missing scan data\n");
				::fprintf(stdout, "        please set scan data file\n");
				return -1;
			}
			fp = ::fopen( scan_file.value, "r" );
			if( !fp ) {
				::fprintf(stdout, " error: fail to open \"%s\"\n", scan_file.value);
				return -1;
			}

			fp_ini = ::fopen("scan-on-init-pos.txt", "w");
			if( !fp_ini ) {
				::fprintf(stdout, " error: fail to open \"%s\"\n", "scan-on-init-pos.txt");
				return -1;
			}

			// first file read (get number of scan data)
			cnt = 0;
			while( !::feof(fp)  ) {
				double x, y;
				// read data file
				ret = ::fscanf(fp, "%lf %lf\n", &x, &y);
				// when read all data, break
				if( ret < 0 ) break;
				cnt++;
			}

			// allocate
			nscandata = cnt;
			scandata = new struct position[nscandata];

			::fseek(fp, 0, SEEK_SET);
			// first file read (get number of scan data)
			cnt = 0;
			while( !::feof(fp)  ) {
				double x, y;
				// read data file
				ret = ::fscanf(fp, "%lf %lf\n", &x, &y);
				// when read all data, break
				if( ret < 0 ) break;
				scandata[cnt].x = x;
				scandata[cnt].y = y;
				// file out scan data on initial position
				::fprintf(fp_ini, "%lf %lf\n",
						scandata[cnt].x * cosv - scandata[cnt].y * sinv + pos_ini.value[0],
						scandata[cnt].x * sinv + scandata[cnt].y * cosv + pos_ini.value[1]);
				cnt++;
			}
			::fclose(fp);
			::fclose(fp_ini);

		} // <--- read scan data file

		// build scan matching map
		ret = gnd::lssmap::build_map( &lssmap, &cmap, 30.0, 10 );
		if( ret < 0 ) {
			::fprintf(stderr, " error: fail to build map\n");
			return 0;
		}

		// build likelihood filed
		ret = gnd::lssmap::build_bmp( &bmp, &lssmap, 0.05 );
		if( ret < 0 ) {
			::fprintf(stderr, " error: fail to build likelihood firled\n");
			return 0;
		}

		// allocate
		particles = new struct particle[ nparticle.value ];

	} // <--- initialize



	{ // ---> operation
		double t = 0, l = 0;

		// timer init (for linux)
		TIMER_INIT();


		// set random seed
		gnd::random_set_seed();

		// show initial position
		::fprintf(stdout, "       init pos = %lf, %lf, %lf\n", pos_ini.value[0],  pos_ini.value[1],  pos_ini.value[2] * 180 / M_PI);
		::fprintf(stdout, "   particle num = %d\n",  nparticle.value);
		::fprintf(stdout, "  scan data num = %d\n",  nscandata);

		::fprintf(stdout, "\n");
		::fprintf(stdout, " => Scan Matching Begin\n");

		// ---> compute particle weight (log likelihood)
		TIMER_BEGIN();
		max = 0;
		max_lh = -DBL_MAX;
		for( int i = 0; i < nparticle.value; i++ ) {
			double cosv, sinv;
			particles[i].x = pos_ini.value[0] + gnd::random_gaussian( error_dist.value[0] );
			particles[i].y = pos_ini.value[1] + gnd::random_gaussian( error_dist.value[1] );
			particles[i].th = pos_ini.value[2] + gnd::random_gaussian( error_dist.value[2] );
			particles[i].lh = 0;

			cosv = cos(particles[i].th);
			sinv = sin(particles[i].th);
			for( int j = 0; j < nscandata; j++ ) {
				uint32_t *v;
				double x = scandata[j].x * cosv - scandata[j].y * sinv + particles[i].x;
				double y = scandata[j].x * sinv + scandata[j].y * cosv + particles[i].y;

				if( (v = bmp.ppointer(x, y)) ) {
					// within map
					particles[i].lh += ::log( (double)*v + DBL_EPSILON );	// log likelihood (use DBL_EPSILON for avoiding 0)
				}
				else {
					// out of map
					particles[i].lh += ::log( DBL_EPSILON );	// log likelihood ( use DBL_EPSILON for avoiding 0)
				}
			}
			if( particles[i].lh > particles[max].lh){
				max = i;
				max_lh = particles[i].lh;
			}
		} // <--- compute particle weight (log likelihood)

		// ---> transform particle weight from log likelihood to likelihood ratio
		for( int i = 0; i < nparticle.value; i++ ) {
			particles[i].lh = ::exp( particles[i].lh - max_lh );
		} // <--- transform particle weight from log likelihood to likelihood ratio
		TIMER_REC(&t, &l);


		TIMER_SHOW(&t, &l);
		::fprintf(stdout, "-------- maximum likelihood particle --------\n");
		::fprintf(stdout, "         pos = %lf, %lf, %lf\n", particles[max].x, particles[max].y, particles[max].th * 180 / M_PI);
		::fprintf(stdout, "\n");

		::fprintf(stdout, " ... Finish\n");

	} // <--- operation


	{ // ---> finalize

		{ // ---> output matched scan data (global coordinate)
			double cosv = cos(particles[max].th);
			double sinv = sin(particles[max].th);
			FILE					*fp_m;					// scan data on matching position (file out)

			::fprintf(stdout, " => file-out matched scan points on global coordinate\n");


			// open output file (scan on matching position)
			fp_m = ::fopen("scan-on-match-pos.txt", "w");
			if( !fp_m ) {
				::fprintf(stdout, " error: fail to open \"%s\"\n", "scan-on-match-pos.txt");
				return -1;
			}

			// ---> fileout
			for( int i = 0; i < nscandata; i++ ) {
				::fprintf(fp_m, "%lf %lf\n",
						scandata[i].x * cosv - scandata[i].y * sinv + particles[max].x,
						scandata[i].x * sinv + scandata[i].y * cosv + particles[max].y );
			} // ---> fileout

			::fprintf(stdout, " ... output \"matched-scan-data.txt\"\n");

			// file close
			::fclose(fp_m);
		} // <--- output matched scan data (global coordinate)

		// file close

		delete[] scandata;
		delete[] particles;

		gnd::lssmap::destroy_counting_map(&cmap);
		gnd::lssmap::destroy_map(&lssmap);
	} // <--- finalize


	return 0;
}
