
#ifndef GND_ROUTE_HPP_
#define GND_ROUTE_HPP_

#include "gnd-multi-io.h"
#include "gnd-multi-math.h"

#include "gnd-queue.hpp"
#include "gnd-lib-error.h"

#include <string.h>
#include <stdint.h>

// ---> type definition
namespace gnd {
	namespace route {
		typedef int32_t seq_id_t;
		typedef struct waypoint_t {
			double x;
			double y;
		} waypoint_t;
		typedef struct waypoint_directional_t : virtual public waypoint_t {
			double theta;
		} waypoint_directional_t;

		typedef struct waypoint_named_t : virtual public waypoint_t {
			seq_id_t seq_id;
			char name[128];
		} waypoint_named_t;

		typedef struct waypoint_directional_named_t :
		public waypoint_named_t, public waypoint_directional_t {
		} waypoint_directional_named_t;




		template <class PathProp>
		struct _path_t {
			waypoint_directional_named_t end;
			double curvature;
			PathProp prop;
		};

		template <class PathProp>
		struct _path_list_t {
			waypoint_named_t start;
			queue< _path_t<PathProp> > paths;

			_path_list_t<PathProp>& operator=(const _path_list_t<PathProp>& src) {
				memcpy(&start, &src.start, sizeof(start));
				paths.clear();
				if( src.paths.size() > 0 ) {
					paths.push_back( src.paths.const_begin(), src.paths.size() );
				}
				return *this;
			}
		};


		template <class PathProp>
		class route {
		public:
			typedef seq_id_t				seq_id_t;
			typedef _path_t<PathProp> 		path_t;
			typedef path_t* 				path_pt;
			typedef const path_t* 			const_path_pt;
			typedef _path_list_t<PathProp>	path_list_t;
		private:
			seq_id_t seq_id_;
			queue< path_list_t > path_list;

		public:
			route() : seq_id_(0) {}

		public:
			int n_waypoints();

		public:
			int add_waypoint(char *name, double x, double y);
			int add_waypoint(double x, double y);

			int set_linepath(char *from, char *to, PathProp *prop);
			int set_linepath(seq_id_t from, seq_id_t to, PathProp *prop);

			int set_linepath_bidirectional(char *a, char *b, PathProp *prop);
			int set_linepath_bidirectional(seq_id_t a, seq_id_t b, PathProp *prop);

		private:
			int get_pathlist_index( char *name_start ) const;
			int get_pathlist_index( seq_id_t seq_id ) const;

		public:
			int copy_path(_path_t<PathProp> *dest, char *from, char *to);
			int copy_path(_path_t<PathProp> *dest, seq_id_t from, seq_id_t to);

			int copy_pathlist(_path_list_t<PathProp> *dest, char *waypoint);
			int copy_pathlist(_path_list_t<PathProp> *dest, seq_id_t waypoint);

			const _path_t<PathProp>* const_path_pointer(char *from, char *to) const;
			const _path_t<PathProp>* const_path_pointer(seq_id_t from, seq_id_t to) const;

			const _path_list_t<PathProp>* cosnt_pathlist_pointer( char *waypoint ) const;
			const _path_list_t<PathProp>* cosnt_pathlist_pointer( seq_id_t waypoint ) const;

		};


		typedef double _pathprop_speed_t;
		typedef struct {
			double speed;
			double left_width;
			double right_width;
			double front_extend;
			double back_extend;
		} _pathprop_area_and_speed_t;

		typedef route< _pathprop_speed_t > 					route_speed_limited;
		typedef route_speed_limited::path_t 				path_speed_limited_t;
		typedef route_speed_limited::path_list_t 			path_list_speed_limited_t;
		typedef route< _pathprop_area_and_speed_t >			route_area_speed_limited;
		typedef route_area_speed_limited::path_t	 		path_area_and_speed_limited_t;
		typedef route_area_speed_limited::path_list_t	 	path_list_area_and_speed_limited_t;
	}
} // <--- type definition






// ---> function definition
namespace gnd {

	template< >
	inline
	int queue< route::path_list_speed_limited_t >::__move__( route::path_list_speed_limited_t* dest, const route::path_list_speed_limited_t* src, uint32_t len)
	{
		uint32_t i = 0;

		for( int i = 0; i < (signed)len; i++  ) {
			dest[i] = src[i];
		}

		return 0;
	}

	template< >
	inline
	int queue< route::path_list_speed_limited_t >::__copy__( route::path_list_speed_limited_t* dest, const route::path_list_speed_limited_t* src, uint32_t len)
	{
		return __move__(dest, src, len);
	}

	template< >
	inline
	int queue< route::path_list_area_and_speed_limited_t >::__move__( route::path_list_area_and_speed_limited_t* dest, const route::path_list_area_and_speed_limited_t* src, uint32_t len)
	{
		uint32_t i = 0;

		for( int i = 0; i < (signed) len; i++  ) {
			dest[i] = src[i];
		}

		return 0;
	}

	template< >
	inline
	int queue< route::path_list_area_and_speed_limited_t >::__copy__( route::path_list_area_and_speed_limited_t* dest, const route::path_list_area_and_speed_limited_t* src, uint32_t len)
	{
		return __move__(dest, src, len);
	}



	namespace route {

		// ---> route class member function
		template <class PathProp>
		inline
		int route<PathProp>::n_waypoints() {
			return path_list.size();
		}

		template <class PathProp>
		inline
		int route<PathProp>::add_waypoint(char *name, double x, double y) {
			waypoint_named_t w;
			_path_list_t<PathProp> pl;

			// already existed
			if( (get_pathlist_index(name) ) >= 0 )  return -1;

			{ // ---> set
				w.seq_id = seq_id_++;
				strcpy(w.name, name);
				w.x = x;
				w.y = y;

				memcpy(&pl.start, &w, sizeof(w));
				pl.paths.clear();
				path_list.push_back(&pl);
			} // <--- set

			return w.seq_id;
		}

		template <class PathProp>
		inline
		int route<PathProp>::add_waypoint( double x, double y ) {
			waypoint_named_t w;
			_path_list_t<PathProp> pl;

			{ // ---> set
				w.seq_id = seq_id_++;
				memset(w.namae, 0, sizeof(w.name));
				w.x = x;
				w.y = y;

				memcpy(&pl.start, &w, sizeof(w));
				pl.paths.clear();
				path_list.push_back(&pl);
			} // <--- set

			return w.seq_id;
		}


		template <class PathProp>
		inline
		int route<PathProp>::set_linepath(char *from, char *to, PathProp *prop) {
			int f, t;
			_path_t<PathProp> p;

			if( (f = get_pathlist_index(from) ) < 0 )  return -1;
			if( (t = get_pathlist_index(to) ) < 0 )  return -1;

			p.end.x = path_list[t].start.x;
			p.end.y = path_list[t].start.y;
			p.end.theta = atan2( p.end.y - path_list[f].start.y, p.end.x - path_list[f].start.x );
			strcpy( p.end.name, path_list[t].start.name);
			p.end.seq_id = path_list[t].start.seq_id;
			p.curvature = 0;
			memcpy(&p.prop, prop, sizeof(PathProp));

			path_list[f].paths.push_back(&p);

			return 0;
		}

		template <class PathProp>
		inline
		int route<PathProp>::set_linepath(seq_id_t from, seq_id_t to, PathProp *prop) {
			int f, t;
			_path_t<PathProp> p;

			if( (f = get_pathlist_index(from) ) < 0 )  return -1;
			if( (t = get_pathlist_index(to) ) < 0 )  return -1;

			p.end.x = path_list[t].start.x;
			p.end.y = path_list[t].start.y;
			p.end.theta = atan2( p.end.y - path_list[f].start.y, p.end.x - path_list[f].start.x );
			strcpy( p.end.name, path_list[t].start.name);
			p.end.seq_id = path_list[t].start.seq_id;
			p.curvature = 0;
			memcpy(&p.prop, prop, sizeof(PathProp));

			path_list[f].paths.push_back(&p);

			return 0;
		}



		template <class PathProp>
		inline
		int route<PathProp>::set_linepath_bidirectional(char *a, char *b, PathProp *prop) {
			if( set_linepath(a,b,prop) < 0 )	return -1;
			if( set_linepath(b,a,prop) < 0 )	return -1;
			return 0;
		}

		template <class PathProp>
		inline
		int route<PathProp>::set_linepath_bidirectional(seq_id_t a, seq_id_t b, PathProp *prop) {
			if( set_linepath(a,b,prop) < 0 )	return -1;
			if( set_linepath(b,a,prop) < 0 )	return -1;
			return 0;
		}





		template <class PathProp>
		inline
		int route<PathProp>::get_pathlist_index( char *name_start ) const {
			int i;

			for( i = 0; i < (signed)path_list.size(); i++ ) {
				if( strcmp( (path_list.const_begin() + i)->start.name, name_start ) == 0  ) break;
			}
			if( i >= (signed) path_list.size() ) {
				return -1;
			}

			return i;
		}

		template <class PathProp>
		inline
		int route<PathProp>::get_pathlist_index( int seq_id ) const {
			int i;

			for( i = 0; i < (signed)path_list.size(); i++ ) {
				if( (path_list.const_begin() + i)->start.seq_id == seq_id  ) break;
			}
			if( i >= (signed) path_list.size() ) {
				return -1;
			}

			return i;
		}




		template <class PathProp>
		inline
		int route<PathProp>::copy_path(_path_t<PathProp> *dest, char *from, char *to) {
			const _path_t<PathProp> *p;
			if( !(p = const_path_pointer(from, to)) ) {
				return -1;
			}

			memcpy(dest, p, sizeof(_path_t<PathProp>) );
			return 0;
		}
		template <class PathProp>
		inline
		int route<PathProp>::copy_path(_path_t<PathProp> *dest, seq_id_t from, seq_id_t to){
			const _path_t<PathProp> *p;
			if( !(p = const_path_pointer(from, to)) ) {
				return -1;
			}

			memcpy(dest, p, sizeof(_path_t<PathProp>) );
			return 0;
		}


		template <class PathProp>
		inline
		int route<PathProp>::copy_pathlist(_path_list_t<PathProp> *dest, char *waypoint) {
			const _path_list_t<PathProp> *p;
			if( !(p = cosnt_pathlist_pointer(waypoint)) ) {
				return -1;
			}

			*dest = *p;
			return 0;
		}

		template <class PathProp>
		inline
		int route<PathProp>::copy_pathlist(_path_list_t<PathProp> *dest, seq_id_t waypoint) {
			const _path_list_t<PathProp> *p;
			if( !(p = cosnt_pathlist_pointer(waypoint)) ) {
				return -1;
			}

			*dest = *p;
			return 0;
		}


		template <class PathProp>
		inline
		const _path_t<PathProp>* route<PathProp>::const_path_pointer(char *from, char *to) const {
			int i;
			const _path_list_t<PathProp> *p;
			if( !(p = cosnt_pathlist_pointer(from)) ) return 0;

			for( i = 0; i < (signed)p->paths.size(); i++ ) {
				if( strcmp( p->paths[i].end.name, to) == 0 ) break;
			}
			if( i >= (signed) p->paths.size() ) {
				return 0;
			}

			return p->paths.const_begin() + i;
		}
		template <class PathProp>
		inline
		const _path_t<PathProp>* route<PathProp>::const_path_pointer(seq_id_t from, seq_id_t to) const {
			int i;
			const _path_list_t<PathProp> *p;
			if( !(p = cosnt_pathlist_pointer(from)) ) return 0;

			for( i = 0; i < (signed)p->paths.size(); i++ ) {
				if( (p->paths.const_begin() + i)->end.seq_id == to ) break;
			}
			if( i >= (signed) p->paths.size() ) {
				return 0;
			}

			return p->paths.const_begin() + i;
		}

		template <class PathProp>
		inline
		const _path_list_t<PathProp>* route<PathProp>::cosnt_pathlist_pointer( char *waypoint ) const {
			int i;
			if ( (i = get_pathlist_index(waypoint)) < 0) return 0;
			return path_list.const_begin() + i;
		}
		template <class PathProp>
		inline
		const _path_list_t<PathProp>* route<PathProp>::cosnt_pathlist_pointer( seq_id_t waypoint ) const {
			int i;
			if ( (i = get_pathlist_index(waypoint)) < 0) return 0;
			return path_list.const_begin() + i;
		}
		// <--- route class member function

	}
} // <--- function definition


#endif

