
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



		struct path_unit {
			waypoint_directional_named_t end;
			double curvature;

		};
		struct path {
			waypoint_named_t start;
			queue< path_unit > pth;

			path& operator=(const path& src) {
				memcpy(&start, &src.start, sizeof(start));
				pth.clear();
				if( src.pth.size() > 0 ) {
					pth.push_back( src.pth.const_begin(), src.pth.size() );
				}
				return *this;
			}
		};
		class path_net {
		public:
			typedef seq_id_t			seq_id_t;
			typedef path_unit			path_unit_t;
			typedef path_unit_t* 		path_unit_pt;
			typedef const path_unit_t* 	const_path_unit_pt;
			typedef path				path_t;
		private:
			seq_id_t seq_id_;
			queue< path_t > net_;

		public:
			path_net() : seq_id_(0) {}

		public:
			int n_waypoints();

		public:
			int add_waypoint(char *name, double x, double y);
			int add_waypoint(double x, double y);

			int set_linepath(char *from, char *to);
			int set_linepath(seq_id_t from, seq_id_t to);

			int set_linepath_bidirectional(char *a, char *b);
			int set_linepath_bidirectional(seq_id_t a, seq_id_t b);

		private:
			int get_pathlist_index( char *name_start ) const;
			int get_pathlist_index( seq_id_t seq_id ) const;

		public:
			int copy_path(path_unit *dest, char *from, char *to);
			int copy_path(path_unit *dest, seq_id_t from, seq_id_t to);

			int copy_pathlist(path *dest, char *waypoint);
			int copy_pathlist(path *dest, seq_id_t waypoint);

			const path_unit* const_path_pointer(char *from, char *to) const;
			const path_unit* const_path_pointer(seq_id_t from, seq_id_t to) const;

			const path* cosnt_pathlist_pointer( char *waypoint ) const;
			const path* cosnt_pathlist_pointer( seq_id_t waypoint ) const;

		};




		template <class PathProp>
		struct path_unit_with_property : public path_unit {
			PathProp prop;
		};
		template <class PathProp>
		struct path_with_property {
			waypoint_named_t start;
			queue< path_unit_with_property<PathProp> > paths;

			path_with_property<PathProp>& operator=(const path_with_property<PathProp>& src) {
				memcpy(&start, &src.start, sizeof(start));
				paths.clear();
				if( src.paths.size() > 0 ) {
					paths.push_back( src.paths.const_begin(), src.paths.size() );
				}
				return *this;
			}
		};
		template <class PathProp>
		class path_net_with_property {
		public:
			typedef seq_id_t							seq_id_t;
			typedef path_unit_with_property<PathProp> 	path_unit_t;
			typedef path_unit_t* 						path_unit_pt;
			typedef const path_unit_t* 					const_path_unit_pt;
			typedef path_with_property<PathProp>		path_t;
		private:
			seq_id_t seq_id_;
			queue< path_t > net_;

		public:
			path_net_with_property() : seq_id_(0) {}

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
			int copy_path(path_unit_with_property<PathProp> *dest, char *from, char *to);
			int copy_path(path_unit_with_property<PathProp> *dest, seq_id_t from, seq_id_t to);

			int copy_pathlist(path_with_property<PathProp> *dest, char *waypoint);
			int copy_pathlist(path_with_property<PathProp> *dest, seq_id_t waypoint);

			const path_unit_with_property<PathProp>* const_path_pointer(char *from, char *to) const;
			const path_unit_with_property<PathProp>* const_path_pointer(seq_id_t from, seq_id_t to) const;

			const path_with_property<PathProp>* cosnt_pathlist_pointer( char *waypoint ) const;
			const path_with_property<PathProp>* cosnt_pathlist_pointer( seq_id_t waypoint ) const;

		};


		typedef struct {
			double limit_translate;
			double limit_rotate;
		} pathprop_speed_t;
		typedef struct {
			double speed;
			double left_width;
			double right_width;
			double front_extend;
			double back_extend;
		} pathprop_area_and_speed_t;

		typedef path_net_with_property< pathprop_speed_t > 				path_net_speed_limited;
		typedef path_net_speed_limited::path_unit_t 					path_unit_speed_limited;
		typedef path_net_speed_limited::path_t 							path_speed_limited;
		typedef path_net_with_property< pathprop_area_and_speed_t >		path_net_area_speed_limited;
		typedef path_net_area_speed_limited::path_unit_t	 			path_unit_area_and_speed_limited;
		typedef path_net_area_speed_limited::path_t						path_area_and_speed_limited;
	}
} // <--- type definition






// ---> function definition
namespace gnd {

	template< >
	inline
	int queue< route::path >::__move__( route::path* dest, const route::path* src, uint32_t len)
	{
		uint32_t i = 0;

		for( int i = 0; i < (signed)len; i++  ) {
			dest[i] = src[i];
		}

		return 0;
	}

	template< >
	inline
	int queue< route::path >::__copy__( route::path* dest, const route::path* src, uint32_t len)
	{
		return __move__(dest, src, len);
	}



	template< >
	inline
	int queue< route::path_speed_limited >::__move__( route::path_speed_limited* dest, const route::path_speed_limited* src, uint32_t len)
	{
		uint32_t i = 0;

		for( int i = 0; i < (signed)len; i++  ) {
			dest[i] = src[i];
		}

		return 0;
	}

	template< >
	inline
	int queue< route::path_speed_limited >::__copy__( route::path_speed_limited* dest, const route::path_speed_limited* src, uint32_t len)
	{
		return __move__(dest, src, len);
	}


	template< >
	inline
	int queue< route::path_area_and_speed_limited >::__move__( route::path_area_and_speed_limited* dest, const route::path_area_and_speed_limited* src, uint32_t len)
	{
		uint32_t i = 0;

		for( int i = 0; i < (signed) len; i++  ) {
			dest[i] = src[i];
		}

		return 0;
	}

	template< >
	inline
	int queue< route::path_area_and_speed_limited >::__copy__( route::path_area_and_speed_limited* dest, const route::path_area_and_speed_limited* src, uint32_t len)
	{
		return __move__(dest, src, len);
	}





	// ---> function definition (path_net)
	namespace route {

		inline
		int path_net::n_waypoints() {
			return net_.size();
		}

		inline
		int path_net::add_waypoint(char *name, double x, double y) {
			waypoint_named_t w;
			path pl;

			// already existed
			if( (get_pathlist_index(name) ) >= 0 )  return -1;

			{ // ---> set
				w.seq_id = seq_id_++;
				strcpy(w.name, name);
				w.x = x;
				w.y = y;

				memcpy(&pl.start, &w, sizeof(w));
				pl.pth.clear();
				net_.push_back(&pl);
			} // <--- set

			return w.seq_id;
		}

		inline
		int path_net::add_waypoint( double x, double y ) {
			waypoint_named_t w;
			path pl;

			{ // ---> set
				w.seq_id = seq_id_++;
				memset(w.name, 0, sizeof(w.name));
				w.x = x;
				w.y = y;

				memcpy(&pl.start, &w, sizeof(w));
				pl.pth.clear();
				net_.push_back(&pl);
			} // <--- set

			return w.seq_id;
		}


		inline
		int path_net::set_linepath(char *from, char *to) {
			int f, t;
			path_unit p;

			if( (f = get_pathlist_index(from) ) < 0 )  return -1;
			if( (t = get_pathlist_index(to) ) < 0 )  return -1;

			p.end.x = net_[t].start.x;
			p.end.y = net_[t].start.y;
			p.end.theta = atan2( p.end.y - net_[f].start.y, p.end.x - net_[f].start.x );
			strcpy( p.end.name, net_[t].start.name);
			p.end.seq_id = net_[t].start.seq_id;
			p.curvature = 0;

			net_[f].pth.push_back(&p);

			return 0;
		}

		inline
		int path_net::set_linepath(seq_id_t from, seq_id_t to) {
			int f, t;
			path_unit p;

			if( (f = get_pathlist_index(from) ) < 0 )  return -1;
			if( (t = get_pathlist_index(to) ) < 0 )  return -1;

			p.end.x = net_[t].start.x;
			p.end.y = net_[t].start.y;
			p.end.theta = atan2( p.end.y - net_[f].start.y, p.end.x - net_[f].start.x );
			strcpy( p.end.name, net_[t].start.name);
			p.end.seq_id = net_[t].start.seq_id;
			p.curvature = 0;

			net_[f].pth.push_back(&p);

			return 0;
		}



		inline
		int path_net::set_linepath_bidirectional(char *a, char *b) {
			if( set_linepath(a,b) < 0 )	return -1;
			if( set_linepath(b,a) < 0 )	return -1;
			return 0;
		}

		inline
		int path_net::set_linepath_bidirectional(seq_id_t a, seq_id_t b) {
			if( set_linepath(a,b) < 0 )	return -1;
			if( set_linepath(b,a) < 0 )	return -1;
			return 0;
		}



		inline
		int path_net::get_pathlist_index( char *name_start ) const {
			int i;

			for( i = 0; i < (signed)net_.size(); i++ ) {
				if( strcmp( (net_.const_begin() + i)->start.name, name_start ) == 0  ) break;
			}
			if( i >= (signed) net_.size() ) {
				return -1;
			}

			return i;
		}

		inline
		int path_net::get_pathlist_index( int seq_id ) const {
			int i;

			for( i = 0; i < (signed)net_.size(); i++ ) {
				if( (net_.const_begin() + i)->start.seq_id == seq_id  ) break;
			}
			if( i >= (signed) net_.size() ) {
				return -1;
			}

			return i;
		}




		inline
		int path_net::copy_path(path_unit *dest, char *from, char *to) {
			const path_unit *p;
			if( !(p = const_path_pointer(from, to)) ) {
				return -1;
			}

			memcpy(dest, p, sizeof(path_unit) );
			return 0;
		}
		inline
		int path_net::copy_path(path_unit *dest, seq_id_t from, seq_id_t to){
			const path_unit *p;
			if( !(p = const_path_pointer(from, to)) ) {
				return -1;
			}

			memcpy(dest, p, sizeof(path_unit) );
			return 0;
		}


		inline
		int path_net::copy_pathlist(path *dest, char *waypoint) {
			const path *p;
			if( !(p = cosnt_pathlist_pointer(waypoint)) ) {
				return -1;
			}

			*dest = *p;
			return 0;
		}

		inline
		int path_net::copy_pathlist(path *dest, seq_id_t waypoint) {
			const path *p;
			if( !(p = cosnt_pathlist_pointer(waypoint)) ) {
				return -1;
			}

			*dest = *p;
			return 0;
		}


		inline
		const path_unit* path_net::const_path_pointer(char *from, char *to) const {
			int i;
			const path_t *p;
			if( !(p = cosnt_pathlist_pointer(from)) ) return 0;

			for( i = 0; i < (signed)p->pth.size(); i++ ) {
				if( strcmp( (p->pth.const_begin() + i)->end.name, to) == 0 ) break;
			}
			if( i >= (signed) p->pth.size() ) {
				return 0;
			}

			return p->pth.const_begin() + i;
		}
		inline
		const path_unit* path_net::const_path_pointer(seq_id_t from, seq_id_t to) const {
			int i;
			const path *p;
			if( !(p = cosnt_pathlist_pointer(from)) ) return 0;

			for( i = 0; i < (signed)p->pth.size(); i++ ) {
				if( (p->pth.const_begin() + i)->end.seq_id == to ) break;
			}
			if( i >= (signed) p->pth.size() ) {
				return 0;
			}

			return p->pth.const_begin() + i;
		}

		inline
		const path* path_net::cosnt_pathlist_pointer( char *waypoint ) const {
			int i;
			if ( (i = get_pathlist_index(waypoint)) < 0) return 0;
			return net_.const_begin() + i;
		}
		inline
		const path* path_net::cosnt_pathlist_pointer( seq_id_t waypoint ) const {
			int i;
			if ( (i = get_pathlist_index(waypoint)) < 0) return 0;
			return net_.const_begin() + i;
		}
		// <--- route class member function

	} // <--- function definition (path_net)



	// ---> function definition (path_net_with_propert)
	namespace route {

		// ---> route class member function
		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::n_waypoints() {
			return net_.size();
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::add_waypoint(char *name, double x, double y) {
			waypoint_named_t w;
			path_with_property<PathProp> pl;

			// already existed
			if( (get_pathlist_index(name) ) >= 0 )  return -1;

			{ // ---> set
				w.seq_id = seq_id_++;
				strcpy(w.name, name);
				w.x = x;
				w.y = y;

				memcpy(&pl.start, &w, sizeof(w));
				pl.paths.clear();
				net_.push_back(&pl);
			} // <--- set

			return w.seq_id;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::add_waypoint( double x, double y ) {
			waypoint_named_t w;
			path_with_property<PathProp> pl;

			{ // ---> set
				w.seq_id = seq_id_++;
				memset(w.name, 0, sizeof(w.name));
				w.x = x;
				w.y = y;

				memcpy(&pl.start, &w, sizeof(w));
				pl.paths.clear();
				net_.push_back(&pl);
			} // <--- set

			return w.seq_id;
		}


		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::set_linepath(char *from, char *to, PathProp *prop) {
			int f, t;
			path_unit_with_property<PathProp> p;

			if( (f = get_pathlist_index(from) ) < 0 )  return -1;
			if( (t = get_pathlist_index(to) ) < 0 )  return -1;

			p.end.x = net_[t].start.x;
			p.end.y = net_[t].start.y;
			p.end.theta = atan2( p.end.y - net_[f].start.y, p.end.x - net_[f].start.x );
			strcpy( p.end.name, net_[t].start.name);
			p.end.seq_id = net_[t].start.seq_id;
			p.curvature = 0;
			memcpy(&p.prop, prop, sizeof(PathProp));

			net_[f].paths.push_back(&p);

			return 0;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::set_linepath(seq_id_t from, seq_id_t to, PathProp *prop) {
			int f, t;
			path_unit_with_property<PathProp> p;

			if( (f = get_pathlist_index(from) ) < 0 )  return -1;
			if( (t = get_pathlist_index(to) ) < 0 )  return -1;

			p.end.x = net_[t].start.x;
			p.end.y = net_[t].start.y;
			p.end.theta = atan2( p.end.y - net_[f].start.y, p.end.x - net_[f].start.x );
			strcpy( p.end.name, net_[t].start.name);
			p.end.seq_id = net_[t].start.seq_id;
			p.curvature = 0;
			memcpy(&p.prop, prop, sizeof(PathProp));

			net_[f].paths.push_back(&p);

			return 0;
		}



		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::set_linepath_bidirectional(char *a, char *b, PathProp *prop) {
			if( set_linepath(a,b,prop) < 0 )	return -1;
			if( set_linepath(b,a,prop) < 0 )	return -1;
			return 0;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::set_linepath_bidirectional(seq_id_t a, seq_id_t b, PathProp *prop) {
			if( set_linepath(a,b,prop) < 0 )	return -1;
			if( set_linepath(b,a,prop) < 0 )	return -1;
			return 0;
		}





		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::get_pathlist_index( char *name_start ) const {
			int i;

			for( i = 0; i < (signed)net_.size(); i++ ) {
				if( strcmp( (net_.const_begin() + i)->start.name, name_start ) == 0  ) break;
			}
			if( i >= (signed) net_.size() ) {
				return -1;
			}

			return i;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::get_pathlist_index( int seq_id ) const {
			int i;

			for( i = 0; i < (signed)net_.size(); i++ ) {
				if( (net_.const_begin() + i)->start.seq_id == seq_id  ) break;
			}
			if( i >= (signed) net_.size() ) {
				return -1;
			}

			return i;
		}




		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::copy_path(path_unit_with_property<PathProp> *dest, char *from, char *to) {
			const path_unit_with_property<PathProp> *p;
			if( !(p = const_path_pointer(from, to)) ) {
				return -1;
			}

			memcpy(dest, p, sizeof(path_unit_with_property<PathProp>) );
			return 0;
		}
		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::copy_path(path_unit_with_property<PathProp> *dest, seq_id_t from, seq_id_t to){
			const path_unit_with_property<PathProp> *p;
			if( !(p = const_path_pointer(from, to)) ) {
				return -1;
			}

			memcpy(dest, p, sizeof(path_unit_with_property<PathProp>) );
			return 0;
		}


		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::copy_pathlist(path_with_property<PathProp> *dest, char *waypoint) {
			const path_with_property<PathProp> *p;
			if( !(p = cosnt_pathlist_pointer(waypoint)) ) {
				return -1;
			}

			*dest = *p;
			return 0;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::copy_pathlist(path_with_property<PathProp> *dest, seq_id_t waypoint) {
			const path_with_property<PathProp> *p;
			if( !(p = cosnt_pathlist_pointer(waypoint)) ) {
				return -1;
			}

			*dest = *p;
			return 0;
		}


		template <class PathProp>
		inline
		const path_unit_with_property<PathProp>* path_net_with_property<PathProp>::const_path_pointer(char *from, char *to) const {
			int i;
			const path_with_property<PathProp> *p;
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
		const path_unit_with_property<PathProp>* path_net_with_property<PathProp>::const_path_pointer(seq_id_t from, seq_id_t to) const {
			int i;
			const path_with_property<PathProp> *p;
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
		const path_with_property<PathProp>* path_net_with_property<PathProp>::cosnt_pathlist_pointer( char *waypoint ) const {
			int i;
			if ( (i = get_pathlist_index(waypoint)) < 0) return 0;
			return net_.const_begin() + i;
		}
		template <class PathProp>
		inline
		const path_with_property<PathProp>* path_net_with_property<PathProp>::cosnt_pathlist_pointer( seq_id_t waypoint ) const {
			int i;
			if ( (i = get_pathlist_index(waypoint)) < 0) return 0;
			return net_.const_begin() + i;
		}
		// <--- route class member function

	} // <--- function definition (path_net_with_propert)

} // <--- function definition


#endif

