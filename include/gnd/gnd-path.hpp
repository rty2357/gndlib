
#ifndef GND_ROUTE_HPP_
#define GND_ROUTE_HPP_

#include "gnd-multi-platform.h"
#include "gnd-multi-io.h"
#include "gnd-multi-math.h"

#include "gnd-util.h"
#include "gnd-lib-error.h"
#include "gnd-queue.hpp"

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <float.h>

// ---> type definition
namespace gnd {
	namespace path {
		typedef int32_t seq_id_t;
		typedef struct waypoint_t {
			double x;
			double y;
		} waypoint_t;
		typedef struct waypoint_directional_t : virtual public waypoint_t {
			double theta;
		} waypoint_directional_t;

		typedef struct waypoint_named_t : virtual public waypoint_t {
			char name[128];
		} waypoint_named_t;

		typedef struct waypoint_directional_named_t :
		public waypoint_named_t, public waypoint_directional_t {
		} waypoint_directional_named_t;



		struct path_unit {
			waypoint_directional_named_t end;
			double curvature;

		};

		template <class PathProp>
		struct path_unit_with_property : public path_unit {
			PathProp prop;
		};
		template <class PathProp>
		struct waypoint_and_paths_with_property {
			waypoint_named_t start;
			queue< path_unit_with_property<PathProp> > path;

			waypoint_and_paths_with_property<PathProp>& operator=(const waypoint_and_paths_with_property<PathProp>& src) {
				memcpy(&start, &src.start, sizeof(start));
				path.clear();
				if( src.path.size() > 0 ) {
					path.push_back( src.path.const_begin(), src.path.size() );
				}
				return *this;
			}
		};


		// path from start to destination
		template <class PathProp>
		struct path_with_property {
			waypoint_named_t start;
			queue< path_unit_with_property<PathProp> > path;
			path_with_property<PathProp>& operator=(const path_with_property<PathProp>& src) {
				memcpy(&start, &src.start, sizeof(start));
				path.clear();
				if( src.path.size() > 0 ) {
					path.push_back( src.path.const_begin(), src.path.size() );
				}
				for( int i = 0; i < path.size(); i++ ) {
					fprintf(stdout, "operator = %s\n", path[i].end.name );
				}
				return *this;
			}
		};

		template <class PathProp>
		class path_net_with_property {
		public:
			typedef PathProp									property_t;
			typedef path_unit_with_property<PathProp> 			path_unit_t;
			typedef path_unit_t* 								path_unit_pt;
			typedef const path_unit_t* 							const_path_unit_pt;
			typedef waypoint_and_paths_with_property<PathProp>	waypoint_and_paths_t;
			typedef path_with_property<PathProp>				path_t;
			typedef waypoint_named_t							waypoint_t;
		private:
			queue< waypoint_and_paths_t > net_;

		public:
			path_net_with_property(){}

		public:
			int clear();

		public:
			int n_waypoints();

		public:
			int add_waypoint(const char *name, double x, double y);
			int erase_waypoint(const char *name);
			int set_waypoint(const char *name, double x, double y);
			int get_waypoint(const char *name, double *x, double *y);
			int rename_waypoint(const char *name, const char *new_name);
		public:
			int add_linepath(const char *from, const char *to, PathProp *prop);
			int add_linepath_bidirectional(const char *a, const char *b, PathProp *prop);
			int erase_path(const char *from, const char *to);
			int erase_path_bidirectional(const char *a, const char *b);
			int set_path_property(const char *from, const char *to, PathProp *prop);
			int get_path_property(const char *from, const char *to, PathProp *prop);

		public:
			int name_waypoint(int index, char *name) const;
			int index_waypoint( const char *name ) const;

		public:
			int copy_path_unit(path_unit_with_property<PathProp> *dest, const char *from, const char *to);
			int copy_path_list(waypoint_and_paths_with_property<PathProp> *dest, const char *waypoint);
			const path_unit_with_property<PathProp>* const_pointer_path_unit(const char *from, const char *to);
			const waypoint_and_paths_with_property<PathProp>* cosnt_pointer_path_list( const char *waypoint ) const;
			path_unit_with_property<PathProp>* pointer_path_unit(const char *from, const char *to);
			waypoint_and_paths_with_property<PathProp>* pointer_path_list( const char *waypoint );

		public:
			int find_path_dijkstra( path_with_property<PathProp>* p, const char *from, const char *to);

		public:
			waypoint_and_paths_with_property<PathProp>& operator[](int i);
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






// ---> function definition (path_net_with_propert)
namespace gnd {

	template< >
	inline
	int queue< path::path_speed_limited >::__move__( path::path_speed_limited* dest, const path::path_speed_limited* src, uint32_t len)
	{
		uint32_t i = 0;

		for( int i = 0; i < (signed)len; i++  ) {
			dest[i] = src[i];
		}

		return 0;
	}

	template< >
	inline
	int queue< path::path_speed_limited >::__copy__( path::path_speed_limited* dest, const path::path_speed_limited* src, uint32_t len)
	{
		return __move__(dest, src, len);
	}




	template< >
	inline
	int queue< path::path_net_speed_limited::waypoint_and_paths_t >::__move__( path::path_net_speed_limited::waypoint_and_paths_t* dest, const path::path_net_speed_limited::waypoint_and_paths_t* src, uint32_t len)
	{
		uint32_t i = 0;

		for( int i = 0; i < (signed)len; i++  ) {
			dest[i] = src[i];
		}

		return 0;
	}

	template< >
	inline
	int queue< path::path_net_speed_limited::waypoint_and_paths_t >::__copy__( path::path_net_speed_limited::waypoint_and_paths_t* dest, const path::path_net_speed_limited::waypoint_and_paths_t* src, uint32_t len)
	{
		return __move__(dest, src, len);
	}



	template< >
	inline
	int queue< path::path_area_and_speed_limited >::__move__( path::path_area_and_speed_limited* dest, const path::path_area_and_speed_limited* src, uint32_t len)
	{
		uint32_t i = 0;

		for( int i = 0; i < (signed) len; i++  ) {
			dest[i] = src[i];
		}

		return 0;
	}

	template< >
	inline
	int queue< path::path_area_and_speed_limited >::__copy__( path::path_area_and_speed_limited* dest, const path::path_area_and_speed_limited* src, uint32_t len)
	{
		return __move__(dest, src, len);
	}


	namespace path {

		// ---> route class member function
		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::clear() {
			net_.clear();
			return 0;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::n_waypoints() {
			return net_.size();
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::add_waypoint(const char *name, double x, double y) {
			waypoint_and_paths_with_property<PathProp> wap;

			// exception: already existed
			fprintf(stdout, "call add_waypoint\n");
			if( index_waypoint(name) >= 0 )  return -1;

			{ // ---> set
				strcpy(wap.start.name, name);
				wap.start.x = x;
				wap.start.y = y;
				wap.path.clear();

				net_.push_back(&wap);
				fprintf(stdout, "  ... %d, %s\n", net_.size() - 1, net_[net_.size() - 1].start.name);
			} // <--- set

			fprintf(stdout, "ok\n");
			return n_waypoints();
		}


		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::erase_waypoint(const char *name) {
			int i, j;
			// exception: not existing
			if( (i = index_waypoint(name)) < 0 )  return -1;

			// erase waypoint
			net_.erase(i);
			// ease path to the waypoint
			for( i = 0; i < (signed)net_.size(); i++ ) {
				for( j = 0; j < (signed)net_[i].path.size(); j++) {
					if( strcmp(net_[i].path[j].end.name, name) == 0 ) {
						net_[i].path.erase(j);
					}
				}
			}

			return n_waypoints();
		}


		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::set_waypoint(const char *name, double x, double y) {
			int i, j;
			path_unit_pt p;
			// exception: not existing
			if( (i = index_waypoint(name)) < 0 )  return -1;

			net_[i].start.x = x;
			net_[i].start.y = y;

			// set path end
			for( j = 0; j < (signed)net_.size(); j++ ) {
				if( (p = pointer_path_unit( net_[j].start.name, name )) ) {
					p->end.x = x;
					p->end.y = y;
				}
			}

			return i;
		}


		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::get_waypoint(const char *name, double *x, double *y) {
			int i;
			// exception: not existing
			if( (i = index_waypoint(name)) < 0 )  return -1;

			*x = net_[i].start.x;
			*y = net_[i].start.y;
			return i;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::rename_waypoint(const char *name, const char *new_name){
			int i, j;
			// exception: not existing
			if( (i = index_waypoint(name)) < 0 )  return -1;

			strcpy(net_[i].start.name, new_name);

			// ease path to the waypoint
			for( i = 0; i < (signed)net_.size(); i++ ) {
				for( j = 0; j < (signed)net_[i].path.size(); j++) {
					if( strcmp(net_[i].path[j].end.name, name) == 0 ) {
						strcpy(net_[i].path[j].end.name, new_name);
						break;
					}
				}
			}
			return i;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::add_linepath(const char *from, const char *to, PathProp *prop) {
			int i, f, t;
			path_unit_with_property<PathProp> p;

			if( (f = index_waypoint(from) ) < 0 )  return -1;
			if( (t = index_waypoint(to) ) < 0 )  return -1;
			// already exist
			for( i = 0; i < (signed)net_[f].path.size(); i++ ) {
				if( strcmp(net_[t].start.name, net_[f].path[i].end.name) == 0 ) {
					return -1;
				}
			}

			p.end.x = net_[t].start.x;
			p.end.y = net_[t].start.y;
			p.end.theta = atan2( p.end.y - net_[f].start.y, p.end.x - net_[f].start.x );
			strcpy( p.end.name, net_[t].start.name);
			p.curvature = 0;
			memcpy(&p.prop, prop, sizeof(PathProp));

			net_[f].path.push_back(&p);

			return 0;
		}



		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::add_linepath_bidirectional(const char *a, const char *b, PathProp *prop) {
			if( add_linepath(a,b,prop) < 0 )	return -1;
			if( add_linepath(b,a,prop) < 0 )	return -1;
			return 0;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::erase_path(const char *from, const char *to) {
			int f, t;
			int i;

			if( (f = index_waypoint(from) ) < 0 || net_[f].start.name[0] == '\0' )  return -1;
			if( (t = index_waypoint(to) ) < 0 || net_[t].start.name[0] == '\0' )  return -1;

			for( i = 0; i < (signed)net_[f].path.size(); i++ ) {
				if( strcmp(net_[f].path[i].end.name, net_[t].start.name) == 0 ) {
					net_[f].path.erase(i);
				}
			}
			return 0;
		}
		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::erase_path_bidirectional(const char *a, const char *b) {
			erase_path(a, b);
			erase_path(b, a);
			return 0;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::get_path_property(const char *from, const char *to, PathProp *prop) {
			path_unit_t *pu;
			if( (pu = pointer_path_unit(from, to)) == 0 ) return -1;
			*prop = pu->prop;
			return 0;
		}

		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::set_path_property(const char *from, const char *to, PathProp *prop) {
			path_unit_t *pu;
			if( (pu = pointer_path_unit(from, to)) == 0 ) return -1;
			pu->prop = *prop;
			return 0;
		}


		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::name_waypoint(int index, char *name) const {
			if( index < 0 || index >= net_.size() ) return -1;
			if( (net_ .const_begin() + index)->start.name[0] == '\0' ) return -1;
			strcpy(name, (net_ .const_begin() + index)->start.name);
			return 0;
		}


		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::index_waypoint( const char *name_start ) const {
			int i;

			for( i = 0; i < (signed)net_.size(); i++ ) {
				if( (net_.const_begin() + i)->start.name[0] == '\0' ) continue;
				if( strcmp( (net_.const_begin() + i)->start.name, name_start ) == 0  ) break;
			}
			if( i >= (signed)net_.size() ) {
				return -1;
			}

			return i;
		}



		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::copy_path_unit(path_unit_with_property<PathProp> *dest, const char *from, const char *to) {
			const path_unit_with_property<PathProp> *p;
			if( !(p = const_pointer_path_unit(from, to)) ) {
				return -1;
			}

			memcpy(dest, p, sizeof(path_unit_with_property<PathProp>) );
			return 0;
		}


		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::copy_path_list(waypoint_and_paths_with_property<PathProp> *dest, const char *waypoint) {
			const waypoint_and_paths_with_property<PathProp> *p;
			if( !(p = cosnt_pointer_path_list(waypoint)) ) {
				return -1;
			}

			*dest = *p;
			return 0;
		}


		template <class PathProp>
		inline
		const waypoint_and_paths_with_property<PathProp>* path_net_with_property<PathProp>::cosnt_pointer_path_list( const char *waypoint ) const {
			int i;
			if ( (i = index_waypoint(waypoint)) < 0) return 0;
			return net_.const_begin() + i;
		}

		template <class PathProp>
		inline
		const path_unit_with_property<PathProp>* path_net_with_property<PathProp>::const_pointer_path_unit(const char *from, const char *to) {
			int i, j;
			if ( (i = index_waypoint(from)) < 0) return 0;

			for( j = 0; j < (signed)net_[i].path.size(); j++ ) {
				if( strcmp( net_[i].path[j].end.name, to) == 0 ) break;
			}
			if( j >= (signed) net_[i].path.size() ) {
				return 0;
			}
			return net_[i].path.const_begin() + j;
		}


		template <class PathProp>
		inline
		waypoint_and_paths_with_property<PathProp>* path_net_with_property<PathProp>::pointer_path_list( const char *waypoint ){
			int i;
			if ( (i = index_waypoint(waypoint)) < 0) return 0;
			return net_.begin() + i;
		}

		template <class PathProp>
		inline
		path_unit_with_property<PathProp>* path_net_with_property<PathProp>::pointer_path_unit(const char *from, const char *to) {
			int i, j;
			if ( (i = index_waypoint(from)) < 0) return 0;

			for( j = 0; j < (signed)net_[i].path.size(); j++ ) {
				if( strcmp( net_[i].path[j].end.name, to) == 0 ) break;
			}
			if( j >= (signed) net_[i].path.size() ) {
				return 0;
			}
			return net_[i].path.begin() + j;
		}


		// ---> search path
		template <class PathProp>
		inline
		int path_net_with_property<PathProp>::find_path_dijkstra( path_with_property<PathProp>* p, const char *start, const char *dest) {
			int i_start;
			int i_dest;
			int i;
			struct node_t {
				double cost;
				bool done;
				queue< int > path;
			};
			queue< node_t > costs;
			node_t ws;

			if( (i_start = index_waypoint(start)) < 0 )	return -1;
			if( (i_dest = index_waypoint(dest)) < 0 )	return -1;

			{ // ---> initialize cost
				costs.resize( n_waypoints() );
				for( i = 0; i < (signed)costs.size(); i++ ) {
					costs[i].cost = DBL_MAX;
					costs[i].done = false;
					costs[i].path.clear();
				}
				// at start, cost is 0
				costs[i_start].cost = 0;
				costs[i_start].path.push_back(&i_start);
			} // <--- initialize cost


			{ // ---> search

				// ---> scanning loop
				while(1) {
					int i_min = -1;
					{ // ---> get minimum cost node
						double min;

						min = DBL_MAX;
						for( i = 0; i < (signed)costs.size(); i++ ) {
							if( !costs[i].done && costs[i].cost < min ) {
								min = costs[i].cost;
								i_min = i;
							}
						}

						// exception: search all waypoint connecting from start
						if( i_min == -1) {
							return -1;
						}
						// exception:  reach to destination
						else if( i_min == i_dest ) {
							// finish
							costs[i_min].path.push_back( &i_min );
							break;
						}
					} // <--- get minimum cost node


					{ // ---> update nodes cost
						double cost_;
						for( i = 0; i < (signed)net_[i_min].path.size(); i++ ) {
							int index_ = index_waypoint(net_[i_min].path[i].end.name);
							if( costs[index_].done ) continue;
							cost_ = costs[i_min].cost +
									sqrt( gnd_square( net_[i_min].path[i].end.x - net_[i_min].start.x ) + gnd_square( net_[i_min].path[i].end.y - net_[i_min].start.y ) );

							// update
							if( cost_ < costs[index_].cost ) {
								costs[index_].cost = cost_;
								costs[index_].path.clear();
								costs[index_].path.push_back( costs[i_min].path.begin(), costs[i_min].path.size() );
								costs[index_].path.push_back( &i_min );
							}

						}
					} // <--- update nodes cost

					costs[i_min].done = true;
				} // <--- scanning loop
			} // <--- search

			{ // ---> packing
				int j = costs[i_dest].path[0];

				{ // ---> set start
					p->start = net_[j].start;
				} // <--- set start

				// set path
				p->path.clear();
				for( i = 1; i < (signed)costs[i_dest].path.size(); i++ ) {
					p->path.push_back( const_pointer_path_unit( net_[ costs[i_dest].path[i-1] ].start.name,  net_[ costs[i_dest].path[i] ].start.name ) );
				}
			} // <--- packing

			return 0;
		}
		// <--- search path

		template <class PathProp>
		inline
		waypoint_and_paths_with_property<PathProp>& path_net_with_property<PathProp>::operator[](int i) {
			return net_[i];
		}

	} // <--- function definition (path_net_with_propert)

} // <--- function definition


#endif

