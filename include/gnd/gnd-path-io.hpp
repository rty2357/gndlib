/*
 * gnd-path-file.hpp
 *
 *  Created on: 2014/12/18
 *      Author: tyamada
 */

#ifndef GND_PATH_IO_HPP_
#define GND_PATH_IO_HPP_


#include "gnd/gnd-path.hpp"
#include "gnd/gnd-config-file.hpp"
#include "gnd/gnd-util.h"

// ---> const definition
namespace gnd {
	namespace path {
		static const gnd::conf::parameter_array<char, 128> Default_Item_Waypoint = {
				"waypoints",
		};
		static const gnd::conf::parameter_array<char, 128> Default_Item_Path_List = {
				"path",
		};
		static const gnd::conf::parameter_array<char, 128> Default_Item_Path_Start = {
				"start",
		};
		static const gnd::conf::parameter_array<char, 128> Default_Item_Path_Destination = {
				"destination",
		};
		static const char Format_Item_Path[] = "from_%s_to_%s";
	}
} // <--- const definition



// <--- function declaration
namespace gnd {
	namespace path {
		static const gnd::conf::param_double Default_Item_Path_Limit_Speed_Translate = {
				"speed_translate",
		};
		static const gnd::conf::param_double Default_Item_Path_Limit_Speed_Rotate = {
				"speed_rotate"
		};

		int fread(const char* fname, path_net_speed_limited *path);
		int fwrite(const char* fname, path_net_speed_limited *path);
	}
} // <--- function declaration



// ---> function definition
namespace gnd {
	namespace path {

		inline
		int fread(const char* fname, path_net_speed_limited *net) {
			gnd::conf::file_stream root;
			gnd::conf::configuration *child;

			if( root.read(fname) < 0 ) {
				return -1;
			}
			net->clear();

			// ---> waypoint
			if( (child = root.child_find(Default_Item_Waypoint.item)) == 0 ) {
				// nothing to do
			}
			else {
				int i = 0;

				for( i = 0; i < (signed) child->nchild(); i++ ) {
					net->add_waypoint(
							(*child)[i].name(),
							atof( (*child)[i][0].value() ),
							atof( (*child)[i][1].value() ) );
				}
			} // <--- waypoint

			// ---> path
			if( (child = root.child_find(Default_Item_Path_List.item)) == 0 ) {
				// nothing to do
			}
			else {
				int i = 0;
				path_net_speed_limited::property_t prop;
				gnd::conf::parameter_array<char, 128> start;
				gnd::conf::parameter_array<char, 128> dest;
				gnd::conf::param_double speed_translate;
				gnd::conf::param_double speed_rotate;

				for( i = 0; i < (signed) child->nchild(); i++ ) {
					memcpy(&start,				&Default_Item_Path_Start,					sizeof(Default_Item_Path_Start));
					memcpy(&dest,				&Default_Item_Path_Destination,				sizeof(Default_Item_Path_Destination));
					memcpy(&speed_translate,	&Default_Item_Path_Limit_Speed_Translate,	sizeof(Default_Item_Path_Limit_Speed_Translate));
					memcpy(&speed_rotate,		&Default_Item_Path_Limit_Speed_Rotate,		sizeof(Default_Item_Path_Limit_Speed_Rotate));

					gnd::conf::get_parameter(&(*child)[i], &start);
					gnd::conf::get_parameter(&(*child)[i], &dest);
					gnd::conf::get_parameter(&(*child)[i], &speed_translate);
					gnd::conf::get_parameter(&(*child)[i], &speed_rotate);

					prop.limit_translate = speed_translate.value;
					prop.limit_rotate = gnd_deg2ang(speed_rotate.value);

					net->add_linepath(start.value, dest.value, &prop );
				}
			} // <--- path

			return 0;
		}


		inline
		int fwrite(const char* fname, path_net_speed_limited *net) {
			typedef gnd::conf::parameter_array<double, 2> param_waypoint_t;
			gnd::conf::file_stream root;

			{ // ---> set waypoint
				int i;
				gnd::conf::configuration *waypoints;
				param_waypoint_t waypoint;

				root.child_push_back(Default_Item_Waypoint.item, 0, 0);
				waypoints = root.child_find(Default_Item_Waypoint.item, 0);
				memset( waypoint.comment, 0, sizeof(waypoint.comment) );
				for( i = 0; i < net->n_waypoints(); i++) {
					strcpy( waypoint.item, (*net)[i].start.name );
					waypoint.value[0] = (*net)[i].start.x;
					waypoint.value[1] = (*net)[i].start.y;

					gnd::conf::set_parameter( waypoints, &waypoint );
				}

			} // <--- set waypoints


			{ // ---> set path list
				int i = 0;
				int j;
				gnd::conf::configuration *path_list;
				gnd::conf::configuration *path;
				gnd::conf::parameter_array<char, 128> start;
				gnd::conf::parameter_array<char, 128> dest;
				gnd::conf::param_double speed_translate;
				gnd::conf::param_double speed_rotate;
				char path_name[128];

				// make path list item
				root.child_push_back(Default_Item_Path_List.item, 0, 0);
				{ // ---> push paths
					path_list = root.child_find(Default_Item_Path_List.item, 0);
					// init parameter items
					memcpy(&start,				&Default_Item_Path_Start,					sizeof(Default_Item_Path_Start));
					memcpy(&dest,				&Default_Item_Path_Destination,				sizeof(Default_Item_Path_Destination));
					memcpy(&speed_translate,	&Default_Item_Path_Limit_Speed_Translate,	sizeof(Default_Item_Path_Limit_Speed_Translate));
					memcpy(&speed_rotate,		&Default_Item_Path_Limit_Speed_Rotate,		sizeof(Default_Item_Path_Limit_Speed_Rotate));

					for(j = 0; j < (signed)net->n_waypoints(); j++) {
						// get start waypoint name
						strcpy(start.value, (*net)[j].start.name);

						for(i = 0; i < (signed)(*net)[j].path.size(); i++){
							// get destination waypoint name
							strcpy(dest.value, (*net)[j].path[i].end.name);

							// path name
							sprintf(path_name, Format_Item_Path, start.value, dest.value);

							// speed set
							speed_translate.value = (*net)[j].path[i].prop.limit_translate;
							speed_rotate.value = gnd_ang2deg( (*net)[j].path[i].prop.limit_rotate );

							// set
							path_list->child_push_back(path_name, 0, 0);
							path = path_list->child_find(path_name, 0 );
							gnd::conf::set_parameter( path, &start );
							gnd::conf::set_parameter( path, &dest );
							gnd::conf::set_parameter( path, &speed_translate );
							gnd::conf::set_parameter( path, &speed_rotate );
						}
					}
				} // <--- push paths
			} // <--- set path list

			if( root.write(fname) < 0) {
				return -1;
			}
			return 0;
		}
	}
} // <--- function definition



#endif /* GND_PATH_IO_HPP_ */
