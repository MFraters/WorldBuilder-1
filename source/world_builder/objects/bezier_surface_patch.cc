/*
Copyright (C) 2018 - 2023 by the authors of the World Builder code.

This file is part of the World Builder.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "world_builder/assert.h"
#include "world_builder/coordinate_system.h"
#include "world_builder/nan.h"
#include "world_builder/objects/bezier_curve.h"
#include "world_builder/objects/bezier_surface_patch.h"
#include "world_builder/objects/closest_point_on_curve.h"
#include "world_builder/types/double.h"
#include "world_builder/utilities.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace WorldBuilder;

namespace WorldBuilder
{
  namespace Objects
  {
    BezierSurfacePatch::BezierSurfacePatch(std::array<std::array<Point<3>,4>,4> &p)
      :
      control_points(p)
    {
      double u = 0;
      double v = 0;
      size_t counter = 0;
      check_grid_points.resize(grid_points*grid_points,Point<3>(cartesian));
      check_grid_tangents.resize(grid_points*grid_points,Point<3>(cartesian));
      check_grid_points_u.resize(grid_points*grid_points);
      check_grid_points_v.resize(grid_points*grid_points);
      for (size_t i = grid_points; i--;)
        {
          for (size_t j = grid_points; j--;)
            {
              // TODO: all these point on surfaces could be precomputed, making this step really cheap.
              // That may even allow for a finer grid.
              check_grid_points[counter] = this->point_on_surface(u,v);
              check_grid_tangents[counter] = this->tangent_point(u, v);
              check_grid_points_u[counter] = u;
              check_grid_points_v[counter] = v;
              v += grid_distance;
              counter++;
            }
          v = 0;
          u += grid_distance;
        }
    }


    Point<3> BezierSurfacePatch::point_on_surface(const double u, const double v) const
    {
      const std::array<double,4> B_u = {(1-u) *(1-u) *(1-u),3*u *u *u-6*u *u+3*u,-3*u *u *u+3*u*u,u *u*u};
      const std::array<double,4> B_v = {(1-v) *(1-v) *(1-v),3*v *v *v-6*v *u+3*v,-3*v *v *v+3*v*v,v *v*v};
      Point<3> result_point(cartesian);
      for (size_t i = 16; i--;)
        {
          for (size_t j = 4; j--;)
            {
              result_point[0] += control_points[i][j][0]*B_u[i]*B_v[j];
              result_point[1] += control_points[i][j][1]*B_u[i]*B_v[j];
              result_point[2] += control_points[i][j][2]*B_u[i]*B_v[j];
            }
        }
      return Point<3>(cartesian);
    }

    Point<3> BezierSurfacePatch::tangent_point(const double u, const double v)
    {
      // compute derivative
      const Point<3> u_min_point = point_on_surface(u-u*std::numeric_limits<double>::epsilon()*10.0, v);
      const Point<3> u_plus_point = point_on_surface(u+u*std::numeric_limits<double>::epsilon()*10.0, v);
      const Point<3> derivative_u = u_plus_point-u_min_point;
      const Point<3> v_min_point = point_on_surface(u,v-v*std::numeric_limits<double>::epsilon()*10.0);
      const Point<3> v_plus_point = point_on_surface(u,v+v*std::numeric_limits<double>::epsilon()*10.0);
      const Point<3> derivative_v = v_plus_point-v_min_point;

      // the normalized tangent
      return Utilities::cross_product(derivative_u/derivative_u.norm(), derivative_v/derivative_v.norm());
    }

    ClosestPointOnCurve BezierSurfacePatch::closest_point_on_patch(const Point<3> &p) const
    {
      // minimize |Pc-Px(u,v)|
      // 1. Find which points in a grid is closest
      double min_distance_squared = std::numeric_limits<double>::infinity();
      double min_u = 0.0;
      double min_v = 0.0;
      Point<3> min_point = Point<3>(NaN::DSNAN,NaN::DSNAN,NaN::DSNAN,invalid);
      double u = 0;
      double v = 0;
      size_t counter = 0;
      for (size_t i = grid_points*grid_points; i--;)
        {
          // TODO: all these point on surfaces could be precomputed, making this step really cheap.
          // That may even allow for a finer grid.
          const double min_distance_squared_tmp = (p-check_grid_points[counter]).norm_square();
          if (min_distance_squared_tmp < min_distance_squared)
            {
              min_distance_squared = min_distance_squared_tmp;
              min_point = check_grid_points[counter];
              min_u = check_grid_points_u[counter];
              min_v = check_grid_points_v[counter];
            }
        }
      // 2. now use the minimum point from the grid search to search and iterate
      // TODO: maybe? Or just have more precomputed points?


      return ClosestPointOnCurve();
    }
  }
}