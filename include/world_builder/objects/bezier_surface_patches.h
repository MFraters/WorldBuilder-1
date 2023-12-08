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

#ifndef WORLD_BUILDER_OBJECTS_BEZIER_SURFACE_PATCHES_H
#define WORLD_BUILDER_OBJECTS_BEZIER_SURFACE_PATCHES_H

#include "world_builder/objects/closest_point_on_curve.h"
#include "world_builder/objects/bezier_surface_patch.h"
#include <array>
#include <vector>

namespace WorldBuilder
{
  template <int dim> class Point;

  namespace Objects
  {
    class BezierCurve;

    /**
     * @brief Class for a Bezier surface patch
     *
     */
    class BezierSurfacePatches
    {
      public:

        /**
         * // TODO
         */
        BezierSurfacePatches() {};
        /**
         * // TODO
         */
        BezierSurfacePatches(std::vector<WorldBuilder::Objects::BezierCurve> contour_curves,
                             std::vector<double> depths,
                             std::vector<std::vector<double> > downward_angle_contraints,
                             std::vector<std::vector<double> > thicknesses,
                             std::vector<std::vector<double> > top_truncation,
                             std::vector<std::vector<Point<2> > > directions);


        /**
         * // TODO
         */
        Point<3> point_on_surface(const size_t i, const size_t j, const double u, const double v) const;


        /**
         * // TODO
         */
        Point<3> tangent_point(const size_t i, const size_t j, const double u, const double v);


        /**
         * @brief Finds the closest point on the curve. If the the closest point
         *        doesn't fall on the segment, return a point with x and y being nan.
         *
         * @param p
         * @return ClosestPointOnCurve
         */
        ClosestPointOnCurve closest_point_on_patch(const Point<3> &p) const;

        /**
         * @brief
         *
         * @param i
         * @param x
         * @return Point<2>
         */
        Point<2> operator()(const size_t i, const double x) const;

        std::vector<std::vector<Objects::BezierSurfacePatch>> patches;

      private:
        std::array<std::array<Point<3>,4>,4> control_points;
        std::vector<Point<3>> check_grid_points;
        std::vector<Point<3>> check_grid_tangents;
        std::vector<double> check_grid_points_u;
        std::vector<double> check_grid_points_v;

    };
  }

}


#endif