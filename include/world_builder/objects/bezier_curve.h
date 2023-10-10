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

#ifndef WORLD_BUILDER_OBJECTS_BEZIER_CURVE_H
#define WORLD_BUILDER_OBJECTS_BEZIER_CURVE_H

#include "world_builder/objects/closest_point_on_curve.h"
#include "world_builder/point.h"
#include <array>
#include <vector>

namespace WorldBuilder
{
  namespace Objects
  {

    /**
     * @brief Class for circle line/spline, including interpolation on it
     *
     */
    class BezierCurve
    {
      public:
        /**
         * @brief Construct a new Bezier Curve object
         *
         * @param p
         * @param angle_constrains
         */
        BezierCurve() {};

        /**
         * @brief Construct a new Bezier Curve object
         *
         * @param p
         * @param angle_constrains
         */
        BezierCurve(const std::vector<Point<2> > &p, const std::vector<double> &angle_constrains = {});

        /**
         * @brief Finds the closest point on the curve. If the the closest point
         *        doesn't fall on the segment, return a point with x and y being nan.
         *
         * @param p
         * @return ClosestPointOnCurve
         */
        ClosestPointOnCurve closest_point_on_curve_segment(const Point<2> &p) const;



        /**
         * @brief // TODO
         *
         * @param i
         * @param x
         * @return Point<2>
         */
        inline
        double approximate_length(const size_t points_per_segment) const
        {
          double arc_length = 0;
          double parametric_distance = 1./(double)points_per_segment;
          for (size_t point_i = 0; point_i < points.size()-1; ++point_i)
            {
              Point<2> previous_point = points[point_i];
              for (size_t segment_point_i = 0; segment_point_i < points_per_segment; ++segment_point_i)
                {
                  // new point
                  const double t = (segment_point_i+1)*parametric_distance;
                  Point<2> current_point = operator()(point_i,t);
                  arc_length += previous_point.distance(current_point);
                  previous_point = current_point;
                }
            }
          return arc_length;
        }


        /**
         * @brief // TODO
         *
         * @param i
         * @param x
         * @return Point<2>
         */
        inline
        std::vector<std::vector<double> > compute_parameter_to_length_map(const size_t points_per_segment) const
        {
          double arc_length = 0;
          double parametric_distance = 1./(double)points_per_segment;
          std::vector<std::vector<double> > parameter_to_length_map = std::vector<std::vector<double> > (points.size()-1);
          for (size_t point_i = 0; point_i < points.size()-1; ++point_i)
            {
              Point<2> previous_point = points[point_i];
              parameter_to_length_map[point_i].emplace_back(0);
              for (size_t segment_point_i = 0; segment_point_i < points_per_segment; ++segment_point_i)
                {
                  // new point
                  const double t = (segment_point_i+1)*parametric_distance;
                  Point<2> current_point = operator()(point_i,t);
                  arc_length += previous_point.distance(current_point);
                  //std::cout << "i = " << segment_point_i << ", t = " << t << ", arc_length = " << arc_length << std::endl;
                  parameter_to_length_map[point_i].emplace_back(arc_length);
                  previous_point = current_point;
                }
            }
          return parameter_to_length_map;
        }


        /**
         * @brief return the parameteric distance for a given arclength distance
         *
         * @param parameter_to_length_map
         * @param arc_length
         * @param distance
         * @return double Parametric distance (t)
         */
        inline
        double distance_to_t(const std::vector<std::vector<double> > &parameter_to_length_map,
                             const double distance) const
        {
          const double arc_length = parameter_to_length_map[parameter_to_length_map.size()-1][parameter_to_length_map[0].size()-1];
          const double n_points = parameter_to_length_map.size();
          const double points_per_segment = parameter_to_length_map[0].size();
          //std::cout << "distance = " << distance << ", arc_length = " << arc_length << std::endl;
          if (distance >= 0. && distance <= arc_length)
            {
              //std::cout << "points_per_segment = " << points_per_segment << ", n_points = " << n_points << std::endl;
              for (size_t point_i = 0; point_i < n_points; ++point_i)
                {
                  //std::cout << "point = " << point_i << ", parameter_to_length_map[point_i][0] = " << parameter_to_length_map[point_i][0] << ", parameter_to_length_map[point_i][points_per_segment-1] = " << parameter_to_length_map[point_i][points_per_segment-1] << std::endl;
                  if (distance >= parameter_to_length_map[point_i][0] && distance <= parameter_to_length_map[point_i][points_per_segment-1])
                    {
                      //std::cout << "flag 1" << std::endl;
                      for (size_t i = 0; i < points_per_segment-1; ++i)
                        {
                          //std::cout << "i = " << i << ", distance = " << distance << ", parameter_to_length_map[i] = " << parameter_to_length_map[point_i][i] << ", parameter_to_length_map[i+1] = " << parameter_to_length_map[point_i][i+1] << std::endl;
                          if (distance > parameter_to_length_map[point_i][i] && distance <= parameter_to_length_map[point_i][i+1])
                            {
                              const double previous_distance = parameter_to_length_map[point_i][i];
                              const double next_distance = parameter_to_length_map[point_i][i+1];
                              const double distance_fraction = (distance-previous_distance)/(next_distance-previous_distance);
                              //std::cout << "distance_fraction = " << distance_fraction << ", next_distance = " << next_distance << ", distance = " << distance << ", previous_distance = " << previous_distance << std::endl;
                              //const double previous_t_value = parameter_to_length_map[i][0];
                              //const double next_t_value = parameter_to_length_map[i+1][0];
                              const double previous_t_value = (double)i / double(points_per_segment-1);
                              const double next_t_value = (double)(i+1) / double(points_per_segment-1);
                              //std::cout << "parameter_to_length_map[i] = " << parameter_to_length_map[point_i][i] << ", parameter_to_length_map[i+1] = " << parameter_to_length_map[point_i][i+1]
                              //          << ", distance_fraction = " << distance_fraction << ", previous_t_value = " << previous_t_value << ", next_t_value = " << next_t_value
                              //          << ", return = " << previous_t_value + distance_fraction * (next_t_value-previous_t_value)<< std::endl;
                              return (double)point_i + previous_t_value + distance_fraction * (next_t_value-previous_t_value);
                            }
                        }
                    }
                }
            }
          return distance/arc_length + (distance > arc_length ? (n_points-1) : 0.0);
        }

        /**
         * @brief // TODO
         *
         * @param i
         * @param x
         * @return Point<2>
         */
        Point<2> operator()(const size_t i, const double x) const;

        /**
         * @brief // TODO
         *
         */
        inline
        const std::vector<Point<2> > &get_points() const
        {
          return points;
        }
        /**
         * @brief // TODO
         *
         */
        inline
        const std::vector<Point<2> > &get_points()
        {
          return points;
        }
        /**
         * @brief // TODO
         *
         */
        inline
        const std::vector<std::array<Point<2>,2 > > &get_control_points() const
        {
          return control_points;
        }
        /**
         * @brief // TODO
         *
         */
        inline
        const std::vector<std::array<Point<2>,2 > > &get_control_points()
        {
          return control_points;
        }

        std::vector<double> angles;
      private:
        std::vector<Point<2> > points;
        std::vector<std::array<Point<2>,2 > > control_points;
        std::vector<double> lengths;

    };
  }

}


#endif