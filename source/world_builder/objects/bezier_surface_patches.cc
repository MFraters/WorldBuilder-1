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
#include "world_builder/objects/bezier_surface_patches.h"
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
    BezierSurfacePatches::BezierSurfacePatches(std::vector<WorldBuilder::Objects::BezierCurve> contour_curves,
                                               std::vector<double> depths,
                                               std::vector<std::vector<double> > downward_angle_contraints,
                                               std::vector<std::vector<double> > thicknesses,
                                               std::vector<std::vector<double> > top_truncation,
                                               std::vector<std::vector<Point<2> > > directions)
    {
      constexpr double degree_to_rad = Consts::PI/180.0;

      // find the contour which has the most points
      // also compute the length of each curve
      double contour_most_points = 0;
      for (const Objects::BezierCurve &contour_i : contour_curves)
        {
          if (contour_i.get_points().size()>contour_most_points)
            {
              contour_most_points = contour_i.get_points().size();
            }
        }
      WBAssertThrow(contour_most_points > 0, "Contours do not have enough points.");

      // create new contours with this many points, distributed more or less evenly
      std::vector<WorldBuilder::Objects::BezierCurve> evenly_distributed_contour_curves;
      std::vector<std::vector<double> > evenly_distributed_angle_contraints;
      for (size_t contour_i = 0; contour_i < contour_curves.size(); ++contour_i)
        {
          // compute the length and create t to length map
          const std::vector<std::vector<double> > parameter_to_length_map = contour_curves[contour_i].compute_parameter_to_length_map(5);
          const double arc_length = parameter_to_length_map[parameter_to_length_map.size()-1][parameter_to_length_map[0].size()-1];

          // create new contours
          std::vector<Point<2>> contour_points;
          contour_points.reserve(contour_most_points);
          std::vector<double> contour_angles;
          contour_angles.reserve(contour_most_points);
          const double segment_length = arc_length/(contour_most_points-1);
          for (size_t point_i = 0; point_i < contour_most_points; ++point_i)
            {
              const double distance = segment_length*point_i;
              const double parameter_t = contour_curves[contour_i].distance_to_t(parameter_to_length_map, distance);
              //std::cout << "distance = " << distance << ", parameter_t = " << parameter_t << std::endl;
              //auto point = contour_curves[contour_i](parameter_t,parameter_t-(size_t)parameter_t);
              //std::cout << "point = " << point << ", contour_most_points = " << contour_most_points << ",segment_length = " << segment_length << std::endl;
              contour_points.emplace_back(contour_curves[contour_i](parameter_t,parameter_t-(size_t)parameter_t));
              contour_angles.emplace_back(downward_angle_contraints[contour_i][(size_t)parameter_t] + (parameter_t-(size_t)parameter_t) * (downward_angle_contraints[contour_i][(size_t)parameter_t+1]-downward_angle_contraints[contour_i][(size_t)parameter_t]));
              std::cout << "contour_points = " << contour_points[contour_points.size()-1] << ", contour_angles = " << contour_angles[contour_points.size()-1] << std::endl;
            }
          evenly_distributed_contour_curves.emplace_back(Objects::BezierCurve(contour_points));
          evenly_distributed_angle_contraints.emplace_back(contour_angles);
        }

      // find the vector pointing to next point down
      patches.resize(evenly_distributed_contour_curves.size()-1);
      for (size_t contour_i = 0; contour_i < evenly_distributed_contour_curves.size()-1; ++contour_i)
        {
          patches[contour_i].resize(contour_most_points-1);
          std::array<std::array<Point<3>,4>,4> patch_points;
          for (size_t point_i = 0; point_i < contour_most_points-1; ++point_i)
            {
              std::cout << "------------------------------" << point_i << "---------------------------------" << std::endl;
              std::array<Point<2>, 4> top_points =
              {
                {
                  evenly_distributed_contour_curves[contour_i].get_points()[point_i],
                  evenly_distributed_contour_curves[contour_i].get_control_points()[point_i][0],
                  evenly_distributed_contour_curves[contour_i].get_control_points()[point_i][1],
                  evenly_distributed_contour_curves[contour_i].get_points()[point_i+1]
                }
              };

              std::array<double,4> top_angles =
              {
                {
                  evenly_distributed_angle_contraints[contour_i][point_i],
                  evenly_distributed_angle_contraints[contour_i][point_i]+1./3.*(evenly_distributed_angle_contraints[contour_i][point_i+1]-evenly_distributed_angle_contraints[contour_i][point_i]),
                  evenly_distributed_angle_contraints[contour_i][point_i]+2./3.*(evenly_distributed_angle_contraints[contour_i][point_i+1]-evenly_distributed_angle_contraints[contour_i][point_i]),
                  evenly_distributed_angle_contraints[contour_i][point_i+1]
                }
              };

              std::cout << "top angle 2 = " << evenly_distributed_angle_contraints[contour_i][point_i]+1./3.*(evenly_distributed_angle_contraints[contour_i][point_i+1]-evenly_distributed_angle_contraints[contour_i][point_i])
                        << ", 1:" << evenly_distributed_angle_contraints[contour_i][point_i] << "+ 1/3.*(" << evenly_distributed_angle_contraints[contour_i][point_i+1] << "-" << evenly_distributed_angle_contraints[contour_i][point_i] << ")" << std::endl;

              std::cout << "top angle 3 = " << evenly_distributed_angle_contraints[contour_i][point_i]+2./3.*(evenly_distributed_angle_contraints[contour_i][point_i+1]-evenly_distributed_angle_contraints[contour_i][point_i])
                        << ", 1:" << evenly_distributed_angle_contraints[contour_i][point_i] << "+ 1/3.*(" << evenly_distributed_angle_contraints[contour_i][point_i+1] << "-" << evenly_distributed_angle_contraints[contour_i][point_i] << ")" << std::endl;

              std::array<Point<2>, 4> bottom_points =
              {
                {
                  evenly_distributed_contour_curves[contour_i+1].get_points()[point_i],
                  evenly_distributed_contour_curves[contour_i+1].get_control_points()[point_i][0],
                  evenly_distributed_contour_curves[contour_i+1].get_control_points()[point_i][1],
                  evenly_distributed_contour_curves[contour_i+1].get_points()[point_i+1]
                }
              };

              std::array<double,4> bottom_angles =
              {
                {
                  evenly_distributed_angle_contraints[contour_i+1][point_i],
                  evenly_distributed_angle_contraints[contour_i+1][point_i]+1./3.*(evenly_distributed_angle_contraints[contour_i+1][point_i+1]-evenly_distributed_angle_contraints[contour_i+1][point_i]),
                  evenly_distributed_angle_contraints[contour_i+1][point_i]+2./3.*(evenly_distributed_angle_contraints[contour_i+1][point_i+1]-evenly_distributed_angle_contraints[contour_i+1][point_i]),
                  evenly_distributed_angle_contraints[contour_i+1][point_i+1]
                }
              };

              std::cout << "bottom angle 2 = " << evenly_distributed_angle_contraints[contour_i+1][point_i]+1./3.*(evenly_distributed_angle_contraints[contour_i+1][point_i+1]-evenly_distributed_angle_contraints[contour_i+1][point_i])
                        << ", 1:" << evenly_distributed_angle_contraints[contour_i+1][point_i] << "+ 1/3.*(" << evenly_distributed_angle_contraints[contour_i+1][point_i+1] << "-" << evenly_distributed_angle_contraints[contour_i+1][point_i] << ")" << std::endl;

              std::cout << "bottom angle 3 = " << evenly_distributed_angle_contraints[contour_i+1][point_i]+2./3.*(evenly_distributed_angle_contraints[contour_i+1][point_i+1]-evenly_distributed_angle_contraints[contour_i+1][point_i])
                        << ", 1:" << evenly_distributed_angle_contraints[contour_i+1][point_i] << "+ 1/3.*(" << evenly_distributed_angle_contraints[contour_i+1][point_i+1] << "-" << evenly_distributed_angle_contraints[contour_i+1][point_i] << ")" << std::endl;


              for (size_t i = 0; i < 4; ++i)
                {
                  // one downard strip of patch
                  const Point<3> up_point = Point<3>(top_points[i][0],top_points[i][1],depths[contour_i],cartesian);
                  const Point<3> down_point = Point<3>(bottom_points[i][0],bottom_points[i][1],depths[contour_i+1],cartesian);

                  patch_points[0][i] = up_point;
                  patch_points[3][i] = down_point;

                  const double simple_3d_length = (down_point-up_point).norm();

                  {
                    // compute the upper control point

                    std::cout << "===================a=======================\nflag 1: i="<< i << ", down_point_2d-up_point_2d = "
                              << bottom_points[i]-top_points[i] << ", bottom_points[i] = " << bottom_points[i] << ", top_points[i] = " << top_points[i]
                              << ", angle = " << top_angles[i] << std::endl;
                    // TODO: This next_point_2d vector should probably be replaced with the normal to the curve at this point.
                    const Point<2> next_point_2d = (bottom_points[i]-top_points[i])/(bottom_points[i]-top_points[i]).norm();
                    const Point<3> first_control_point_3d_unrotated = Point<3>(next_point_2d[0]*simple_3d_length*0.1,next_point_2d[1]*simple_3d_length*0.1,0,cartesian);
                    // the vector_first_control_point_3d is the control piont for 0 degree dip.
                    // Now rotate with the needed angle.
                    // TODO: implement spherical by changing the depth in spherical coordinates and then
                    // convert to cartesian and subtract the point by the point which has a changed depth.
                    Point<3> y_axis = Point<3>(0,0,1,cartesian);
                    Point<3> x_axis = first_control_point_3d_unrotated/first_control_point_3d_unrotated.norm();
                    const double angle = top_angles[i]*degree_to_rad;

                    // compute location of point in 2D plane
                    const Point<2> first_control_point_2d_origin(x_axis * first_control_point_3d_unrotated,
                                                                 y_axis * first_control_point_3d_unrotated,
                                                                 cartesian);
                    // rotate 2D point with the provided angle.
                    const Point<2> first_control_point_2d_origin_rotated = Point<2>(
                                                                             first_control_point_2d_origin[0]*cos(angle)-first_control_point_2d_origin[1]*sin(angle),
                                                                             first_control_point_2d_origin[0]*sin(angle)-first_control_point_2d_origin[1]*cos(angle),
                                                                             cartesian);
                    // now move it back to 3d
                    const Point<3> first_control_point_3d = up_point + x_axis*first_control_point_2d_origin_rotated[0] + y_axis*first_control_point_2d_origin_rotated[1];


                    std::cout << "control_point_3d_unrotated = " << first_control_point_3d_unrotated
                              << ", control_point_2d_origin = " << first_control_point_2d_origin
                              <<  ", control_point_2d_origin_rotated = " << first_control_point_2d_origin_rotated
                              << ", control_point_3d = " << first_control_point_3d
                              << ", x:y axis = " << x_axis << ":" << y_axis << std::endl;
                    patch_points[1][i] = first_control_point_3d;
                  }

                  {
                    // compute the bottom control point

                    std::cout << "====================b======================\nflag 1: up_point_2d-down_point_2d = " << top_points[i]-bottom_points[i]
                              << ", bottom_points[i] = " << bottom_points[i] << ", top_points[i] = " << top_points[i]
                              << ", angle = " << bottom_angles[i] << std::endl;
                    // TODO: This next_point_2d vector should probably be replaced with the normal to the curve at this point.
                    const Point<2> next_point_2d = (top_points[i]-bottom_points[i])/(top_points[i]-bottom_points[i]).norm();
                    const Point<3> second_control_point_3d_unrotated = Point<3>(next_point_2d[0]*simple_3d_length*0.1,next_point_2d[1]*simple_3d_length*0.1,0,cartesian);
                    // the vector_first_control_point_3d is the control piont for 0 degree dip.
                    // Now rotate with the needed angle.
                    // TODO: implement spherical by changing the depth in spherical coordinates and then
                    // convert to cartesian and subtract the point by the point which has a changed depth.
                    Point<3> y_axis = Point<3>(0,0,1,cartesian);
                    Point<3> x_axis = second_control_point_3d_unrotated/second_control_point_3d_unrotated.norm();
                    const double angle = -bottom_angles[i]*degree_to_rad;

                    // compute location of point in 2D plane
                    const Point<2> second_control_point_2d_origin(x_axis * second_control_point_3d_unrotated,
                                                                  y_axis * second_control_point_3d_unrotated,
                                                                  cartesian);
                    // rotate 2D point with the provided angle.
                    const Point<2> second_control_point_2d_origin_rotated = Point<2>(
                                                                              second_control_point_2d_origin[0]*cos(angle)-second_control_point_2d_origin[1]*sin(angle),
                                                                              second_control_point_2d_origin[0]*sin(angle)-second_control_point_2d_origin[1]*cos(angle),
                                                                              cartesian);
                    // now move it back to 3d
                    const Point<3> second_control_point_3d = down_point + x_axis*second_control_point_2d_origin_rotated[0] + y_axis*second_control_point_2d_origin_rotated[1];


                    std::cout << "control_point_3d_unrotated = " << second_control_point_3d_unrotated
                              << ", control_point_2d_origin = " << second_control_point_2d_origin
                              <<  ", control_point_2d_origin_rotated = " << second_control_point_2d_origin_rotated
                              << ", control_point_3d = " << second_control_point_3d
                              << ", x:y axis = " << x_axis << ":" << y_axis << std::endl;
                    patch_points[2][i] = second_control_point_3d;
                  }
                }

              //for (size_t i = 0; i < 4; ++i)
              //  {
              //    for (size_t j = 0; j < 4; ++j)
              //      {
              //        WBAssert(patches[contour_i][point_i].control_points[i][j].get_coordinate_system() != invalid, "coord is not set for " << i << ":" << j);
              //        std::cout << patches[contour_i][point_i].control_points[i][j] << " | ";
              //      }
              //    std::cout << std::endl;
              //  }
            }
          patches[contour_i].emplace_back(Objects::BezierSurfacePatch(patch_points));
        }
    }


    Point<3> BezierSurfacePatches::point_on_surface(const size_t i, const size_t j, const double u, const double v) const
    {

      return patches[i][j].point_on_surface(u,v);
    }

    Point<3> BezierSurfacePatches::tangent_point(const size_t i, const size_t j, const double u, const double v)
    {
      return patches[i][j].tangent_point(u,v);
    }

    ClosestPointOnCurve BezierSurfacePatches::closest_point_on_patch(const Point<3> &p) const
    {
      double smallest_distance = std::numeric_limits<double>::infinity();
      ClosestPointOnCurve closest_closest_point;
      for (size_t i =0; i < patches.size(); ++i)
        {
          for (size_t j = 0; j < patches[i].size(); ++i)
            {
              const ClosestPointOnCurve closest_point = patches[i][j].closest_point_on_patch(p);
              if (closest_point.distance < smallest_distance)
                {
                  closest_closest_point =  closest_point;
                }
            }
        }

      WBAssert(smallest_distance < std::numeric_limits<double>::infinity(), "clost point not found.");

      return closest_closest_point;
    }

  }
}