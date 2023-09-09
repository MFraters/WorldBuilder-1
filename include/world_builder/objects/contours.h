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

#ifndef WORLD_BUILDER_OBJECTS_CONTOURS_H
#define WORLD_BUILDER_OBJECTS_CONTOURS_H

//#include "world_builder/utilities.h"

#include <memory>
#include <vector>

#include "world_builder/nan.h"
#include "world_builder/objects/bezier_curve.h"
#include "world_builder/coordinate_systems/interface.h"
#include "world_builder/objects/natural_coordinate.h"
#include "world_builder/point.h"

namespace WorldBuilder
{
  namespace Objects
  {
    struct DistanceInterpolationData
    {
      DistanceInterpolationData() : signed_distance_from_bezier_surface(std::numeric_limits<double>::infinity()),
        distance_along_surface(std::numeric_limits<double>::infinity()),
        curve_above_index(NaN::DSNAN),
        curve_above_section_index(NaN::DSNAN),
        curve_above_interplation_fraction(NaN::DSNAN),
        curve_below_index(NaN::DSNAN),
        curve_below_section_index(NaN::DSNAN),
        curve_below_interplation_fraction(NaN::DSNAN),
        curve_local_interpolation_fraction(NaN::DSNAN) {}
      /**
       * @brief A distance to the Bezier surface. If only_positive is false, the distance is signed
       *        Up is negative and down is postive. It emulates a local depth system.
       */
      double signed_distance_from_bezier_surface;

      /**
       * @brief A distance along the Bezier surface to the first contour.
       */
      double distance_along_surface;

      size_t curve_above_index;
      size_t curve_above_section_index;
      double curve_above_interplation_fraction;
      size_t curve_below_index;
      size_t curve_below_section_index;
      double curve_below_interplation_fraction;
      double curve_local_interpolation_fraction;

      //double curve_above_interpolation_fraction;
      //double curve_below_interpolation_fraction;
      //double curve_local_interpolation_fraction;
      ///**
      // * @brief
      // *
      // */
      //size_t curve_above_index;
      //size_t curve_above_section_index;
      ////double curve_above_section_fraction;
      //size_t curve_below_index;
      //size_t curve_below_section_index;
      //double curve_below_section_fraction;
      ////size_t curve_local_index;
      ////size_t curve_local_section_index;
      //double curve_local_section_fraction;
    };

    template <class A, class B, class C>
    class Contours
    {
      public:
        /**
         * Constructor to create an empty contours.
         */
        Contours();

        /**
         * Constructor to create a contours from value at points object output.
         */
        Contours(const std::vector<std::vector<Point<2> > > &points,
                 const std::vector<double> depths,
                 const std::vector<std::vector<double> > &thicknesses,
                 const std::vector<std::vector<double> > &top_truncation,
                 const double start_radius,
                 const std::vector<std::vector<double> > &angle_contraints,
                 const std::vector<std::vector<Point<2> > > &directions,
                 std::vector<std::vector<std::vector<std::shared_ptr<A> > > > temperature_systems,
                 std::vector<std::vector<std::vector<std::shared_ptr<B> > > > composition_systems,
                 std::vector<std::vector<std::vector<std::shared_ptr<C> > > > grains_systems);

        /**
         * Copy constructor
         */
        Contours(Contours const &other);

        /**
         * Destructor
         */
        ~Contours();

        std::pair<Point<3>,Point<3> >
        compute_cross_section_axes(Point<3> origin,
                                   Point<3> x_direction,
                                   Point<3> y_direction) const;

        DistanceInterpolationData
        distance_interpolation_data(const Point<3> &check_point_cartesian,
                                    const Objects::NaturalCoordinate &check_point_natural,
                                    const Point<2> &reference_point,
                                    const std::unique_ptr<CoordinateSystems::Interface> &coordinate_system,
                                    const std::vector<std::vector<double> > &interpolation_properties,
                                    const double start_radius,
                                    const bool only_positive) const;

        /**
         * @brief Stores the splines for each depth contour
         *
         */
        std::vector<WorldBuilder::Objects::BezierCurve> contour_curves;

        std::vector<std::vector<Point<2> > > points; // TODO: probably do no need to store this one.
        std::vector<double> depths;
        std::vector<std::vector<double> > angle_contraints;
        std::vector<std::vector<double> > thicknesses;
        std::vector<std::vector<double> > top_truncation;
        std::vector<std::vector<Point<2> > > directions;
        double start_radius;

        std::vector<std::vector<std::vector<std::shared_ptr<A> > >> temperature_systems;
        std::vector<std::vector<std::vector<std::shared_ptr<B> > >> composition_systems;
        std::vector<std::vector<std::vector<std::shared_ptr<C> > >> grains_systems;

        /**
         * @brief Store the max thickness on each curve.
         *
         */
        std::vector<double> max_thickness_on_curve;

        /**
         * @brief stores the distance along the Bezier surface to reach this point.
         *
         */
        std::vector<std::vector<double> > distance_along_surface;
        /**
         * @brief stores the starting and ending depth range of each contour interval
         *
         * This means that the first entry is the depth[i]-thickness[i] and
         * the second entry is depth[i+1]+thickness[i+1]
         */
        //std::vector<std::pair<double,double> > depth_ranges;


        /**
         * @brief Stores which grid points a single grid point connects to the contour below it.
         *
         * Example layout: [
         * contour1:[connections for point 1:[0,1], connections for point 2: [2]],
         * contour2:[connections for point 1:[0], connections for point 2: [1,2]]
         * ]
         */
        //std::vector<std::vector<std::vector<unsigned int> > > connectivity;



      private:
    };


//    class Contours
//    {
//      public:
//        /**
//         * Constructor to create an empty contours.
//         */
//        Contours();
//
//        /**
//         * Constructor to create a contours from value at points object output.
//         */
//        Contours(std::vector<std::pair<double,std::vector<double>>> &contours);
//
//        /**
//         * Returns the value of the contours at the check point.
//         */
//        double local_value(const Point<2> check_point) const;
//
//        /**
//         * Wether the contours is a constant value or not. This is used for optimalization.
//         */
//        bool constant_value;
//
//        /**
//         * The minimum value of all provided points.
//         */
//        double minimum;
//
//        /**
//         * The maximum value of all provided points.
//         */
//        double maximum;
//
//        /**
//         * The KD tree which stores the centroids of all triangles and an index to the triangle points
//         * and values stored in the triangles member variable.
//         */
//        KDTree::KDTree tree;
//
//        /**
//         * Stores the triangles as a list of three points.
//         */
//        std::vector<std::array<std::array<double,3>,3> > triangles;
//
//        /**
//         * @brief A structure to store information on each element of the grid.
//         *
//         */
//        struct GridPoint
//        {
//          double dip;
//          double min_distance_along_slab;
//        };
//
//        /**
//         * @brief A variable to store the information on the grid.
//         *
//         */
//        std::vector<std::vector<GridPoint> > grid;
//
//        /**
//         * @brief Stores which grid points a single grid point connects to the contour below it.
//         *
//         * Example layout: [
//         * contour1:[connections for point 1:[0,1], connections for point 2: [2]],
//         * contour2:[connections for point 1:[0], connections for point 2: [1,2]]
//         * ]
//         */
//        std::vector<std::vector<std::vector<unsigned int> > > connectivity;
//
//        /**
//         * @brief Stores the splines for each depth contour
//         *
//         */
//        std::vector<std::array<Utilities::interpolation,2> > contour_splines;
//
//      private:
//        /**
//         * Test whether a point is in a triangle. If that is the case is stores the interpolated
//         * value of the tirangle into `interpolated_value` and returns true.
//         */
//        bool in_triangle(const std::array<std::array<double,3>,3> &points,
//                         const Point<2> check_point,
//                         double &interpolate_value) const;
//    };
  }

}

#endif