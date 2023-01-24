/*
Copyright (C) 2018 - 2021 by the authors of the World Builder code.

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
#include "world_builder/objects/bezier_curve.h"
#include "world_builder/objects/cubic_spline.h"
#include "world_builder/assert.h"
#include "world_builder/nan.h"
#include "world_builder/utilities.h"

#include <sstream>
#include <cstddef>
#include <iomanip>

using namespace WorldBuilder;

namespace WorldBuilder
{
  namespace Objects
  {
    BezierCurve::BezierCurve(const std::vector<Point<2> > &p, const std::vector<double> &angle_constrains_input)
    {
      points = p;
      const size_t n_points = p.size();
      control_points.resize(n_points-1, {p[0],p[0]});
      lengths.resize(n_points-1,NaN::DSNAN);
      angles.resize(n_points,NaN::DSNAN);
      std::vector<double> angle_constrains = angle_constrains_input;
      angle_constrains.resize(n_points,NaN::DSNAN);

      // if no angle is provided, compute the angle as the average angle bewteen the previous and next point.
      // The first angle points at the second point and the last angle points at the second to last point.
      // The check points are set at a distance of 1/10th the line length from the point in the direction of the angle.
      if (std::isnan(angle_constrains[0]))
        {
          Point<2> P1P2 = points[1]-points[0];
          angles[0] = atan2(P1P2[1],P1P2[0]);
        }
      else
        {
          angles[0] = angle_constrains[0];
        }

      for (size_t p_i = 1; p_i < n_points-1; ++p_i)
        {
          bool clockwise = true;

          // first determine the angle
          if (std::isnan(angle_constrains[p_i]))
            {
              // get the average angle
              const Point<2> P1P2 = points[p_i-1]-points[p_i];
              const Point<2> P3P2 = points[p_i+1]-points[p_i];

              double diff_angle = atan2(P1P2[0]*P3P2[1] - P1P2[1]*P3P2[0] , P1P2[0]*P3P2[0] + P1P2[1]*P3P2[1]);
              const double angle_p1p2 = atan2(P1P2[1],P1P2[0]);
              const double angle_p3p1 = atan2(P3P2[1],P3P2[0]);
              const double average_angle = (angle_p1p2 + angle_p3p1)*0.5;
              angles[p_i] = average_angle;

              // if direction changes between clockwise and anti-clockwise
              // adjust direction by 90 degrees.
              //  if (diff_angle < 0 && !clockwise)
              //    {
              //      clockwise=true;
              //      angles[p_i] -= Consts::PI*0.5;
              //    }
              //  else if (diff_angle > 0 && clockwise)
              //    {
              //      clockwise=false;
              //      angles[p_i] -= Consts::PI*0.5;
              //    }
              angles[p_i] -= Consts::PI*0.5;
            }
          else
            {
              angles[p_i] = angle_constrains[p_i];
            }

        }

      if (std::isnan(angle_constrains[n_points-1]))
        {
          Point<2> P1P2 = points[n_points-2]-points[n_points-1];
          Point<2> P2P1 = points[n_points-1]-points[n_points-2];
          angles[n_points-1] =  atan2(P1P2[1],P1P2[0]);

          //std::cout << "atan2(P1P2[0],P1P2[1]) degree = " << atan2(P1P2[0],P1P2[1])*(180./Consts::PI) << ", atan2(P1P2[1],P1P2[0]) degree = " << atan2(P1P2[1],P1P2[0])*(180./Consts::PI) << std::endl;
          //std::cout << "atan2(P2P1[0],P2P1[1]) degree = " << atan2(P2P1[0],P2P1[1])*(180./Consts::PI) << ", atan2(P2P1[1],P2P1[0]) degree = " << atan2(P2P1[1],P2P1[0])*(180./Consts::PI) << std::endl;
          //std::cout << "angles[n_points-1]=" << angles[n_points-1]*(180./Consts::PI) << std::endl;;
        }
      else
        {
          angles[n_points-1] = angle_constrains[n_points-1];
        }


      // next determine the location of the control points
      // the location of the control point is 1/10th p1p2 distance in the direction of the angle.
      // make sure the angle is pointing away from the next point, e.g.
      // the check point is on the other side of the of the line p1p2 compared to p3.
      const double fraction_of_length = 0.2;
      {
        const Point<2> &p1 = points[0];
        const Point<2> &p2 = points[1];
        const Point<2> &p3 = points[2];
        const double length = (points[0]-points[1]).norm(); // can be squared
        control_points[0][0][0] = cos(angles[0])*length*fraction_of_length+p1[0];
        control_points[0][0][1] = sin(angles[0])*length*fraction_of_length+p1[1];
        control_points[0][1][0] = cos(angles[1])*length*fraction_of_length+p2[0];
        control_points[0][1][1] = sin(angles[1])*length*fraction_of_length+p2[1];
        {
          const int side_of_line_1 =  (p1[0] - p2[0]) * (control_points[0][1][1] - p1[1])
                                      - (p1[1] - p2[1]) * (control_points[0][1][0] - p1[0])
                                      < 0 ? -1.0 : 1.0;
          const int side_of_line_2 =  (p1[0] - p2[0]) * (p3[1] - p1[1])
                                      - (p1[1] - p2[1]) * (p3[0] - p1[0])
                                      < 0 ? -1.0 : 1.0;
          if (side_of_line_1 == side_of_line_2)
            {
              //std::cout << "switching a" << std::endl;
              // use a 180 degree rotated angle to create this control_point
              control_points[0][1][0] = cos(angles[1]+Consts::PI)*length*fraction_of_length+p2[0];
              control_points[0][1][1] = sin(angles[1]+Consts::PI)*length*fraction_of_length+p2[1];
            }
        }
      }
      Point<2> a = 3.*control_points[0][0]-3.*control_points[0][1]+points[0+1]-points[0];
      Point<2> b = 3.*points[0] - 6.*control_points[0][0]+3.*control_points[0][1];
      Point<2> c = -3.*points[0] + 3.*control_points[0][0];
      Point<2> d = points[0];
      //std::cout << "parametric plot (" << a[0]*(180/Consts::PI) << "*t^3 + " << b[0]*(180/Consts::PI) << "*t^2 + " << c[0]*(180/Consts::PI) << "*t + " << d[0]*(180/Consts::PI) << ","
      //              << a[1]*(180/Consts::PI) << "*t^3 + " << b[1]*(180/Consts::PI) << "*t^2 + " << c[1]*(180/Consts::PI) << "*t + " << d[1]*(180/Consts::PI) << ") for t=0 to 1" <<std::endl;

      //std::cout << "parametric plot (" << a[0] << "*t^3 + " << b[0] << "*t^2 + " << c[0] << "*t + " << d[0] << ","
      //          << a[1] << "*t^3 + " << b[1] << "*t^2 + " << c[1] << "*t + " << d[1] << ") for t=0 to 1" <<std::endl;

      for (size_t p_i = 1; p_i < n_points-1; ++p_i)
        {
          const Point<2> &p1 = points[p_i];
          const Point<2> &p2 = points[p_i+1];
          const Point<2> &p3 = points[p_i+2];
          const double length = (points[p_i]-points[p_i+1]).norm(); // can be squared
          control_points[p_i][0][0] = cos(angles[p_i])*length*fraction_of_length+p1[0];
          control_points[p_i][0][1] = sin(angles[p_i])*length*fraction_of_length+p1[1];

          {
            const int side_of_line_1 =  (p1[0] - p2[0]) * (control_points[p_i-1][1][1] - p1[1])
                                        - (p1[1] - p2[1]) * (control_points[p_i-1][1][0] - p1[0])
                                        < 0 ? -1.0 : 1.0;
            const int side_of_line_2 =  (p1[0] - p2[0]) * (control_points[p_i][0][1] - p1[1])
                                        - (p1[1] - p2[1]) * (control_points[p_i][0][0] - p1[0])
                                        < 0 ? -1.0 : 1.0;
            if (side_of_line_1 == side_of_line_2)
              {
                //std::cout << "switching b: old = " << control_points[p_i][0] << std::endl;
                // use a 180 degree rotated angle to create this control_point
                control_points[p_i][0][0] = cos(angles[p_i]+Consts::PI)*length*fraction_of_length+p1[0];
                control_points[p_i][0][1] = sin(angles[p_i]+Consts::PI)*length*fraction_of_length+p1[1];
              }
          }

          control_points[p_i][1][0] = cos(angles[p_i+1])*length*fraction_of_length+points[p_i+1][0];
          control_points[p_i][1][1] = sin(angles[p_i+1])*length*fraction_of_length+points[p_i+1][1];
          //std::cout << "angles[p_i+1]= " << angles[p_i+1]*(180./Consts::PI) << std::endl;
          if (p_i+1 < n_points-1)
            {
              const int side_of_line_1 =  (p1[0] - p2[0]) * (control_points[p_i][1][1] - p1[1])
                                          - (p1[1] - p2[1]) * (control_points[p_i][1][0] - p1[0])
                                          < 0 ? -1.0 : 1.0;
              const int side_of_line_2 =  (p1[0] - p2[0]) * (p3[1] - p1[1])
                                          - (p1[1] - p2[1]) * (p3[0] - p1[0])
                                          < 0 ? -1.0 : 1.0;
              if (side_of_line_1 == side_of_line_2)
                {
                  // use a 180 degree rotated angle to create this control_point
                  //std::cout << "switching c: old = " << control_points[p_i][1] << std::endl;
                  control_points[p_i][1][0] = cos(angles[p_i+1]+Consts::PI)*length*fraction_of_length+p2[0];
                  control_points[p_i][1][1] = sin(angles[p_i+1]+Consts::PI)*length*fraction_of_length+p2[1];
                }
            }

          Point<2> a = 3.*control_points[p_i][0]-3.*control_points[p_i][1]+points[p_i+1]-points[p_i];
          Point<2> b = 3.*points[p_i] - 6.*control_points[p_i][0]+3.*control_points[p_i][1];
          Point<2> c = -3.*points[p_i] + 3.*control_points[p_i][0];
          Point<2> d = points[p_i];
          //std::cout << "1: " << points[p_i] << ", 2: " << control_points[p_i][0] << ", 3: " << control_points[p_i][1] << ", 4: " << points[p_i+1] << std::endl;
          //std::cout << "parametric plot (" << a[0]*(180/Consts::PI) << "*t^3 + " << b[0]*(180/Consts::PI) << "*t^2 + " << c[0]*(180/Consts::PI) << "*t + " << d[0]*(180/Consts::PI) << ","
          //          << a[1]*(180/Consts::PI) << "*t^3 + " << b[1]*(180/Consts::PI) << "*t^2 + " << c[1]*(180/Consts::PI) << "*t + " << d[1]*(180/Consts::PI) << ") for t=0 to 1" <<std::endl;
          //std::cout << "parametric plot (" << a[0] << "*t^3 + " << b[0] << "*t^2 + " << c[0] << "*t + " << d[0] << ","
          //          << a[1] << "*t^3 + " << b[1] << "*t^2 + " << c[1] << "*t + " << d[1] << ") for t=0 to 1" <<std::endl;
        }
      //std::cout << "==============================" << std::endl;

    }


    Point<2>
    BezierCurve::operator()(const size_t i, const double t) const
    {
      return (1-t)*(1-t)*(1-t)*points[i] + 3*(1-t)*(1-t)*t*control_points[i][0] + 3.*(1-t)*t*t*control_points[i][1]+t*t*t*points[i+1]; //(1-t)*(1-t)*points[i] + 2*t*(1-t)*control_points[i] + t*t*points[i+1];
    }

    double
    BezierCurve::arc_length(const Point<2> &a, const Point<2>  &b, const Point<2> &c) const
    {
      // compute the arc length of the Bezier curve
      // see https://gamedev.stackexchange.com/questions/6009/bezier-curve-arc-length
      Point<2> v = a;
      Point<2> w = b;
      v[0] = 2.*(b[0] - a[0]);
      v[1] = 2.*(b[1] - a[1]);
      w[0] = c[0] - 2.*b[0] + a[0];
      w[1] = c[1] - 2.*b[1] + a[1];

      const double uu = 4.*(w[0]*w[0] + w[1]*w[1]);

      // check if all points are in a line
      if (uu < 0.00001 || (a-b).norm() < 1e-6 || (b-c).norm() < 1e-6 || fabs(a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])) <= 1e-9)
        {
          return std::sqrt((c[0] - a[0])*(c[0] - a[0]) + (c[1] - a[1])*(c[1] - a[1]));
        }

      const double vv = 4.*(v[0]*w[0] + v[1]*w[1]);
      const double ww = v[0]*v[0] + v[1]*v[1];

      const double t1 = (2.*std::sqrt(uu*(uu + vv + ww)));
      const double t2 = 2.*uu+vv;
      const double t3 = vv*vv - 4.*uu*ww;
      const double t4 = (2.*std::sqrt(uu*ww));

      WBAssert(!std::isnan(((t1*t2 - t3*std::log(t2+t1) -(vv*t4 - t3*std::log(vv+t4))) / (8*std::pow(uu, 1.5)))), "result is nan: " << ((t1*t2 - t3*std::log(t2+t1) -(vv*t4 - t3*std::log(vv+t4))) / (8*std::pow(uu, 1.5)))
               << ", t1=" << t1 << ", t2=" << t2 << ", t3=" << t3 << ", t4=" << t4 << ", v = " << v << ", w=" << w << ", vv=" << vv << ", uu=" << uu
               << ", a=" << a << ", b=" << b << ", c=" << c
               << ", std::log(t2+t1)=" << std::log(t2+t1) << ", std::log(vv+t4)=" << std::log(vv+t4) << ", std::pow(uu, 1.5)=" << std::pow(uu, 1.5) << ", std::log(vv+t4)=" << std::log(vv+t4));

      return ((t1*t2 - t3*std::log(t2+t1) -(vv*t4 - t3*std::log(vv+t4))) / (8*std::pow(uu, 1.5)));
    }


    double
    BezierCurve::arc_length(const size_t index, const double fraction) const
    {
      /*// This method uses a similair approach as https://malczak.info/blog/quadratic-bezier-curve-length
      // but instead of the full length, we integrate the full equaion (https://www.wolframalpha.com/input?i=integrate+sqrt%28A*t%5E2%2BB*t%2Bc%29+dt)
      // leaving t in there. Then we compute the integration constant by stating that the length at t=0 should
      // be zero.
      auto dt = points[index]-2.*control_points[index]+points[index+1];
      auto ddt = 2.*control_points[index]-2.*points[index];
      const double a = 4*(dt[0]*dt[0]+dt[1]*dt[1]);
      const double c = ddt[0]*ddt[0]+ddt[1]*ddt[1];

      if (a < 5e-4 * c || points[index] == control_points[index] || control_points[index] == points[index+1])
        {
          // all points are on a line
          return sqrt(((points[index+1]-points[index])*fraction)*((points[index+1]-points[index])*fraction));//std::sqrt((points[index][0] + dx1*t)*(points[index][0] + dx1*t) + (points[index][1] + dy1*t)*(points[index][1] + dy1*t));
        }

      const double b = 4*(dt[0]*ddt[0]+dt[1]*ddt[1]);

      const double inv_4a = 1./(4*a);
      const double u = (b*b)*inv_4a;
      const double k = (4*a*c-b*b)*inv_4a;
      const double sqrt_a = sqrt(a);
      const double sqrt_c = sqrt(c);
      const double inv_8_sqrt_a_a_a = 1./(8.*sqrt(a*a*a));
      const double sqrt_c_t_b_at = sqrt(c+fraction*(b+a*fraction));
      double x = fraction*sqrt_a+sqrt(u);

      // todo: optimize
      const double integral = ((b+2.*a*fraction)*sqrt_c_t_b_at)*inv_4a - ((b*b-4.*a*c)*log(b+2.*a*fraction+2.*sqrt_a*sqrt_c_t_b_at))*inv_8_sqrt_a_a_a;
      const double constant = (b*sqrt_c)*inv_4a - ((b*b-4.*a*c)*log(b+2.*sqrt_a*sqrt_c))*inv_8_sqrt_a_a_a;

      return integral-constant;*/
      return 0;
    }


    ClosestPointOnCurve
    BezierCurve::closest_point_on_curve_segment(const Point<2> &check_point) const
    {
      ClosestPointOnCurve closest_point_on_curve;
      const Point<2> &cp = check_point;
      double min_squared_distance = std::numeric_limits<double>::infinity();
      if (check_point.get_coordinate_system() == CoordinateSystem::cartesian)
        {
          for ( size_t cp_i = 0; cp_i < control_points.size(); ++cp_i)
            {
              std::stringstream output;
              const Point<2> &p1 = points[cp_i];
              const Point<2> &p2 = points[cp_i+1];
              min_squared_distance = std::min(std::min(min_squared_distance,(check_point-p1).norm_square()),(check_point-p1).norm_square());

              // Getting an estimate for where the closest point is with a linear approximation
              Point<2> P1P2 = p2-p1;
              Point<2> P1Pc = check_point-p1;

              double est =  std::min(1.,std::max(0.,(P1Pc*P1P2) / (P1P2*P1P2))); // est=estimate of solution
              bool found = false;

              Point<2> a = 3.*control_points[cp_i][0]-3.*control_points[cp_i][1]+points[cp_i+1]-points[cp_i];
              Point<2> b = 3.*points[cp_i] - 6.*control_points[cp_i][0]+3.*control_points[cp_i][1];
              Point<2> c = -3.*points[cp_i] + 3.*control_points[cp_i][0];
              Point<2> d = points[cp_i];

              Point<2> estimate_point = a*est*est*est+b*est*est+c*est+d;
              double min_squared_distance_cartesian_temp = estimate_point.cheap_relative_distance_cartesian(cp);

              for (size_t newton_i = 0; newton_i < 150; newton_i++)
                {
                  // based on https://stackoverflow.com/questions/2742610/closest-point-on-a-cubic-bezier-curve
                  a = 3.*control_points[cp_i][0]-3.*control_points[cp_i][1]+points[cp_i+1]-points[cp_i];
                  b = 3.*points[cp_i] - 6.*control_points[cp_i][0]+3.*control_points[cp_i][1];
                  c = -3.*points[cp_i] + 3.*control_points[cp_i][0];
                  d = points[cp_i];
#ifndef NDEBUG
                  output << "  wolfram alpha: (" << a[0] << "*x^3+" << b[0] << "*x^2+"<< c[0] << "*x+" << d[0] << "-" << cp[0] << ")^2+(" << a[1] << "*x^3+" << b[1] << "*x^2+"<< c[1] << "*x+" << d[1] << "-" << cp[1] << ")^2 with x=" << est << std::endl;
#endif
                  estimate_point = a*est*est*est+b*est*est+c*est+d;
                  const double squared_distance_cartesian = estimate_point.cheap_relative_distance_cartesian(cp);

                  const double squared_distance_cartesian_derivative = 2.0*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])
                                                                       + 2.0*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]);
                  const double squared_distance_cartesian_second_derivative  = 2.0*(6.0*a[0]*est + 2.0*b[0])*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])
                                                                               + 2.0*(3.0*a[0]*est*est + 2.0*b[0]*est + c[0])*(3.0*a[0]*est*est + 2.0*b[0]*est + c[0])
                                                                               + 2.0*(6.0*a[1]*est + 2.0*b[1])*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])
                                                                               + 2.0*(3.0*a[1]*est*est + 2.0*b[1]*est + c[1])*(3.0*a[1]*est*est + 2.0*b[1]*est + c[1]) ;

                  // the local minimum is where  squared_distance_cartesian_derivative=0 and squared_distance_cartesian_derivative>=0
                  const double update = std::min(0.5,std::max(-0.5,squared_distance_cartesian_derivative/std::fabs(squared_distance_cartesian_second_derivative)));//std::min(0.25,std::max(0.25,squared_distance_cartesian_derivative/std::fabs(squared_distance_cartesian_second_derivative)));
                  double line_search = 1.;
                  double est_test = est-update*line_search;
                  double squared_distance_cartesian_test = squared_distance_cartesian;
                  double squared_distance_cartesian_test_previous = squared_distance_cartesian;
                  double squared_distance_cartesian_derivative_test = squared_distance_cartesian_derivative;
                  double line_search_step = 2./3.;

                  for (unsigned int i = 0; i < 10; i++)
                    {
                      est_test = est-update*line_search;
                      estimate_point = a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d;

                      squared_distance_cartesian_test = (a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-cp[0])*(a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-cp[0])
                                                        +(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-cp[1])*(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-cp[1]);

#ifndef NDEBUG
                      squared_distance_cartesian_derivative_test = 2.0*(3.0*a[0]*est_test*est_test+2.0*b[0]*est_test+c[0])*(a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-cp[0])
                                                                   + 2.0*(3.0*a[1]*est_test*est_test+2.0*b[1]*est_test+c[1])*(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-cp[1]);
                      double squared_distance_cartesian_second_derivative_test = 2.0*(6.0*a[0]*est_test + 2.0*b[0])*(a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-cp[0])
                                                                                 + 2.0*(3.0*a[0]*est_test*est_test + 2.0*b[0]*est_test + c[0])*(3.0*a[0]*est_test*est_test + 2.0*b[0]*est_test + c[0])
                                                                                 + 2.0*(6.0*a[1]*est_test + 2.0*b[1])*(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-cp[1])
                                                                                 + 2.0*(3.0*a[1]*est_test*est_test + 2.0*b[1]*est_test + c[1])*(3.0*a[1]*est_test*est_test + 2.0*b[1]*est_test + c[1]) ;
                      output << "    i: " << cp_i << ", ni: " << newton_i<< ", lsi: " << i << ", line_search_step=" << line_search_step << ": squared_distance_cartesian_test = " << squared_distance_cartesian_test << ", diff= " << squared_distance_cartesian_test-squared_distance_cartesian << ", tests: " << (squared_distance_cartesian_test_previous < squared_distance_cartesian ? "true" : "false") << ":" << (squared_distance_cartesian_test > squared_distance_cartesian_test_previous ? "true" : "false") << ", est_test=" << est_test << ", update=" << update << ", ls=" << line_search << ", up*ls=" << update *line_search << ", test deriv =" << squared_distance_cartesian_derivative_test  << ", test upate=" << squared_distance_cartesian_derivative_test/fabs(squared_distance_cartesian_second_derivative_test) << ", p1=" << p1 << ", p2= " << p2 << ", poc= " << a *est_test *est_test *est_test+b *est_test *est_test+c *est_test+d << ", cp= " <<  check_point << ", ds:" << ((a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d)-check_point).norm_square() << ":" << min_squared_distance_cartesian_temp << ", diff = " << squared_distance_cartesian_test-min_squared_distance_cartesian_temp<< std::endl;
#endif
                      if (i > 0 && (squared_distance_cartesian_test > squared_distance_cartesian_test_previous))
                        {
                          if (squared_distance_cartesian_test_previous-squared_distance_cartesian < 0)
                            {
                              line_search *= 1/line_search_step;
                              break;
                            }
                          if (i> 1)
                            {
                              line_search *= (1/line_search_step)*(1/line_search_step);
                              est_test = est-update*line_search;
                              estimate_point = a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d;

                              squared_distance_cartesian_test_previous = (a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-cp[0])*(a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-cp[0])
                                                                         +(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-cp[1])*(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-cp[1]);
                              line_search_step = std::min(line_search_step*(11./10.),0.95);
                              continue;
                            }
                        }
                      squared_distance_cartesian_test_previous = squared_distance_cartesian_test;

                      line_search *= line_search_step;
                    }

                  est -= update*line_search;

                  min_squared_distance_cartesian_temp =  (a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])
                                                         +(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]);

                  if (std::fabs(update) < 1e-4)
                    {
                      found = true;
                      if (min_squared_distance_cartesian_temp < min_squared_distance)
                        {
                          if (est >= -1e-8 && est-1. <= 1e-8)
                            {
                              min_squared_distance = squared_distance_cartesian;
                              const Point<2> point_on_curve = a*est*est*est+b*est*est+c*est+d;

                              // the sign is rotating the derivative by 90 degrees.
                              // When moving in the direction of increasing t, left is positve and right is negative.
                              //const double &t = est;
                              //const Point<2> point_on_curve = a*est*est*est+b*est*est+c*est+d;
                              //(1-t)*(1-t)*(1-t)*points[cp_i] + 3*(1-t)*(1-t)*t *control_points[cp_i][0] + 3.*(1-t)*t *t *control_points[cp_i][1]+t *t *t *points[cp_i+1];
                              //const Point<2> derivative_point = ((-2+2*est)*points[cp_i] + (2-4*est)*control_points[cp_i] + 2*est*points[cp_i+1]);
                              //const Point<2> second_derivative_point = ((2)*points[cp_i] + (-4)*control_points[cp_i] + 2*points[cp_i+1]);
                              // https://www.wolframalpha.com/input?i=d%2Fdt+%281-t%29*%281-t%29*%281-t%29*a+%2B+3*%281-t%29*%281-t%29*t*b+%2B+3.*%281-t%29*t*t*c%2Bt*t*t*d
                              const Point<2> derivative_point = points[cp_i]*((6.-3.*est)*est-3.) + control_points[cp_i][0]*(est*(9*est-12)+3)
                                                                + control_points[cp_i][1]*(6.-9.*est)*est + points[cp_i+1]*3.*est*est;
                              //https://www.wolframalpha.com/input?i=d%2Fdt+d%2Fdt+%281-t%29*%281-t%29*%281-t%29*a+%2B+3*%281-t%29*%281-t%29*t*b+%2B+3.*%281-t%29*t*t*c%2Bt*t*t*d
                              const Point<2> second_derivative_point = points[cp_i]*(6-6*est) + control_points[cp_i][0]*(18*est-12)
                                                                       + control_points[cp_i][1]*-18*est+ points[cp_i+1]*6*est;

                              Point<2> tangent_point = derivative_point - point_on_curve;
                              // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
                              const double dot_product = (tangent_point*(check_point-point_on_curve));
                              const double sign = dot_product < 0. ? -1. : 1.;
                              tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

                              closest_point_on_curve.distance = sign*std::sqrt(min_squared_distance);
                              closest_point_on_curve.parametric_fraction = est;
                              closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
                              closest_point_on_curve.index = cp_i;
                              closest_point_on_curve.point = point_on_curve;
                              Point<2> normal = point_on_curve;
                              {
                                Point<2> derivative = a*est*est+b*est+c;
                                normal=derivative;
                                double normal_size = derivative.norm();
                                normal[0] = derivative[1]/normal_size;
                                normal[1] = -derivative[0]/normal_size;
                              }
                              closest_point_on_curve.normal = normal;
                            }
                        }
                      break;
                    }
                }
              //std::cout << output.str();
              WBAssertThrow(found, "Could not find a good solution. " << output.str());
            }
        }
      else
        {
          for ( size_t cp_i = 0; cp_i < control_points.size(); ++cp_i)
            {
              const Point<2> &p1 = points[cp_i];
              const Point<2> &p2 = points[cp_i+1];
              // Getting an estimate for where the closest point is with a linear approximation
              Point<2> P1P2 = p2-p1;
              Point<2> P1Pc = check_point-p1;

              double est =  std::min(1.,std::max(0.,(P1Pc*P1P2) / (P1P2*P1P2)));
              bool found = false;
              std::stringstream output;
              Point<2> a = 3.*control_points[cp_i][0]-3.*control_points[cp_i][1]+points[cp_i+1]-points[cp_i];
              Point<2> b = 3.*points[cp_i] - 6.*control_points[cp_i][0]+3.*control_points[cp_i][1];
              Point<2> c = -3.*points[cp_i] + 3.*control_points[cp_i][0];
              Point<2> d = points[cp_i];

              Point<2> estimate_point = a*est*est*est+b*est*est+c*est+d;
              const double cos_cp_lat = cos(cp[1]);
              double cos_lat = cos(estimate_point[1]);
              double sin_d_long_h = sin((estimate_point[0]-cp[0])*0.5);
              double sin_d_lat_h = sin((estimate_point[1]-cp[1])*0.5);
              double min_squared_distance_cartesian_temp = sin_d_lat_h*sin_d_lat_h+sin_d_long_h*sin_d_long_h*cos_cp_lat*cos_lat;
#ifndef NDEBUG
              output << "cp_i=" << cp_i << ", init est = " << est << ", min_squared_distance = " << min_squared_distance << ", min_squared_distance_cartesian_temp: " << min_squared_distance_cartesian_temp << ", p1: " << p1 << ", p2: " << p2 << std::endl;
              output  << std::setprecision(5) << "  wolfram: sin((" << a[1] << "*x^3+" << b[1] << "*x^2+"<< c[1] << "*x+" << d[1] << "-" << cp[1] << ")*.5)^2+sin((" << a[0] << "*x^3+" << b[0] << "*x^2+"<< c[0] << "*x+" << d[0] << "-" << cp[0] << ")*.5)^2*cos(" << cp[1] << ")*cos(" << a[1] << "*x^3+" << b[1] << "*x^2+"<< c[1] << "*x+" << d[1] << "-" << cp[1] << ") with x=" << est << std::endl;
              output  << std::setprecision(10) << "  python: y=np.sin((" << a[1] << "*x**3+" << b[1] << "*x**2+"<< c[1] << "*x+" << d[1] << "-" << cp[1] << ")*.5)**2+np.sin((" << a[0] << "*x**3+" << b[0] << "*x**2+"<< c[0] << "*x+" << d[0] << "-" << cp[0] << ")*.5)**2*np.cos(" << cp[1] << ")*np.cos(" << a[1] << "*x**3+" << b[1] << "*x**2+"<< c[1] << "*x+" << d[1] << "-" << cp[1] << "); x=" << est << std::endl;
#endif
              for (size_t newton_i = 0; newton_i < 150; newton_i++)
                {
                  // based on https://stackoverflow.com/questions/2742610/closest-point-on-a-cubic-bezier-curve
                  estimate_point = a*est*est*est+b*est*est+c*est+d;

                  cos_lat = cos(estimate_point[1]);
                  sin_d_long_h = sin((estimate_point[0]-cp[0])*0.5);
                  sin_d_lat_h = sin((estimate_point[1]-cp[1])*0.5);
                  const double cos_d_lat = cos(estimate_point[1]-cp[1]);
                  const double squared_distance_cartesian_full = sin((a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])*0.5)*sin((a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])*0.5)+sin((a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])*0.5)*sin((a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])*0.5)*cos(cp[1])*cos(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]);
                  const double squared_distance_cartesian = sin_d_lat_h*sin_d_lat_h+sin_d_long_h*sin_d_long_h*cos_cp_lat*cos_d_lat;

                  double sin_dlat = sin(estimate_point[1]-cp[1]);
                  double cos_dlong = cos(estimate_point[0]-cp[0]);
                  double cos_dlong_h = cos(0.5*(estimate_point[0]-cp[0]));
                  double cos_dlat_h = cos(0.5*(estimate_point[1]-cp[1]));
                  double deriv_long = (3.0*a[0]*est*est+2.0*b[0]*est+c[0]);
                  double deriv_lat = (3.0*a[1]*est*est+2.0*b[1]*est+c[1]);

                  const double squared_distance_cartesian_derivative = cos_cp_lat*(-deriv_lat)*sin_d_long_h*sin_d_long_h*sin_dlat+cos_cp_lat*deriv_long*sin_d_long_h*cos_dlong_h*cos_d_lat+deriv_lat*sin_d_lat_h*cos_dlat_h;//cos_cp_lat*(-(3.0*a[1]*est*est+2.0*b[1]*est+c[1]))*sin((estimate_point[0]-cp[0])*0.5)*sin((estimate_point[0]-cp[0])*0.5)*sin_dlat+cos_cp_lat*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*sin((estimate_point[0]-cp[0])*0.5)*cos(0.5*(estimate_point[0]-cp[0]))*cos_lat+(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*sin(0.5*(estimate_point[1]-cp[1]))*cos(0.5*(estimate_point[1]-cp[1]));;//cos_cp_lat*(-deriv_lat)*sin_d_long_h*sin_d_long_h*sin_dlat+cos_cp_lat*deriv_long*sin_d_long_h*cos_dlong_h*cos_d_lat+deriv_lat*sin_d_lat_h*cos_dlat_h;
                  const double squared_distance_cartesian_second_derivative = cos_cp_lat*cos_d_lat*(-0.5*deriv_long*deriv_long*sin_d_long_h*sin_d_long_h+0.5*deriv_long*deriv_long*cos_dlong_h*cos_dlong_h+(6.0*a[0]*est+2.0*b[0])*sin_d_long_h*cos_dlong_h)+cos_cp_lat*sin_d_long_h*sin_d_long_h*(deriv_lat*deriv_lat*(-cos_d_lat)-(6.0*a[1]*est+2.0*b[1])*sin_dlat)-2.0*cos_cp_lat*deriv_long*deriv_lat*sin_d_long_h*cos_dlong_h*sin_dlat-0.5*deriv_lat*deriv_lat*sin_d_lat_h*sin_d_lat_h+0.5*deriv_lat*deriv_lat*cos_dlat_h*cos_dlat_h+(6.0*a[1]*est+2.0*b[1])*sin_d_lat_h*cos_dlat_h;

#ifndef NDEBUG
                  const double squared_distance_cartesian_derivative_full = cos(cp[1])*(-(3.0*a[1]*est*est+2.0*b[1]*est+c[1]))*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])+cos(cp[1])*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])+(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*sin(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))*cos(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]));
                  const double squared_distance_cartesian_second_derivative_full = cos(cp[1])*cos(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])*(-0.5*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))+0.5*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))+(6.0*a[0]*est+2.0*b[0])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])))+cos(cp[1])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*((3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(-cos(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))-(6.0*a[1]*est+2.0*b[1])*sin(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))-2.0*cos(cp[1])*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])-0.5*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*sin(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))*sin(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))+0.5*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*cos(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))*cos(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))+(6.0*a[1]*est+2.0*b[1])*sin(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))*cos(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]));
                  output <<"sqd = " << squared_distance_cartesian <<":" << squared_distance_cartesian_full << ", diff=" << squared_distance_cartesian-squared_distance_cartesian_full << ", sqdd: " << squared_distance_cartesian_derivative <<":" << squared_distance_cartesian_derivative_full << ", diff="<< squared_distance_cartesian_derivative-squared_distance_cartesian_derivative_full << ", sqdd: " << squared_distance_cartesian_second_derivative << ":" << squared_distance_cartesian_second_derivative_full << ", diff= " << squared_distance_cartesian_second_derivative-squared_distance_cartesian_second_derivative_full << ", est: " << est << std::endl;
#endif
                  // the local minimum is where  squared_distance_cartesian_derivative=0 and squared_distance_cartesian_derivative>=0
                  const double update = std::min(0.5,std::max(-0.5,squared_distance_cartesian_derivative/std::fabs(squared_distance_cartesian_second_derivative)));//std::min(0.25,std::max(0.25,squared_distance_cartesian_derivative/std::fabs(squared_distance_cartesian_second_derivative)));
                  double line_search = 1.;
                  double est_test = est-update*line_search;
                  double squared_distance_cartesian_test = squared_distance_cartesian;
                  double squared_distance_cartesian_test_previous = squared_distance_cartesian;
                  double squared_distance_cartesian_derivative_test = squared_distance_cartesian_derivative;
                  double line_search_step = 2./3.;

                  for (unsigned int i = 0; i < 10; i++)
                    {
                      est_test = est-update*line_search;
                      estimate_point = a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d;

                      cos_lat = cos(estimate_point[1]);
                      sin_d_long_h = sin((estimate_point[0]-cp[0])*0.5);
                      sin_d_lat_h = sin((estimate_point[1]-cp[1])*0.5);
                      squared_distance_cartesian_test = sin_d_lat_h*sin_d_lat_h+sin_d_long_h*sin_d_long_h*cos_cp_lat*cos(estimate_point[1]-cp[1]);

#ifndef NDEBUG
                      sin_dlat = sin(estimate_point[1]-cp[1]);
                      cos_dlong = cos(estimate_point[0]-cp[0]);
                      deriv_long = (3.0*a[0]*est_test*est_test+2.0*b[0]*est_test+c[0]);
                      deriv_lat = (3.0*a[1]*est_test*est_test+2.0*b[1]*est_test+c[1]);
                      squared_distance_cartesian_derivative_test = cos_cp_lat*(-deriv_lat)*sin_d_long_h*sin_d_long_h*sin_dlat+cos_cp_lat*deriv_long*sin_d_long_h*cos_dlong_h*cos_d_lat+deriv_lat*sin_d_lat_h*cos_dlat_h;
                      double squared_distance_cartesian_second_derivative_test = cos_cp_lat*cos_d_lat*(-0.5*deriv_long*deriv_long*sin_d_long_h*sin_d_long_h+0.5*deriv_long*deriv_long*cos_dlong_h*cos_dlong_h+(6.0*a[0]*est_test+2.0*b[0])*sin_d_long_h*cos_dlong_h)+cos_cp_lat*sin_d_long_h*sin_d_long_h*(deriv_lat*deriv_lat*(-cos_d_lat)-(6.0*a[1]*est_test+2.0*b[1])*sin_dlat)-2.0*cos_cp_lat*deriv_long*deriv_lat*sin_d_long_h*cos_dlong_h*sin_dlat-0.5*deriv_lat*deriv_lat*sin_d_lat_h*sin_d_lat_h+0.5*deriv_lat*deriv_lat*cos_dlat_h*cos_dlat_h+(6.0*a[1]*est_test+2.0*b[1])*sin_d_lat_h*cos_dlat_h;
                      output << "    i: " << cp_i << ", ni: " << newton_i<< ", lsi: " << i << ", line_search_step=" << line_search_step << ": squared_distance_cartesian_test = " << squared_distance_cartesian_test << ", diff= " << squared_distance_cartesian_test-squared_distance_cartesian << ", tests: " << (squared_distance_cartesian_test_previous < squared_distance_cartesian ? "true" : "false") << ":" << (squared_distance_cartesian_test > squared_distance_cartesian_test_previous ? "true" : "false") << ", est_test=" << est_test << ", update=" << update << ", ls=" << line_search << ", up*ls=" << update *line_search << ", test deriv =" << squared_distance_cartesian_derivative_test  << ", test upate=" << squared_distance_cartesian_derivative_test/fabs(squared_distance_cartesian_second_derivative_test) << ", p1=" << p1 << ", p2= " << p2 << ", poc= " << a *est_test *est_test *est_test+b *est_test *est_test+c *est_test+d << ", cp= " <<  check_point << ", ds:" << ((a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d)-check_point).norm_square() << ":" << min_squared_distance_cartesian_temp << ", diff = " << squared_distance_cartesian_test-min_squared_distance_cartesian_temp<< std::endl;
#endif
                      if (i > 0 && (squared_distance_cartesian_test > squared_distance_cartesian_test_previous))
                        {
                          if (squared_distance_cartesian_test_previous-squared_distance_cartesian < 0)
                            {
                              line_search *= 1/line_search_step;
                              break;
                            }
                          if (i> 1)
                            {
                              line_search *= (1/line_search_step)*(1/line_search_step);//6./5.;
                              est_test = est-update*line_search;
                              estimate_point = a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d;

                              cos_lat = cos(estimate_point[1]);
                              sin_d_long_h = sin((estimate_point[0]-cp[0])*0.5);
                              sin_d_lat_h = sin((estimate_point[1]-cp[1])*0.5);
                              squared_distance_cartesian_test_previous = sin_d_lat_h*sin_d_lat_h+sin_d_long_h*sin_d_long_h*cos_cp_lat*cos(estimate_point[1]-cp[1]);
                              line_search_step = std::min(line_search_step*(11./10.),0.95);
                              continue;
                            }
                        }
                      squared_distance_cartesian_test_previous = squared_distance_cartesian_test;

                      line_search *= line_search_step;
                    }
#ifndef NDEBUG
                  output << "    i: " << cp_i << ", ni: " << newton_i<< ", est= " << est-update *line_search << ", ls:" << line_search << ": squared_distance_cartesian_test = " << squared_distance_cartesian_test << ", diff= " << squared_distance_cartesian_test-squared_distance_cartesian << std::endl;
#endif
                  est -= update*line_search;

                  estimate_point = a*est*est*est+b*est*est+c*est+d;

                  cos_lat = cos(estimate_point[1]);
                  sin_d_long_h = sin((estimate_point[0]-cp[0])*0.5);
                  sin_d_lat_h = sin((estimate_point[1]-cp[1])*0.5);

                  min_squared_distance_cartesian_temp = sin_d_lat_h*sin_d_lat_h+sin_d_long_h*sin_d_long_h*cos_cp_lat*cos(estimate_point[1]-cp[1]);

                  if (std::fabs(update) < 1e-4)
                    {
                      found = true;
                      if (min_squared_distance_cartesian_temp < min_squared_distance)
                        {
                          if (est >= -1e-8 && est-1. <= 1e-8)
                            {
                              min_squared_distance = squared_distance_cartesian;
                              const Point<2> point_on_curve = a*est*est*est+b*est*est+c*est+d;

                              // the sign is rotating the derivative by 90 degrees.
                              // When moving in the direction of increasing t, left is positve and right is negative.
                              //const double &t = est;
                              //const Point<2> point_on_curve = a*est*est*est+b*est*est+c*est+d;
                              //(1-t)*(1-t)*(1-t)*points[cp_i] + 3*(1-t)*(1-t)*t *control_points[cp_i][0] + 3.*(1-t)*t *t *control_points[cp_i][1]+t *t *t *points[cp_i+1];
                              //const Point<2> derivative_point = ((-2+2*est)*points[cp_i] + (2-4*est)*control_points[cp_i] + 2*est*points[cp_i+1]);
                              //const Point<2> second_derivative_point = ((2)*points[cp_i] + (-4)*control_points[cp_i] + 2*points[cp_i+1]);
                              // https://www.wolframalpha.com/input?i=d%2Fdt+%281-t%29*%281-t%29*%281-t%29*a+%2B+3*%281-t%29*%281-t%29*t*b+%2B+3.*%281-t%29*t*t*c%2Bt*t*t*d
                              const Point<2> derivative_point = points[cp_i]*((6.-3.*est)*est-3.) + control_points[cp_i][0]*(est*(9*est-12)+3)
                                                                + control_points[cp_i][1]*(6.-9.*est)*est + points[cp_i+1]*3.*est*est;
                              //https://www.wolframalpha.com/input?i=d%2Fdt+d%2Fdt+%281-t%29*%281-t%29*%281-t%29*a+%2B+3*%281-t%29*%281-t%29*t*b+%2B+3.*%281-t%29*t*t*c%2Bt*t*t*d
                              const Point<2> second_derivative_point = points[cp_i]*(6-6*est) + control_points[cp_i][0]*(18*est-12)
                                                                       + control_points[cp_i][1]*-18*est+ points[cp_i+1]*6*est;

                              Point<2> tangent_point = derivative_point - point_on_curve;
                              // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
                              const double dot_product = (tangent_point*(check_point-point_on_curve));
                              const double sign = dot_product < 0. ? -1. : 1.;
                              tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

                              closest_point_on_curve.distance = sign*std::sqrt(min_squared_distance);
                              closest_point_on_curve.parametric_fraction = est;
                              closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
                              closest_point_on_curve.index = cp_i;
                              closest_point_on_curve.point = point_on_curve;
                              Point<2> normal = point_on_curve;
                              {
                                Point<2> derivative = a*est*est+b*est+c;
                                normal=derivative;
                                double normal_size = derivative.norm();
                                normal[0] = derivative[1]/normal_size;
                                normal[1] = -derivative[0]/normal_size;
                              }
                              closest_point_on_curve.normal = normal;
                            }
                        }
                      break;
                    }
                }
              WBAssertThrow(found, "Could not find a good solution. " << output.str());
            }

        }
      return closest_point_on_curve;
    }


    ClosestPointOnCurve
    BezierCurve::closest_point_on_curve(const Point<2> &check_point) const
    {
      ClosestPointOnCurve closest_point_on_curve = this->closest_point_on_curve_segment(check_point);
      if (std::isnan(closest_point_on_curve.point[0]))
        {
          // The closest point is one of the edges, find which one
          const double distance_to_first_point = (points[0]-check_point).norm();
          const double distance_to_last_point = (points[points.size()-1]-check_point).norm();
          if (distance_to_first_point < distance_to_last_point)
            {
              closest_point_on_curve.distance = NaN::DSNAN; // the distance computed above is not following the curve. Could be computed if needed.
              closest_point_on_curve.index = 0;
              closest_point_on_curve.parametric_fraction = 0;
              closest_point_on_curve.interpolation_fraction = 0;
              closest_point_on_curve.point = points[0];
            }
          else
            {
              closest_point_on_curve.distance = NaN::DSNAN; // the distance computed above is not following the curve. Could be computed if needed.
              closest_point_on_curve.index = points.size()-1;
              closest_point_on_curve.parametric_fraction = 0;
              closest_point_on_curve.interpolation_fraction = 0;
              closest_point_on_curve.point = points[points.size()-1];
            }
        }
      return closest_point_on_curve;
    }


    double *
    BezierCurve::solve_cubic_equation_real(const double a_original,const double b_original,const double c_original,const double d_original, double real_roots[4])
    {
      // todo: look at https://stackoverflow.com/questions/2003465/fastest-numerical-solution-of-a-real-cubic-polynomial
      //real_roots = {NaN::DSNAN,NaN::DSNAN,NaN::DSNAN,0.};
      //real_roots.reserve(3);
      //constexpr double one_third = 1./3.;
      size_t index = real_roots[3];
      constexpr double tolerance = 1e-10;
      if (std::abs(a_original) <= tolerance)
        {
          if (std::abs(b_original) <= tolerance)
            {
              // linear equation
              const double &a = c_original;
              const double &b = d_original;
              real_roots[index] = -b/a;
              index++;
              real_roots[3] = index;
              return real_roots;
            }
          // quadratic equation
          const double &a = b_original;
          const double &b = c_original;
          const double &c = d_original;
          const double discriminant = b*b -4.*a*c;
          if (std::abs(discriminant) <= tolerance)
            {
              real_roots[index] = (-b+sqrt(discriminant))/(2*a);
              index++;
              real_roots[3] = index;
              return real_roots;
            }
          else if (discriminant > 0)
            {
              real_roots[index] = (-b + sqrt(discriminant))/(2.*a);
              real_roots[index+1] = (-b - sqrt(discriminant))/(2.*a);
              index += 2;
              real_roots[3] = index;
              return real_roots;
            }
          return real_roots;
        }
      else
        {
          // based on https://quarticequations.com/Cubic.pdf
          // divide by a
          //constexpr double one_third = 1./3.;
          const double b = b_original/a_original; // a_2
          const double c = c_original/a_original; // a_1
          const double d = d_original/a_original; //a_0

          const double q = (c/3.)-(b*b/9.);
          const double r = (c*b-3*d)/6. - b*b*b/27.;
          const double discriminant = r*r+q*q*q;

          if (discriminant > 0)
            {
              // only one real solution
              const double A = std::pow(std::abs(r) + sqrt(discriminant),1./3.);
              const double t = r >= 0 ? A-q/A : q/A-A;
              real_roots[0] = t-b/3.;
              index++;
            }
          else
            {
              // three real solutions
              // std::pow(-q,3./2.) == sqrt(-q*q*q)
              const double phi_1 = std::abs(q) <= tolerance ? 0 : (acos(r/sqrt(-q*q*q))/3.);
              const double phi_2 = phi_1 - 2.*Consts::PI/3.;
              const double phi_3 = phi_1 + 2.*Consts::PI/3.;
              const double sqrt_q_3 = 2*sqrt(-q);
              const double b_t_one_third = b/3.;
              const double value_1 = sqrt_q_3 * FT::cos(phi_1)-b_t_one_third;
              const double value_2 = sqrt_q_3 * FT::cos(phi_2)-b_t_one_third;
              const double value_3 = sqrt_q_3 * FT::cos(phi_3)-b_t_one_third;

              real_roots[0] = value_1;
              index++;
              if (std::abs(value_1 - value_2) > tolerance)
                {
                  real_roots[index] = value_2;
                  index++;
                }
              // we don't have to check value 1 and 3 because z3 <= z2 <= z1
              // so if z2 and z3 are not equal, z1 and z3 are not either.
              if (std::abs(value_2 - value_3) > tolerance)
                {
                  real_roots[index] = value_3;
                  index++;
                }
            }
        }
      real_roots[3] = index;
      return real_roots;
    }

  }
}