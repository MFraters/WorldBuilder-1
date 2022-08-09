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

using namespace WorldBuilder;

namespace WorldBuilder
{
  namespace Objects
  {

    double
    BezierCurve::arc_length(const Point<2> &a, const Point<2>  &b, const Point<2> &c) const
    {

      // compute the arc length of the Bezier curve
      // see https://gamedev.stackexchange.com/questions/6009/bezier-curve-arc-length
      Point<2> v = a;
      Point<2> w = b;
      v[0] = 2*(b[0] - a[0]);
      v[1] = 2*(b[1] - a[1]);
      w[0] = c[0] - 2*b[0] + a[0];
      w[1] = c[1] - 2*b[1] + a[1];

      const double uu = 4*(w[0]*w[0] + w[1]*w[1]);

      if (uu < 0.00001)
        {
          return std::sqrt((c[0] - a[0])*(c[0] - a[0]) + (c[1] - a[1])*(c[1] - a[1]));
        }

      const double vv = 4*(v[0]*w[0] + v[1]*w[1]);
      const double ww = v[0]*v[0] + v[1]*v[1];

      const double t1 = (2*std::sqrt(uu*(uu + vv + ww)));
      const double t2 = 2*uu+vv;
      const double t3 = vv*vv - 4*uu*ww;
      const double t4 = (2*std::sqrt(uu*ww));

      return ((t1*t2 - t3*std::log(t2+t1) -(vv*t4 - t3*std::log(vv+t4))) / (8*std::pow(uu, 1.5)));
    }

    /*double blen()
    {
      const double p0x = 0.;
      const double p0y = 4.;
      const double p1x = 0.787577;
      const double p1y = 3.49613;
      const double p2x = 1.;
      const double p2y = 3.;
      const double ax = p0x - 2*p1x + p2x;
      const double ay = p0y - 2*p1y + p2y;
      const double bx = 2*p1x - 2*p0x;
      const double by = 2*p1y - 2*p0y;
      double A = 4*(ax*ax + ay*ay);
      double B = 4*(ax*bx + ay*by);
      double C = bx*bx + by*by;

      std::cout << "a: " << ax << " " << ay << ", b: " << bx << " " << by
                << ", A: " << A << ", B: " << B << ", C: " << C << ", ";

      double Sabc = 2*sqrt(A+B+C);
      double A_2 = sqrt(A);
      double A_32 = 2*A*A_2;
      double C_2 = 2*sqrt(C);
      double BA = B/A_2;

      return ( A_32*Sabc +
               A_2*B*(Sabc-C_2) +
               (4*C*A-B*B)*log( (2*A_2+BA+Sabc)/(BA+C_2) )
             )/(4*A_32);
    };*/

    double
    BezierCurve::arc_length(const Point<2> &P0, const Point<2>  &Pc, const Point<2> &P1, const double t) const
    {

      //std::cout << "blen = " << blen() << std::endl;

      // This method uses a similair approach as https://malczak.info/blog/quadratic-bezier-curve-length
      // but instead of the full length, we integrate the full equaion (https://www.wolframalpha.com/input?i=integrate+sqrt%28A*t%5E2%2BB*t%2Bc%29+dt)
      // leaving t in there. Then we compute the integration constant by stating that the length at t=0 should
      // be zero.
      auto dt = P0-2.*Pc+P1;
      auto ddt = 2.*Pc-2.*P0;
      //std::cout << "P0: " << P0 << ", Pc: " << Pc << ", P1: " << P1 << ", dt: " << dt << ", ddt: " << ddt << std::endl;
      const double a = 4*(dt[0]*dt[0]+dt[1]*dt[1]);//dt*dt;//
      const double c = ddt[0]*ddt[0]+ddt[1]*ddt[1];//ddt*ddt;//

      if (a < 5e-4 * c)//(c < 0.0000001)
        {
          // all points are on a line
          // todo: chec that this is actually correct
          WBAssert(false,"todo");
          std::cout << "everyting in a line!" << std::endl;
          return -1;//std::sqrt((P0[0] + dx1*t)*(P0[0] + dx1*t) + (P0[1] + dy1*t)*(P0[1] + dy1*t));
        }

      const double b = 4*(dt[0]*ddt[0]+dt[1]*ddt[1]);//2*ddt*dt;////at*bt;

      const double u = (b*b)/(4*a);
      const double k = (4*a*c-b*b)/(4*a);
      double x = t*sqrt(a)+sqrt(u);

      //std::cout << " with a=" << a << ", b=" << b << ", c=" << c << ", k = " << k << ", u = " << u << ", x = " << x << std::endl;

      const double integral = ((b+2.*a*t)*sqrt(c+t*(b+a*t)))/(4*a) - ((b*b-4.*a*c)*log(b+2.*a*t+2.*sqrt(a)*sqrt(c+t*(b+a*t))))/(8.*pow(a,(3./2.)));
      //const double integral = ((2.*a*t+b)*sqrt(t*(a*t+b)+c))/(4*a) - ((b*b-4.*a*c)*log(2.*sqrt(a)*sqrt(t*(a*t+b)+c)+2.*a*t+b))/(8*pow(a,1.5));
      //const double integral = 0.5 * (x * sqrt(k + x*x) + k * log(x + sqrt(k + x*x)));//0.5*(x*sqrt(x*x+k)+k*log(x+sqrt(x*x+k)));
      //// set the constant so that t=0 returns length zero
//
      //const double x0 = 0.*sqrt(a)+sqrt(u);
      //const double integral0 = 0.5*(x0*sqrt(x0*x0+k)+k*log(x0+sqrt(x0*x0+k)));
      //const double constant = k*log(sqrt(k));
      //x = sqrt(u);
      //const double constant =  0.5 * (x * sqrt(k + x*x) + k * log(x + sqrt(k + x*x)));//((b+2.*a*z)*sqrt(c+z*(b+a*z)))/(4*a) - ((b*b-4.*a*c)*log(b + 2.*a*z + 2.*sqrt(a)*sqrt(c + z*(b + a*z))))/(8.*pow(a,(3./2.)));
      const double z = 0.;
      const double constant = ((b+2.*a*z)*sqrt(c+z*(b+a*z)))/(4*a) - ((b*b-4.*a*c)*log(b+2.*a*z+2.*sqrt(a)*sqrt(c+z*(b+a*z))))/(8.*pow(a,(3./2.)));
      //const double constant = (b*sqrt(c))/(4*a) - (b*b-4*a*c*log(2*sqrt(a)*sqrt(c)+b))/(8*pow(a,1.5));

      //std::cout << "integral = " << integral << ", constant = " << constant << ", min = " << integral-constant << std::endl;

      return integral-constant;


      /*
            const double sabc = sqrt(a+b+c);
            const double a2 = pow(a,-0.5);
            const double a32 = a2*a2*a2;
            const double c2 = 2.*sqrt(c);
            const double ba_c2 = b* a2 + c2;

            const double v0 = 0.25*a2*a2*b*(2.*sabc - c2) + sabc;

            return v0 + 0.25* a32 * (4.0 * c * a - b*b)
                    *std::log(((2.0*a+b)*a2 + 2.0*sabc)/ba_c2);*/
      /*
      // todo: catch case all points are on a line
      const double dx1 = Pc[0]-P0[0];
      const double dx2 = P1[0]-Pc[0];
      const double dy1 = Pc[1]-P0[1];
      const double dy2 = P1[1]-Pc[1];

      const double a = (dx2-dx1)*(dx2-dx1)+(dy2-dy1)*(dy2-dy1);

      if (a < 0.0000001)
        {
          // all points are on a line
          // todo: chec that this is actually correct
          std::cout << "everyting in a line!" << std::endl;
          return std::sqrt((P0[0] + dx1*t)*(P0[0] + dx1*t) + (P0[1] + dy1*t)*(P0[1] + dy1*t));
        }
      const double b = dx1*(dx2-dx1)+(dy1*(dy2-dy1));
      const double c = dx1*dx1+dy1*dy1;
      std::cout << std::endl << "P0[0]: "  << P0[0] << ", Pc[0]: " << Pc[0]
                << ", dx1: " << dx1 << ", dy1: " << dy1
                << ", dx2: " << dx2 << ", dy2: " << dy2
                << ", c: " << c << ", b: " << b << ", a: " << a << ", b/a: " << b/a
                << ", part 1:" << (t+(b/a))*sqrt(c+2*b*t+a*t*t) << ", part 2: " << ((a*c-b*b)/std::pow(a,1.5))*asinh((a*t+b)/sqrt(a*c-b*b))
                << ", r:";

      return (t+(b/a))*sqrt(c+2*b*t+a*t*t)+((a*c-b*b)/std::pow(a,1.5))*asinh((a*t+b)/sqrt(a*c-b*b));
      */
    }



    BezierCurve::BezierCurve(const std::vector<Point<2> > &p, std::vector<double> &angle_constrains)
    {
      points = p;
      // first compute the factors for a monotome spline
      const size_t n_points = p.size();
      std::vector<double> x(n_points);
      std::vector<double> y(n_points);
      for (size_t i = 0; i < n_points; ++i)
        {
          x[i] = p[i][0];
          y[i] = p[i][1];
          std::cout << i << ": (" << x[i] << ":" << y[i] << ")" << std::endl;
        }
      CubicSpline x_spline(x);
      std::cout << "flag 1.1" << std::endl;
      CubicSpline y_spline(y);
      std::cout << "flag 1.2" << std::endl;
      //x_spline.set_points(x);
      //y_spline.set_points(y);

      // resize all vectors
      control_points.resize(n_points-1, p[0]);
      lengths.resize(n_points-1,std::numeric_limits<double>::signaling_NaN());
      angles.resize(n_points,std::numeric_limits<double>::signaling_NaN());
      if (std::isnan(angle_constrains[0]))
        {
          const double value_x = x_spline.m[0][3];
          const double value_y = y_spline.m[0][3];
          const double derivative_x = x_spline.m[0][2];
          const double derivative_y = y_spline.m[0][2];
          //const double value_x_q = x_spline.m[0][0]*0.25*0.25*0.25 + x_spline.m[0][1]*0.25*0.25 + x_spline.m[0][2]*0.25 + x_spline.m[0][3];
          //const double value_y_q = y_spline.m[0][0]*0.25*0.25*0.25 + y_spline.m[0][1]*0.25*0.25 + y_spline.m[0][2]*0.25 + y_spline.m[0][3];
          const double value_x_h = x_spline.m[0][0]*0.5*0.5*0.5 + x_spline.m[0][1]*0.5*0.5 + x_spline.m[0][2]*0.5 + x_spline.m[0][3];
          const double value_y_h = y_spline.m[0][0]*0.5*0.5*0.5 + y_spline.m[0][1]*0.5*0.5 + y_spline.m[0][2]*0.5 + y_spline.m[0][3];
          //const double value_x_3q = x_spline.m[0][0]*0.75*0.75*0.75 + x_spline.m[0][1]*0.75*0.75 + x_spline.m[0][2]*0.75 + x_spline.m[0][3];
          //const double value_y_3q = y_spline.m[0][0]*0.75*0.75*0.75 + y_spline.m[0][1]*0.75*0.75 + y_spline.m[0][2]*0.75 + y_spline.m[0][3];
          angles[0] = std::atan2(value_y_h-value_y,value_x_h-value_x);//std::atan2(derivative_y,derivative_x);
          std::cout << "0" << ": value_x:y = " << value_x << ":" << value_y << ",  value_h x:y = " << value_x_h << ":" << value_y_h << ", angle = " << angles[0] << " (" << angles[0] *180./M_PI << ")" << ", derivative_x = " << derivative_x
                    << ", x a:b:c:d = " << x_spline.m[0][0] << ":" << x_spline.m[0][1] << ":" << x_spline.m[0][2]  << ":" << x_spline.m[0][1]
                    << ", derivative_y = " << derivative_y
                    << ", y a:b:c:d = " << y_spline.m[0][0] << ":" << y_spline.m[0][1] << ":" << y_spline.m[0][2]  << ":" << y_spline.m[0][1]
                    << std::endl;
          // compute p1 and p2 relative to p0 to compute the angles.
          //angles[0] = std::atan2(p[1][1]-p[0][1],p[1][0]-p[0][0]);
          //const double angle_p0p1 = std::atan2(p[1][1]-p[0][1],p[1][0]-p[0][0]);
          //const double angle_p0p2 = std::atan2(p[2][1]-p[0][1],p[2][0]-p[0][0]);
          //angles[0] = angle_p0p1-0.5*(angle_p0p2-angle_p0p1); //2.*angle_p0p1-angle_p0p2;//
          //std::cout << "angles[0] = " << angles[0] << " (" << angles[0] *180./M_PI << ")" << ", angle_p0p1 = " << angle_p0p1  << "(" << angle_p0p1 *180./M_PI << ")"<< ", angle_p0p2 = " << angle_p0p2  << "(" << angle_p0p2 *180./M_PI << ")"<< std::endl;
        }
      else
        {
          // NOTE: start angle assumes slabs or faults going down, which means they should provide a negative angle to get expected behavoir
          angles[0] = angle_constrains[0];
        }

      for (size_t i = 0; i < n_points-2; ++i)
        {
          // first compute next angle:
          if (std::isnan(angle_constrains[i+1]))
            {
              const double value_x = x_spline.m[i+1][3];
              const double value_y = y_spline.m[i+1][3];
              const double derivative_x = x_spline.m[i+1][2];
              const double derivative_y = y_spline.m[i+1][2];
              const double value_x_h = x_spline.m[i+1][0]*0.5*0.5*0.5 + x_spline.m[i+1][1]*0.5*0.5 + x_spline.m[i+1][2]*0.5 + x_spline.m[i+1][3];
              const double value_y_h = y_spline.m[i+1][0]*0.5*0.5*0.5 + y_spline.m[i+1][1]*0.5*0.5 + y_spline.m[i+1][2]*0.5 + y_spline.m[i+1][3];
              //const double value_x = x_spline.m[i+1][0]*(i+1)*(i+1)*(i+1) + x_spline.m[i+1][1]*(i+1)*(i+1) + x_spline.m[i+1][2]*(i+1) + x_spline.m[i+1][3];
              //const double value_y = y_spline.m[i+1][0]*(i+1)*(i+1)*(i+1) + y_spline.m[i+1][1]*(i+1)*(i+1) + y_spline.m[i+1][2]*(i+1) + y_spline.m[i+1][3];
              //const double derivative_x = 3*x_spline.m[i+1][0]*(i+1)*(i+1) + 2*x_spline.m[i+1][1]*(i+1) + x_spline.m[i+1][2];
              //const double derivative_y = 3*y_spline.m[i+1][0]*(i+1)*(i+1) + 2*y_spline.m[i+1][1]*(i+1) + y_spline.m[i+1][2];
              angles[i+1] = std::atan2(value_y_h-value_y,value_x_h-value_x);
              std::cout << i+1 << ": value_x:y = " << value_x << ":" << value_y << ",  value_h x:y = " << value_x_h << ":" << value_y_h << ", angle = " << angles[i+1] << " (" << angles[i+1] *180./M_PI << ")" << ", derivative_x = " << derivative_x
                        << ", x a:b:c:d = " << x_spline.m[i+1][0] << ":" << x_spline.m[i+1][1] << ":" << x_spline.m[i+1][2]  << ":" << x_spline.m[i+1][1]
                        << ", derivative_y = " << derivative_y
                        << ", y a:b:c:d = " << y_spline.m[i+1][0] << ":" << y_spline.m[i+1][1] << ":" << y_spline.m[i+1][2]  << ":" << y_spline.m[i+1][1]
                        << std::endl;
              //const double &x = p[0][0];
              //const double derivative = 3*m[0][0]*x*x + 2*m[0][1]*x + m[0][2];
              //angles[i] = std::atan2(derivative,1.);
              // compute p1 and p2 relative to p0 to compute the angles.
              //if ( i < n_points-2)
              //  {
              //    //angles[i+1] = std::atan2(p[i+2][1]-p[i][1],p[i+2][0]-p[i][0]);
              //
              //    const double angle_p0p1 = std::atan2(p[i+2][1]-p[i+1][1],p[i+2][0]-p[i+1][0]);
              //    const double angle_p0p2 = std::atan2(p[i+3][1]-p[i+1][1],p[i+3][0]-p[i+1][0]);
              //    angles[i+1] = angle_p0p1-0.5*(angle_p0p2-angle_p0p1);//2.*angle_p0p1-angle_p0p2;
              //    std::cout << "angles[i+1] = " << angles[i+1] << " (" << angles[i+1] *180./M_PI << ")" << ", angle_p0p1 = " << angle_p0p1  << "(" << angle_p0p1 *180./M_PI << ")"<< ", angle_p0p2 = " << angle_p0p2  << "(" << angle_p0p2 *180./M_PI << ")"<< std::endl;
              //  }
              //else
              //  {
              //    //angles[i+1] = std::atan2(p[i+1][1]-p[i][1],p[i+1][0]-p[i][0]);
              //    const double angle_p0p1 = std::atan2(p[i][1]-p[i+1][1],p[i][0]-p[i+1][0]);
              //    const double angle_p0p2 = std::atan2(p[i-1][1]-p[i+1][1],p[i-1][0]-p[i+1][0]);
              //    angles[i+1] = 2.*angle_p0p1-angle_p0p2;
              //    std::cout << "angles[i+1] = " << angles[i+1] << " (" << angles[i+1] *180./M_PI << ")" << ", angle_p0p1 = " << angle_p0p1  << "(" << angle_p0p1 *180./M_PI << ")"<< ", angle_p0p2 = " << angle_p0p2  << "(" << angle_p0p2 *180./M_PI << ")"<< std::endl;
              //  }
            }
          else
            {
              // NOTE: start angle doesn't assumes slabs or faults going down, which means they should provide a negative angle to get expected behavoir
              angles[i+1] = angle_constrains[i+1];
            }

          // now find the controll point: where the two lines intersect:
          // p0.x
          // if angles are the same, the control point is either on the line or at infinity. Put it at P[i] for now
          //WBAssert((angles[i]-angles[i+1]) < std::numeric_limits<double>::epsilon()*10., "angles are the same, so the lines are parallel. Can't find a ")
          std::cout << i << ": ";
          if (std::abs(fmod((fmod(angles[i]-angles[i+1],180.) + 180.), 180.)) < std::numeric_limits<double>::epsilon()*10.)
            {
              control_points[i] = p[i];
              std::cout << ", flag 1: "<< angles[i] << " (" << angles[i] *180./M_PI << "), "<< angles[i+1] << " (" << angles[i+1] *180./M_PI << ")";
            }
          else
            {
              std::cout << ", flag 2";
              const double &x0 = p[i][0];
              const double &y0 = p[i][1];
              const double &x1 = p[i+1][0];
              const double &y1 = p[i+1][1];
              if (std::abs(fmod((fmod(angles[i], 180.) + 180.), 180.) - 90.) < std::numeric_limits<double>::epsilon()*10.)
                {
                  std::cout << ", flag 2.1";
                  // vertical line at x = x0
                  control_points[i][0] = x0;
                  control_points[i][1] = std::tan(angles[i+1]) * (x0-x1) + y1;
                }
              else if (std::abs(fmod((fmod(angles[i+1], 180.) + 180.), 180.) - 90.) < std::numeric_limits<double>::epsilon()*10.)
                {
                  std::cout << ", flag 2.2";
                  // vertical line at x = x0
                  control_points[i][0] = x1;
                  control_points[i][1] = std::tan(angles[i]) * (x1-x0) + y0;
                }
              std::cout << ", flag 2.3: angles i: " << angles[i] << " (" << angles[i]*180./M_PI  << "), angles i+1: " << angles[i+1] << " (" << angles[i+1]*180./M_PI  << ")";
              const double m0 = std::tan(angles[i]); // Line 0: y = m0 (x - x0) + y0
              const double m1 = std::tan(angles[i+1]); // Line 1: y = m1 (x - x1) + y1
              const double control_x = ((m0 * x0 - m1 * x1) - (y0 - y1)) / (m0 - m1);
              std::cout << ", m0 = " << m0 << ", m1 = " << m1 << ", m0 - m1 = " << m0 - m1 << ", control point = " << control_x << ":" << m0 * (control_x - x0) + y0;
              control_points[i][0] = control_x;
              control_points[i][1] = m0 * (control_x - x0) + y0;

              double t = 0.25;
              auto pq = (1-t)*(1-t)*p[i] + 2*t*(1-t)*control_points[i] + t*t*p[i+1];
              t = 0.5;
              auto ph = (1-t)*(1-t)*p[i] + 2*t*(1-t)*control_points[i] + t*t*p[i+1];
              t = 0.75;
              auto p3q = (1-t)*(1-t)*p[i] + 2*t*(1-t)*control_points[i] + t*t*p[i+1];
              std::cout << "length = "
                        << sqrt((p[i][0]-pq[0])*(p[i][0]-pq[0])+(p[i][1]-pq[1])*(p[i][1]-pq[1]))
                        + sqrt((ph[0]-pq[0])*(ph[0]-pq[0])+(ph[1]-pq[1])*(ph[1]-pq[1]))
                        + sqrt((ph[0]-p3q[0])*(ph[0]-p3q[0])+(ph[1]-p3q[1])*(ph[1]-p3q[1]))
                        + sqrt((p[i+1][0]-p3q[0])*(p[i+1][0]-p3q[0])+(p[i+1][1]-p3q[1])*(p[i+1][1]-p3q[1])) << std::endl;

              // compute length of segment
              lengths[i] = arc_length(p[i],control_points[i],p[i+1]);

              std::cout << std::endl << " --> total arc length = " << arc_length(p[i],control_points[i],p[i+1])
                        << ", t=1 arc length = " << arc_length(p[i],control_points[i],p[i+1],1.)
                        << ", t=0.5 arc length = " << arc_length(p[i],control_points[i],p[i+1],0.5)
                        << ", t=0 arc length = " << arc_length(p[i],control_points[i],p[i+1],0.) << std::endl;

              //// Based on http://geomalgorithms.com/a02-_lines.html.
              //// if the control point is doesn't project on the line
              //// rotate is 180 degrees.
              //const Point<2> BSP_ESP = p[i+1] - p[i];
              //const Point<2> BSP_CP = control_points[i] - p[i];
//
              //const double c1 = BSP_ESP * BSP_CP;
              //const double c2 = BSP_ESP * BSP_ESP;
//
              //if (c1 < 0 || c2 < c1)
              //  {
              //      control_points[i][0] = 2 * x1 - control_points[i][0];
              //      control_points[i][1] = 2.*y1 - control_points[i][1];
              //  }
            }
          //std::cout << std::endl;
        }
      std::cout << "-----" << x_spline.m.size() << " -- " << n_points-1 << std::endl;
      if (std::isnan(angle_constrains[n_points-2]))
        {
          const double value_x = x_spline.m[n_points-2][3];
          const double value_y = y_spline.m[n_points-2][3];
          const double derivative_x = x_spline.m[n_points-2][2];
          const double derivative_y = y_spline.m[n_points-2][3];
          // have the angle be pointing the the previous halfpoint instead of the next one
          const double value_x_h = x_spline.m[n_points-2][0]*0.5*0.5*0.5 + x_spline.m[n_points-2][1]*0.5*0.5 + x_spline.m[n_points-2][2]*0.5 + x_spline.m[n_points-2][3];
          const double value_y_h = y_spline.m[n_points-2][0]*0.5*0.5*0.5 + y_spline.m[n_points-2][1]*0.5*0.5 + y_spline.m[n_points-2][2]*0.5 + y_spline.m[n_points-2][3];
          angles[n_points-2] = std::atan2(value_y_h-value_y,value_x_h-value_x);//std::atan2(derivative_y,derivative_x);
          std::cout << n_points-1 << ": value_x:y = " << value_x << ":" << value_y << ",  value_h x:y = " << value_x_h << ":" << value_y_h << ", angle = " << angles[n_points-1] << " (" << angles[n_points-1] *180./M_PI << ")" << ", derivative_x = " << derivative_x
                    << ", x a:b:c:d = " << x_spline.m[n_points-2][0] << ":" << x_spline.m[n_points-2][1] << ":" << x_spline.m[n_points-2][2]  << ":" << x_spline.m[n_points-2][1]
                    << ", derivative_y = " << derivative_y
                    << ", y a:b:c:d = " << y_spline.m[n_points-2][0] << ":" << y_spline.m[n_points-2][1] << ":" << y_spline.m[n_points-2][2]  << ":" << y_spline.m[n_points-2][1]
                    << std::endl;
          // compute p1 and p2 relative to p0 to compute the angles.
          //angles[0] = std::atan2(p[1][1]-p[0][1],p[1][0]-p[0][0]);
          //const double angle_p0p1 = std::atan2(p[1][1]-p[0][1],p[1][0]-p[0][0]);
          //const double angle_p0p2 = std::atan2(p[2][1]-p[0][1],p[2][0]-p[0][0]);
          //angles[0] = angle_p0p1-0.5*(angle_p0p2-angle_p0p1); //2.*angle_p0p1-angle_p0p2;//
          //std::cout << "angles[0] = " << angles[0] << " (" << angles[0] *180./M_PI << ")" << ", angle_p0p1 = " << angle_p0p1  << "(" << angle_p0p1 *180./M_PI << ")"<< ", angle_p0p2 = " << angle_p0p2  << "(" << angle_p0p2 *180./M_PI << ")"<< std::endl;
        }
      else
        {
          // NOTE: start angle assumes slabs or faults going down, which means they should provide a negative angle to get expected behavoir
          angles[n_points-2] = angle_constrains[n_points-2];
        }

      if (std::isnan(angle_constrains[n_points-1]))
        {
          const double value_x = x_spline.m[n_points-2][0] + x_spline.m[n_points-2][1] + x_spline.m[n_points-2][2] + x_spline.m[n_points-2][3];
          const double value_y = y_spline.m[n_points-2][0] + y_spline.m[n_points-2][1] + y_spline.m[n_points-2][2] + y_spline.m[n_points-2][3];
          const double derivative_x = x_spline.m[n_points-2][0] + x_spline.m[n_points-2][1] + x_spline.m[n_points-2][2];
          const double derivative_y = y_spline.m[n_points-2][0] + y_spline.m[n_points-2][1] + y_spline.m[n_points-2][3];
          // have the angle be pointing the the previous halfpoint instead of the next one
          const double value_x_h = x_spline.m[n_points-2][0]*0.5*0.5*0.5 + x_spline.m[n_points-2][1]*0.5*0.5 + x_spline.m[n_points-2][2]*0.5 + x_spline.m[n_points-2][3];
          const double value_y_h = y_spline.m[n_points-2][0]*0.5*0.5*0.5 + y_spline.m[n_points-2][1]*0.5*0.5 + y_spline.m[n_points-2][2]*0.5 + y_spline.m[n_points-2][3];
          angles[n_points-1] = std::atan2(value_y_h-value_y,value_x_h-value_x);//std::atan2(derivative_y,derivative_x);
          std::cout << n_points << ": value_x:y = " << value_x << ":" << value_y << ",  value_h x:y = " << value_x_h << ":" << value_y_h << ", angle = " << angles[n_points-1] << " (" << angles[n_points-1] *180./M_PI << ")" << ", derivative_x = " << derivative_x
                    << ", x a:b:c:d = " << x_spline.m[n_points-2][0] << ":" << x_spline.m[n_points-2][1] << ":" << x_spline.m[n_points-2][2]  << ":" << x_spline.m[n_points-2][1]
                    << ", derivative_y = " << derivative_y
                    << ", y a:b:c:d = " << y_spline.m[n_points-2][0] << ":" << y_spline.m[n_points-2][1] << ":" << y_spline.m[n_points-2][2]  << ":" << y_spline.m[n_points-2][1]
                    << std::endl;
          // compute p1 and p2 relative to p0 to compute the angles.
          //angles[0] = std::atan2(p[1][1]-p[0][1],p[1][0]-p[0][0]);
          //const double angle_p0p1 = std::atan2(p[1][1]-p[0][1],p[1][0]-p[0][0]);
          //const double angle_p0p2 = std::atan2(p[2][1]-p[0][1],p[2][0]-p[0][0]);
          //angles[0] = angle_p0p1-0.5*(angle_p0p2-angle_p0p1); //2.*angle_p0p1-angle_p0p2;//
          //std::cout << "angles[0] = " << angles[0] << " (" << angles[0] *180./M_PI << ")" << ", angle_p0p1 = " << angle_p0p1  << "(" << angle_p0p1 *180./M_PI << ")"<< ", angle_p0p2 = " << angle_p0p2  << "(" << angle_p0p2 *180./M_PI << ")"<< std::endl;
        }
      else
        {
          // NOTE: start angle assumes slabs or faults going down, which means they should provide a negative angle to get expected behavoir
          angles[n_points-1] = angle_constrains[n_points-1];
        }
      size_t i = n_points-2;
      std::cout << i << ": ";
      if (std::abs(fmod((fmod(angles[i]-angles[i+1],180.) + 180.), 180.)) < std::numeric_limits<double>::epsilon()*10.)
        {
          control_points[i] = p[i];
          std::cout << ", flag 1: "<< angles[i] << " (" << angles[i] *180./M_PI << "), "<< angles[i+1] << " (" << angles[i+1] *180./M_PI << ")";
        }
      else
        {
          std::cout << ", flag 2";
          const double &x0 = p[i][0];
          const double &y0 = p[i][1];
          const double &x1 = p[i+1][0];
          const double &y1 = p[i+1][1];
          if (std::abs(fmod((fmod(angles[i], 180.) + 180.), 180.) - 90.) < std::numeric_limits<double>::epsilon()*10.)
            {
              std::cout << ", flag 2.1";
              // vertical line at x = x0
              control_points[i][0] = x0;
              control_points[i][1] = std::tan(angles[i+1]) * (x0-x1) + y1;
            }
          else if (std::abs(fmod((fmod(angles[i+1], 180.) + 180.), 180.) - 90.) < std::numeric_limits<double>::epsilon()*10.)
            {
              std::cout << ", flag 2.2";
              // vertical line at x = x0
              control_points[i][0] = x1;
              control_points[i][1] = std::tan(angles[i]) * (x1-x0) + y0;
            }
          std::cout << ", flag 2.3: angles i: " << angles[i] << " (" << angles[i]*180./M_PI  << "), angles i+1: " << angles[i+1] << " (" << angles[i+1]*180./M_PI  << ")";
          const double m0 = std::tan(angles[i]); // Line 0: y = m0 (x - x0) + y0
          const double m1 = std::tan(angles[i+1]); // Line 1: y = m1 (x - x1) + y1
          const double control_x = ((m0 * x0 - m1 * x1) - (y0 - y1)) / (m0 - m1);
          std::cout << ", m0 = " << m0 << ", m1 = " << m1 << ", m0 - m1 = " << m0 - m1 << ", control point = " << control_x << ":" << m0 * (control_x - x0) + y0;
          control_points[i][0] = control_x;
          control_points[i][1] = m0 * (control_x - x0) + y0;

          //// Based on http://geomalgorithms.com/a02-_lines.html.
          //// if the control point is doesn't project on the line
          //// rotate is 180 degrees.
          //const Point<2> BSP_ESP = p[i+1] - p[i];
          //const Point<2> BSP_CP = control_points[i] - p[i];
//
          //const double c1 = BSP_ESP * BSP_CP;
          //const double c2 = BSP_ESP * BSP_ESP;
//
          //if (c1 < 0 || c2 < c1)
          //  {
          //      control_points[i][0] = 2 * x1 - control_points[i][0];
          //      control_points[i][1] = 2.*y1 - control_points[i][1];
          //  }
        }
      //std::cout << std::endl;

      // compute length of segment
      lengths[n_points-2] = arc_length(p[n_points-2],control_points[n_points-2],p[n_points-1]);
      std::cout << p[n_points-2] << " - " << control_points[n_points-2] << " - " << p[n_points-1] << std::endl;


    }


    Point<2>
    BezierCurve::operator()(const size_t i, const double t) const
    {
      return (1-t)*(1-t)*points[i] + 2*t*(1-t)*control_points[i] + t*t*points[i+1];
    }

    ClosestPointOnCurve
    BezierCurve::closest_point_on_curve(const Point<2> &check_point) const
    {
      // go through each section and find all roots in domain 0 to 1 and choose the smallest
      ClosestPointOnCurve closest_point_on_curve;
      closest_point_on_curve.distance = std::numeric_limits<double>::infinity();
      closest_point_on_curve.fraction = std::numeric_limits<double>::signaling_NaN();
      for ( size_t i = 0; i < control_points.size(); ++i)
        {
          // compute a,b,c and d in the cubic equation describing the distance from point p to the local quadratic Bezier curve.
          // using https://blog.gludion.com/2009/08/distance-to-quadratic-bezier-curve.html
          // todo: I should also take a look at: https://iquilezles.org/articles/distfunctions2d/
          const Point<2> A = control_points[i]-points[i];
          const Point<2> B = points[i+1]-control_points[i]-A;
          const double a = B*B;
          const double b = 3*A*B;
          const double c = 2.*A*A+(points[i]-check_point)*B;
          const double d = (points[i]-check_point)*A;
          const std::vector<double> real_roots = this->solve_cubic_equation_real(a,b,c,d);
          //std::cout << "real roots: ";
          for (size_t root_i = 0; root_i < real_roots.size(); ++root_i)
            {
              if (real_roots[root_i] >= 0. && real_roots[root_i] <= 1.)
                {
                  const Point<2> point_on_curve = (*this)(i,real_roots[root_i]);
                  //std::cout << root_i << " = " << real_roots[root_i] << ", poc=" << point_on_curve << ", d = " << point_on_curve.distance(check_point) << "; ";
                  if (point_on_curve.distance(check_point) < std::abs(closest_point_on_curve.distance))
                    {
                      // the sign is rotating the derivative by 90 degrees.
                      // When moving in the direction of increasing t, left is positve and right is negative.
                      const double &t = real_roots[root_i];
                      const Point<2> derivative_point = ((-2+2*t)*points[i] + (2-4*t)*control_points[i] + 2*t*points[i+1]);
                      const Point<2> second_derivative_point = ((2)*points[i] + (-4)*control_points[i] + 2*points[i+1]);
                      Point<2> tangent_point = derivative_point - point_on_curve;
                      // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
                      const double dot_product = (tangent_point*(check_point-point_on_curve));
                      const double sign = dot_product < 0. ? -1. : 1.;
                      tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());
                      //std::cout << "derivative point = " <<  derivative_point  << ", velocity: " << derivative_point.norm() << ", second_derivative_point = " << second_derivative_point
                      //          << ", tangent point? = " << tangent_point << " -- " << tangent_point+point_on_curve << ", sign = " << sign << ", dot product = " << dot_product << ",  distance before: " << point_on_curve.distance(check_point) << std::endl;
                      closest_point_on_curve.distance = sign*point_on_curve.distance(check_point);
                      closest_point_on_curve.fraction = real_roots[root_i];
                      closest_point_on_curve.index = i;
                      closest_point_on_curve.point = point_on_curve;
                    }
                }
              else
                {
                  //std::cout << root_i << " = " << real_roots[root_i] << "; ";
                }
            }
          //std::cout << std::endl ;
        }
      return closest_point_on_curve;
    }


    std::vector<double>
    BezierCurve::solve_cubic_equation_real(const double a_original,const double b_original,const double c_original,const double d_original)
    {

      std::vector<double> real_roots;
      if (std::abs(a_original) <= std::numeric_limits<double>::epsilon()*10.)
        {
          const double &a = b_original;
          const double &b = c_original;
          const double &c = d_original;
          const double discriminant = b*b -4.*a*c;
          if (std::abs(discriminant) <= std::numeric_limits<double>::epsilon()*10.)
            {
              real_roots.emplace_back((-b+sqrt(discriminant))/(2*a));
              return real_roots;
            }
          else if (discriminant > 0)
            {
              real_roots.emplace_back((-b + sqrt(discriminant))/(2.*a));
              real_roots.emplace_back((-b - sqrt(discriminant))/(2.*a));
              return real_roots;
            }
          return real_roots;
        }
      else
        {
          // based on https://quarticequations.com/Cubic.pdf
          // divide by a
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
              real_roots.emplace_back(t-b/3.);
            }
          else
            {
              // three real solutions
              const double theta = std::abs(q) <= std::numeric_limits<double>::epsilon()*10. ? 0 : acos(r/std::pow(-q,3./2.));
              const double phi_1 = theta/3.;
              const double phi_2 = phi_1 - 2.*Consts::PI/3.;
              const double phi_3 = phi_1 + 2.*Consts::PI/3.;
              const double sqrt_q_3 = 2*sqrt(-q);
              const double value_1 = sqrt_q_3 *cos(phi_1)-b/3.;
              const double value_2 = sqrt_q_3 *cos(phi_2)-b/3.;
              const double value_3 = sqrt_q_3 *cos(phi_3)-b/3.;

              real_roots.emplace_back(value_1);
              if (std::abs(value_1 - value_2) > std::numeric_limits<double>::epsilon()*10.)
                {
                  real_roots.emplace_back(value_2);
                }
              // we don't have to check value 1 and 3 because z3 <= z2 <= z1
              // so if z2 and z3 are not equal, z1 and z3 are not either.
              if (std::abs(value_2 - value_3) > std::numeric_limits<double>::epsilon()*10.)
                {
                  real_roots.emplace_back(value_3);
                }
            }
        }
      return real_roots;
    }

  }
}