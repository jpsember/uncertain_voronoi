package gvorn;

import base.*;
import testbed.*;

/**
 */
public final class PlaneCurve implements Globals {

  /**
   * Construct a curve with x, y coefficients exchanged.
   * @return Polyn
   */
  private PlaneCurve exchangeVars() {
    return new PlaneCurve(c(3), c(4), c(5), c(1), c(2), c(0));
  }

  public double c(int index) {
    return poly.c(index);
  }

  /**
   * Construct a curve for a line
   * @param s0 endpoint #0
   * @param s1 endpoint #1
   */
  public PlaneCurve(FPoint2 s0, FPoint2 s1) {
    final double[] args = new double[6];
    args[2] = s1.y - s0.y;
    args[1] = s0.x - s1.x;
    args[0] = s0.x * (s0.y - s1.y) + s0.y * (s1.x - s0.x);
    poly = new Polyn(args);
  }

  /**
   * Construct a curve for a conic
   *
   * @param xSquared : coefficient for x*x, etc.
   * @param xy double
   * @param ySquared double
   * @param x double
   * @param y double
   * @param constant double
   * @return Polyn
   */
  public PlaneCurve(double xSquared, double xy, double ySquared, double x,
      double y, double constant) {
    final double[] args = new double[6];
    args[5] = xSquared;
    args[4] = xy;
    args[3] = ySquared;
    args[2] = x;
    args[1] = y;
    args[0] = constant;
    poly = new Polyn(args);
    //  simplifyConicToLine();
  }

  //  /**
  //   * Detect if conic is a line, and if so, change to a line
  //   */
  //  private void simplifyConicToLine() {
  //    do {
  //      // see notes March 22, (1)
  //
  //      // If conic is actually a line, then we must have
  //      //      D^2 - 4AF = 0
  //      // and  E^2 - 4CF = 0
  //
  //      double A = c(5), C = c(3), D = c(2), E = c(1), F = c(0);
  //
  //      final double EPS = 1e-10;
  //      if (Math.abs(D * D - 4 * A * F) > EPS)
  //        break;
  //      if (Math.abs(E * E - 4 * C * F) > EPS)
  //        break;
  //      if (Math.abs(A) <= EPS) {
  //        poly = new Polyn(0, 1, -E / (2 * C));
  //      } else if (Math.abs(C) <= EPS) {
  //        poly = new Polyn(1, 0, -D / (2 * A));
  //      } else {
  //
  //      }
  //
  //    } while (false);
  //
  //  }

  public static void main(String[] args) {
    try {
      if (true) {
        PlaneCurve c = new PlaneCurve(3, 1, -2, 3, -7, 15);
        final DArray pts = new DArray();
        c.calcExtremeY(pts);
        Streams.out.println(pts);
      } else {
        final PlaneCurve[] curves = { new PlaneCurve(3, 1, -2, 3, -7, 15),
            new PlaneCurve(2, -2, 1, -4, 3, -6), };

        DArray a = new DArray();
        for (int i = 0; i < curves.length; i += 2) {
          try {
            findIntersect(curves[i], curves[i + 1], a);
          } catch (TBError e) {
            System.out.println(e.toString());
          }
        }
      }
    } catch (TBError e) {
      System.out.println(e.toString());
    }
  }

  public String toString() {
    return toString(false);
  }

  public String toString(boolean allDigits) {
    StringBuilder sb = new StringBuilder();
    sb.append("PlaneCurve ");

    if (allDigits) {
      sb.append("\n");
    }

    for (int i = 5; i >= 0; i--) {
      if (allDigits) {
        if (i < 5) {
          sb.append("\n + ");
        } else {
          sb.append("   ");
        }
      } else {
        if (i < 5) {
          sb.append(" + ");
        }
      }
      if (allDigits) {
        String dg = Double.toString(c(i));
        if (dg.charAt(0) != '-') {
          sb.append(' ');
        }
        sb.append(dg);
      } else {
        sb.append(Tools.f(c(i), 8, 2));
      }

      sb.append(' ');

      {
        switch (i) {
        case 0:
          break;
        case 1:
          sb.append('y');
          break;
        case 2:
          sb.append('x');
          break;
        case 3:
          sb.append("y2");
          break;
        case 4:
          sb.append("xy");
          break;
        case 5:
          sb.append("x2");
          break;
        }

      }
    }
    return sb.toString();
  }

  /**
   * Find points of intersection between curve and horizontal line
   * @param y : y-coordinate of line
   * @param ipts : intersection points returned here
   */
  public void findHorizontalIntersect(double y, DArray ipts) {
    ipts.clear();

    // Solve curve for y, producing quadratic in x
    Polyn p = solveForY(y);

    DArray r = new DArray();
    p.solve(r);

    for (int i = 0; i < r.size(); i++) {
      double x = r.getDouble(i);
      ipts.add(new FPoint2(x, y));
    }
  }

  /**
   * Find intersection points of two curves.
   * @param c1 : first curve
   * @param c2 : second curve
   * @param ipts : FPoint2 objects are returned in this array
   */
  public static void findIntersect(PlaneCurve c1, PlaneCurve c2, DArray ipts) {
    final boolean db = false;
    if (db) {
      System.out.println("findIntersect for\nc1=" + c1.toString(true) + "\nc2="
          + c2.toString(true));
    }
    ipts.clear();

    for (int pass = 0; pass < 3; pass++) {
      try {
        boolean swapResults = false;
        do {

          // both y^2 not zero:

          // if ysquared term is not near zero in either curve,
          // we can use the long method above.

          if (!(Polyn.nearZero(c1.c(3)) || Polyn.nearZero(c2.c(3)))) {
            if (db) {
              System.out.println(" Neither C1 nor C2 is near zero;");
            }
            intersectQuadratics(c1, c2, ipts);
            break;
          }

          if (db) {
            System.out.println("c1=\n" + c1.toString(true) + "\nc2=\n"
                + c2.toString(true));
          }

          // both x^2 not zero:

          // try xsquared term.  If nonzero, we can exchange x & y.
          if (!(Polyn.nearZero(c1.c(5)) || Polyn.nearZero(c2.c(5)))) {
            if (db) {
              System.out.println("... long method, exchanging vars");
            }
            PlaneCurve d1 = c1.exchangeVars();
            PlaneCurve d2 = c2.exchangeVars();
            intersectQuadratics(d1, d2, ipts);
            swapResults = true;
            break;
          }
          if (db) {
            System.out.println("c1 degree=" + c1.poly.degree() + ", c2="
                + c2.poly.degree());
          }

          // all x^2, xy, y^2 zero:
          if (c1.poly.degree() <= 2 && c2.poly.degree() <= 2) {
            double D = c1.c(2), E = c1.c(1), F = c1.c(0);
            double J = c2.c(2), K = c2.c(1), L = c2.c(0);

            double numer = J * F - D * L, denom = D * K - J * E;
            if (db) {
              System.out.println("numer=" + numer + "\ndenom=" + denom);
            }
            if (!Polyn.nearZero(denom)) {
              double y = numer / denom;
              double x = (-E * y - F) / D;
              if (db) {
                System.out.println(" x=" + x + "\n y=" + y);
              }
              ipts.add(new FPoint2(x, y));
            }
            break;
          }

          if (c1.poly.degree() <= 2) {
            // solve c1 for x or y
            if (!(Polyn.nearZero(c1.c(2)))) {
              // get expression for x in first eqn
              double D = c1.c(2), E = c1.c(1), F = c1.c(0);
              double j = -F / D, k = -E / D;
              // get quadratic for y from second eqn
              double a = c2.c(5), b = c2.c(4), c = c2.c(3), d = c2.c(2), e = c2
                  .c(1), f = c2.c(0);
              Polyn q = new Polyn(a * k * k + b * k + c, 2 * a * j * k + b * j
                  + d * k + e, a * j * j + d * j + f);
              DArray rts = new DArray();
              q.solve(rts);
              for (int i = 0; i < rts.size(); i++) {
                double y = rts.getDouble(i);
                double x = j + k * y;
                ipts.add(new FPoint2(x, y));
              }
            }
            break;
          } else if (c2.poly.degree() <= 2) {
            PlaneCurve d1 = c1.exchangeVars();
            PlaneCurve d2 = c2.exchangeVars();
            intersectQuadratics(d1, d2, ipts);
            swapResults = true;
            break;
          }

          throw new FPError("short method not supported yet for\n"
              + c1.toString(true) + "\nand\n" + c2.toString(true));
        } while (false);

        if (swapResults) {
          // swap x,y coordinates of answers
          for (int i = 0; i < ipts.size(); i++) {
            FPoint2 pt = (FPoint2) ipts.get(i);
            pt.setLocation(pt.y, pt.x);
          }
        }
      } catch (FPError err) {
        if (pass == 0) {
          PlaneCurve tmp = c1;
          c1 = c2;
          c2 = tmp;
          continue;
        } else if (pass == 1) {

          // use new solver
          PolynOfPolyn p1 = c1.ppolyn(), p2 = c2.ppolyn();
          PolynOfPolyn.solve(p1, p2, ipts);
          break;
        }
        {
          System.out.println("c1=" + c1.toString(true) + "\nc2="
              + c2.toString(true));
          throw new FPError("Failed on pass two:\n" + c1 + "\n" + c2 + "\n"
              + err);
        }
      }
      break;
    }
  }

  /**
   * Find intersection points of two curves
   * @param p1 : first curve
   * @param p2 : second curve
   * @param ipts : (must be empty!): intersection points
   *  are returned here
   */
  private static DArray intersectQuadratics(PlaneCurve p1, PlaneCurve p2,
      DArray ipts) {

    final boolean db = false;

    if (ipts == null)
      ipts = new DArray();
    ipts.clear();

    if (db && true) {
      System.out.println("intersectQuadratics:");
      System.out.println("p1=" + p1.toString(true));
      System.out.println("p2=" + p2.toString(true));
    }

    // give the coefficients symbolic names
    double A1 = p1.c(5), B1 = p1.c(4), C1 = p1.c(3), D1 = p1.c(2), E1 = p1.c(1), F1 = p1
        .c(0);
    double A2 = p2.c(5), B2 = p2.c(4), C2 = p2.c(3), D2 = p2.c(2), E2 = p2.c(1), F2 = p2
        .c(0);

    // Determine symbolic coefficients for change of variables
    // to express y in terms of y'.

    // ************************************************
    // C1 must not be near zero for this to be accurate.
    // ************************************************

    //    Polyn.warnIfNearZero(C1, "C1");

    double G = 1 / (2 * C1);
    double H = -B1 / (2 * C1);
    double J = -E1 / (2 * C1);

    // determine symbolic coefficients for the first curve
    // expressed in terms of x and y', by substituting the
    // expression for y in terms of y' into the second curve equation

    double K = A2 + H * (B2 + C2 * H);
    double L = G * (B2 + 2 * C2 * H);
    double M = C2 * G * G;
    double N = J * (B2 + 2 * C2 * H) + D2 + E2 * H;
    double P = G * (2 * C2 * J + E2);
    double Q = J * (C2 * J + E2) + F2;

    //                           2
    // coefficients to express y'   in terms of x:

    double R = B1 * B1 - 4 * A1 * C1, S = 2 * B1 * E1 - 4 * C1 * D1, T = E1
        * E1 - 4 * C1 * F1;

    // determine symbol coefficients to express y' in terms of x.
    double U = -(K + M * R), V = -(M * S + N), W = -(M * T + Q);

    // get coefficients to express y in terms of x.
    double a = U - B1 * L, b = V - B1 * P - E1 * L, d = W - E1 * P, e = 2 * C1
        * L, f = 2 * C1 * P;

    // come up with another equation for y, by subtracting
    // a multiple of the second curve from the first (to get
    // rid of the y squared term) and solving for y.

    // ************************************************
    // C2 must not be near zero for this to be accurate
    // ************************************************

    Polyn xp = new Polyn(A1 * e * e + a * B1 * e + C1 * a * a,

    2 * A1 * e * f + a * f * B1 + b * e * B1 + 2 * C1 * a * b + D1 * e * e + E1
        * a * e,

    A1 * f * f + b * f * B1 + d * e * B1 + C1 * b * b + 2 * C1 * a * d + 2 * D1
        * e * f + E1 * (a * f + b * e) + F1 * e * e,

    d * f * B1 + 2 * C1 * b * d + D1 * f * f + E1 * (b * f + d * e) + 2 * e * f
        * F1,

    C1 * d * d + E1 * d * f + F1 * f * f);

    if (db) {
      System.out.println(" quartic polynomial for x =\n" + xp.toString(true));
    }

    // find roots of the polynomial for x.  These become
    // candidates for the intersection points

    DArray xRoots = new DArray();
    xp.solve(xRoots);
    if (db) {
      System.out.println(" # candidate roots=" + xRoots.size());
    }

    //    // candidate x,y points to be filtered later
    //      DArray cdPoints = new DArray();
    //    cdPoints.clear();

    for (int ind = 0; ind < xRoots.size(); ind++) {
      double x = xRoots.getDouble(ind);

      if (db) {
        System.out.println(" candidate x = " + x);
      }

      // We have two formulas for y:
      //
      // (16)  y = (a x2 + bx + d) / (ex + f)
      //
      // (18)  y = (h x2 + jx + n) / (ix + m)
      //
      // They seem to both have zero denominators at the same time.
      //
      // So stick with (16); if denominator is zero, then
      // convert curve into a quadratic by plugging in x as a constant,
      // then solve the resulting quadratic in y.

      double denom = e * x + f;

      if (Polyn.isZero(denom)) {
        if (db) {
          System.out.println("denominator is zero; e=" + e + ", f=" + f
              + ", x=" + x);
        }

        // construct a quadratic in y.
        Polyn qy = new Polyn(C1, B1 * x + E1, (A1 * x + D1) * x + F1);

        if (db) {
          System.out.println("quadratic in y =\n" + qy.toString(true));
        }
        final DArray qSoln = new DArray();
        qy.solve(qSoln);
        ptLoop: for (int i = 0; i < qSoln.size(); i++) {
          FPoint2 pt = new FPoint2(x, qSoln.getDouble(i));
          if (db) {
            System.out.println("  root #" + i + " pt=" + pt);
          }

          // don't add the point if it's essentially the same as another.
          for (int j = 0; j < ipts.size(); j++) {
            FPoint2 tp = (FPoint2) ipts.get(j);
            if (db) {
              System.out.println("  comparing to previous " + tp);
            }
            if (Math.abs(pt.x - tp.x) + Math.abs(pt.y - tp.y) < .002) {
              if (db) {
                System.out.println(" same as root #" + j + ", skipping");
              }
              continue ptLoop;
            }
          }
          ipts.add(pt);
        }
        //        continue;
      } else {
        //      double numer = ( (a * x) + b) * x + d;
        double numer = a * x * x + b * x + d;
        double y = numer / denom;
        if (db) {
          System.out.println("           y = " + y);
        }
        ipts.add(new FPoint2(x, y));
      }
    }
    return ipts;
    //    // Now process each candidate x,y point; remove duplicates
    //
    //    if (db) {
    //      System.out.println("# candidate points=" + cdPoints.size());
    //    }
    //    if (cdPoints.size() > 4) {
    //      throw new FPError(
    //          "*Too many candidate points; quadratic filtering problem");
    //    }
    //
    //    for (int i = 0; i < cdPoints.size(); i++) {
    //      FPoint2 pt = (FPoint2) cdPoints.get(i);
    //      if (db) {
    //        System.out.println(" pt=" + pt);
    //      }
    //      double x = pt.x, y = pt.y;
    //
    //      // evaluate x & y in both curves and verify that it is zero
    //      double eval1 = p1.eval(x, y);
    //      double eval2 = p2.eval(x, y);
    //
    //      if (db) {
    //        System.out.println("  eval on c1 = " + eval1);
    //        System.out.println("  eval on c2 = " + eval2);
    //      }
    //
    //      // try a new approach: don't eliminate candidates based on how
    //      // far their roots are from zero.  The user can filter these
    //      // if so desired.
    //      if (false) {
    //        if (!Polyn.isZero(eval1, .005) || !Polyn.isZero(eval2, .005)) {
    //          if (db) {
    //            System.out.println(" not sufficiently close to zero...");
    //          }
    //          continue;
    //        }
    //      }
    //      if (db) {
    //        System.out.println("  ...adding this point");
    //      }
    //
    //      ipts.add(pt);
    //    }

  }

  /**
   * Evaluate point on a curve
   * @param x : x coordinate
   * @param y : y coordinate
   * @return double; should be zero for a point on the curve
   */
  public double eval(double x, double y) {
    return x * (c(5) * x + c(4) * y + c(2)) + y * (c(3) * y + c(1)) + c(0);
  }

  /**
   * Get quadratic in x by substituting a value for y
   * @param y : value to substitute for y
   * @return quadratic
   */
  public Polyn solveForY(double y) {
    return new Polyn(
    // x^2
        c(5),
        // x
        c(4) * y + c(2),
        // 1
        c(3) * y * y + c(1) * y + c(0));
  }

  /**
   * Get quadratic in y by substituting a value for x
   * @param x : value to substitute for x
   * @return quadratic
   */
  public Polyn solveForX(double x) {
    final boolean db = false;
    if (db) {
      System.out.println("solveForX, x=" + x);
    }
    if (db) {
      System.out.println(" c=" + toString(true));
    }
    Polyn p = new Polyn(
    // y^2
        c(3),
        // y
        c(4) * x + c(1),
        // 1
        (c(5) * x + c(2)) * x + c(0));
    if (db) {
      System.out.println(" p=" + p.toString(true));
    }
    return p;
  }

  /**
   * Calculate the extreme points of the curve, vertically;
   * these are the points where the slope of the curve is zero
   * @param pts : extreme points are returned here as FPoint2's
   */
  public void calcExtremeY(DArray pts) {
    calcExtreme(pts, true);
  }

  /**
   * Calculate the extreme points of the curve, horizontally;
   * these are the points where the slope of the curve is infinite
   * @param pts : extreme points are returned here as FPoint2's
   */
  public void calcExtremeX(DArray pts) {
    calcExtreme(pts, false);
  }

  /**
   * Calculate extreme points on the curve, for either x or y
   * @param pts : extreme points are returned here
   * @param yFlag : if true, calculates points where slope is zero; else
   *  where slope is infinite
   */
  private void calcExtreme(DArray pts, boolean yFlag) {
    final boolean db = false;
    if (db) {
      System.out.println("calcExtreme for\n" + this.toString(true));
    }
    pts.clear();

    double A, B, C, D, E, F;
    if (yFlag) {
      A = c(5);
      B = c(4);
      C = c(3);
      D = c(2);
      E = c(1);
      F = c(0);
    } else {
      A = c(3);
      B = c(4);
      C = c(5);
      D = c(1);
      E = c(2);
      F = c(0);
    }

    if (db) {
      System.out.println(" A=" + A + " B=" + B + " C=" + C + " D=" + D + " E="
          + E + " F=" + F);
    }

    boolean quadY = (Math.abs(B) < Math.abs(2 * A));

    Polyn q;

    if (db) {
      System.out.println(" Solving quadratic in " + (quadY ? "y" : "x"));
    }
    if (quadY) {
      q = new Polyn(4 * A * C - B * B, 4 * A * E - 2 * B * D, 4 * A * F - D * D);
    } else {
      q = new Polyn(A * (4 * A * C - B * B), 2 * A * (2 * C * D - B * E), C * D
          * D + B * (B * F - D * E));

    }

    if (db) {
      System.out.println(" q=" + q.toString(true));
    }

    DArray roots = new DArray();
    q.solve(roots);

    if (db) {
      System.out.println(" roots=" + Tools.d(roots));
    }

    for (int i = 0; i < roots.size(); i++) {
      double x, y;
      if (quadY) {
        y = roots.getDouble(i);
        x = -(B * y + D) / (2 * A);
      } else {
        x = roots.getDouble(i);
        y = (-2 * A * x - D) / B;
      }
      FPoint2 pt = yFlag ? new FPoint2(x, y) : new FPoint2(y, x);
      if (!pt.isValid()) {
        continue;
      }
      pts.add(pt);
      if (db) {
        System.out.println("soln " + i + " is " + pt);
      }
    }
  }

  /**
   */
  private PolynOfPolyn ppolyn() {
    if (ppoly == null) {
      ppoly = new PolynOfPolyn();
      ppoly.set(0, new Polyn(poly.c(3), poly.c(1), poly.c(0)));
      ppoly.set(1, new Polyn(poly.c(4), poly.c(2)));
      ppoly.set(2, new Polyn(poly.c(5)));
    }
    return ppoly;
  }

  /**
   * Find points of intersection between curve and line
   * @param p0 FPoint2 : first point on line
   * @param p1 FPoint2 : second point on line
   * @param tvals : array of t values of intersection returned here, where
   *   pt is p0 + t(p1-p0); sorted by t
   */
  public void findLineIntersect(FPoint2 p0, FPoint2 p1, DArray tvals) {
    // get into form
    //   Gx + Hy + J = 0
    double g = p1.x - p0.x, h = p1.y - p0.y;

    double G = h, H = -g, J = g * p0.y - p0.x * h;

    double A = c(5), B = c(4), C = c(3), D = c(2), E = c(1), F = c(0);

    DArray r = new DArray();
    if (Math.abs(H) < Math.abs(G)) {
      double K = -H / G, L = -J / G;

      double c2 = (A * K + B) * K + C, c1 = K * (2 * A * L + D) + B * L + E, c0 = L
          * (A * L + D) + F;
      //
      //      double c2 = B*K + C,
      //          c1 = A*K*K+2*A*K*L + B*L + D*K + E,
      //          c0= A*L*L+D*L+F;
      Polyn yp = new Polyn(c2, c1, c0);
      yp.solve(r);
      for (int i = 0; i < r.size(); i++) {
        double y = r.getDouble(i);
        tvals.addDouble((y - p0.y) / (p1.y - p0.y));
      }
    } else {
      double K = -G / H, L = -J / H;
      double c2 = (C * K + B) * K + A, c1 = K * (2 * C * L + E) + B * L + D, c0 = L
          * (C * L + E) + F;
      //      double c2 = B*K + A,
      //          c1 = C*K*K+2*C*K*L + B*L + E*K + D,
      //          c0= C*L*L+E*L+F;
      Polyn xp = new Polyn(c2, c1, c0);
      xp.solve(r);
      for (int i = 0; i < r.size(); i++) {
        double x = r.getDouble(i);
        tvals.addDouble((x - p0.x) / (p1.x - p0.x));
      }
    }

    for (int i = 0; i < tvals.size(); i++) {
      for (int j = i + 1; j < tvals.size(); j++) {
        if (tvals.getDouble(i) > tvals.getDouble(j))
          tvals.swap(i, j);
      }
    }

  }

  private Polyn poly;

  private PolynOfPolyn ppoly;
}
