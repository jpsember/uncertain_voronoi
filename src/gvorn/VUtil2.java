package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class VUtil2 {

  public static void setTraceArray(DArray trace) {
    VUtil2.trace = trace;
  }

  private static DArray trace;

  private static int trace(Object obj) {
    int ret = 0;
    if (trace != null) {
      ret = trace.size();
      trace.add(obj);
    }
    return ret;
  }
  private static void untrace(int i) {
    if (trace != null) {
      trace.removeRange(i, trace.size());
    }
  }
  /**
   * Construct Voronoi diagram
   * @param sites array of VornSites
   * @return array of VornBisectors representing edges of cells
   */
  public static DArray[] buildVornDiagram(DArray sites) {

    for (int i = 0; i < sites.size(); i++)
      ((VornSite) sites.get(i)).setId(i);

    DArray[] bList = buildBisectors(sites);

    for (int i = 0; i < sites.size(); i++) {
      buildVornCell((VornSite) sites.get(i), bList, i);
      if (false) {
        Tools.warn("skipping remaining cells");
        break;
      }
    }
    return bList;
  }

  /**
   * Build bisectors for Voronoi sites
   * @param sites array of VornSites
   * @return array of bisector lists; list[i] is a DArray of VornBisectors 
   *  associated with site #i
   */
  private static DArray[] buildBisectors(DArray sites) {
    DArray[] bList = new DArray[sites.size()];
    for (int i = 0; i < bList.length; i++)
      bList[i] = new DArray();

    for (int si = 0; si < sites.size(); si++) {
      VornSite siteI = (VornSite) sites.get(si);
      for (int sj = si + 1; sj < sites.size(); sj++) {
        VornSite siteJ = (VornSite) sites.get(sj);

        DArray bl = siteI.buildBisector(siteJ);
        for (int i = 0; i < bl.size(); i++) {
          VornBisector bis = (VornBisector) bl.get(i);
          bList[si].add(bis);
          bList[sj].add(bis);
        }
      }
    }
    return bList;
  }

  /**
   * Construct Voronoi cell
   * @param bList array of bisector lists
   * @param siteIndex index of cell to construct
   */
  private static void buildVornCell(VornSite site, DArray[] bList, int siteIndex) {

    final boolean db = true;
    final boolean db2 = false;

    if (db && T.update())
      T.msg("buildVornCell  siteIndex=" + siteIndex
          + T.show(site, MyColor.cRED));

    DArray bl = bList[siteIndex];
    for (int i = 0; i < bl.size(); i++) {
      VornBisector bi = (VornBisector) bl.get(i);

      VornSite iOther = bi.otherSite(site);

      int tr0 = trace(bi);
      trace(iOther);

      if (db2 && T.update())
        T.msg("bisector i#" + i);

      for (int j = i + 1; j < bl.size(); j++) {
        VornBisector bj = (VornBisector) bl.get(j);
        VornSite jOther = bj.otherSite(site);

        // if sites of both bisectors are the same, we
        // can assume no clipping needs doing at this level
        if (iOther == jOther)
          continue;

        int tr1 = trace(bj);
        trace(jOther);

        if (db2 && T.update())
          T.msg("bisector j#" + j);

        DArray isects = VUtil2.intersect(bi, bj, true);

        if (db2 && T.update()) {
          for (int q = 0; q < isects.size(); q += 2) {
            double t = isects.getDouble(q);
            T.show(bi.pointAt(t));
          }
          T.msg("isects=" + DArray.toString(isects.toDoubleArray()));
        }

        if (isects.size() == 0) {

          // One bisector may be clipping the other out completely.
          // Let pi, pj be arbitrary points from bisectors.

          // If   d(pi, sj) > d(pi, si)
          // and  d(pj, si) < d(pj, sj),  clip out j's bisector

          // If   d(pj, si) > d(pj, sj)
          // and  d(pi, sj) < d(pi, si),  clip out i's bisector

          if (bi.nComponents() > 0 && bj.nComponents() > 0) {

            if (db && T.update())
              T.msg("testing for obscuring" + T.show(bi) + T.show(bj));
            FPoint2 pi = arbPt(bi);
            FPoint2 pj = arbPt(bj);
            boolean flag1 = jOther.distanceFrom(pi) > iOther.distanceFrom(pi);
            boolean flag2 = iOther.distanceFrom(pj) > jOther.distanceFrom(pj);

            if (flag1 && !flag2) {
              if (db && T.update())
                T.msg("clipping all of " + bj);
              bj.clipAll();
            } else if (!flag1 && flag2) {
              if (db && T.update())
                T.msg("clipping all of " + bi);
              bi.clipAll();
            }
          }
        } else {

          for (int pass = 0; pass < 2; pass++) {
            VornBisector b1 = (pass == 0) ? bi : bj;
            VornBisector b2 = (pass == 0) ? bj : bi;
            DArray tv = new DArray();
            tv.addDouble(VornBisector.TMIN);
            tv.addDouble(VornBisector.TMAX);
            for (int k = pass; k < isects.size(); k += 2) {
              tv.addDouble(isects.getDouble(k));
            }
            tv.sort(DArray.COMPARE_DOUBLES);

            if (db2 && T.update())
              T.msg("pass=" + pass);

            for (int seg = 1; seg < tv.size(); seg++) {
              double t0 = tv.getDouble(seg - 1);
              double t1 = tv.getDouble(seg);

              //              final double RNG = 200;
              //              if (Math.min(Math.abs(t0), Math.abs(t1)) < RNG) {
              //                t0 = MyMath.clamp(t0, -RNG, RNG);
              //                t1 = MyMath.clamp(t1, -RNG, RNG);
              //
              //              }

              double t;
              {
                final double R = 60;
                double f0 = t0;
                double f1 = t1;

                if (f0 < -R && f1 > -R)
                  f0 = -R;
                if (f1 > R && f0 < R)
                  f1 = R;
                t = (f0 + f1) / 2;
              }

              //              double tStep = t0 - t1;
              //              tStep = MyMath.clamp(tStep, -4, 4);
              //
              //              double t = t1 + tStep;

              FPoint2 testPt = b1.pointAt(t);
              if (db2 && T.update())
                T.msg(T.show(b1.pointAt(t1)) + "test t=" + Tools.f(t)
                    + " isect t=" + Tools.f(t1));

              boolean clip = false;
              {

                double di = b2.otherSite(site).distanceFrom(testPt);
                double dj = site.distanceFrom(testPt);
                if (db2 && T.update())
                  T.msg("dist from site " + b2.otherSite(site) + "=" + di);

                if (db2 && T.update())
                  T.msg("dist from cell site  " + site + "=" + dj);
                //
                //                double scale = Math.max(di, dj);
                //                scale = Math.max(1, scale);
                //                if ((dj - di) / scale > 1e-6)
                //                  clip = true;

                if (dj - di > Math.abs(dj) * 1e-2) //Math.abs(tStep) / 60)
                  clip = true;
                if (db && T.update())
                  T.msg(T.show(b1) + T.show(b2) + "test point" + T.show(testPt)
                      + ", d(" + b2.otherSite(site) + ")=" + di + " d(" + site
                      + ")=" + dj + ", clip=" + clip);

              }
              if (clip) {
                if (db && T.update()) {
                  final double R = 100;
                  double cc0 = MyMath.clamp(t0, -R, R);
                  double cc1 = MyMath.clamp(t1, -R, R);
                  for (int q = 0; q < 100; q++)
                    T.show(b1.pointAt(cc0 + (cc1 - cc0) * (q / 100.0)));

                  T.msg("clipping out " + Tools.f(t0) + "..." + Tools.f(t1));
                }

                b1.clipOut(t0, t1);

                if (db2 && T.update()) {
                  T.msg("after clipping"
                      + T.show(b1, MyColor.cRED, Globals.STRK_THICK, -1)
                      + T.show(b1.pointAt(t1)));
                }

              }
            }
          }
        }
        untrace(tr1);
      }
      untrace(tr0);

    }
    if (db && T.update())
      T.msg("built voronoi cell" + T.show(bl));
  }

  private static FPoint2 arbPt(VornBisector b) {
    if (b.nComponents() == 0)
      throw new IllegalArgumentException();
    double t0 = b.cStart(0);
    double t1 = b.cEnd(0);
    double tSamp = MyMath.clamp(0, t0, t1);
    return b.pointAt(tSamp);
  }

  /**
   * Calculate intersection points of a pair of bisectors
   * @param b1 
   * @param b2
   * @param applyFilter if true, filters out points that are not equidistant to the
   * relevant sites
   * @return 2n doubles, representing parameters of intersection points
   *   of b1 and b2
   */
  public static DArray intersect(VornBisector b1, VornBisector b2,
      boolean applyFilter) {

    final boolean db = false;
    if (db && T.update())
      T.msg("intersect VornBisectors\nb1=" + b1 + "\nb2=" + b2
          + T.show(b1, MyColor.cRED, Globals.STRK_THICK, -1)
          + T.show(b2, MyColor.cPURPLE, Globals.STRK_THICK, -1));

    DArray iPts = new DArray();
    
    for (int i = 0; i < b1.nPieces(); i++) {
      IPlaneCurve ci = b1.curve(i);
      for (int j = 0; j < b2.nPieces(); j++) {
        IPlaneCurve cj = b2.curve(j);
        DArray qPts = intersectCurves(ci,cj);
        iPts.addAll(qPts);
        //b1.getPlaneCurve(), b2.getPlaneCurve());
      }
    }
    
    DArray jPts = new DArray();

    if (db && T.update())
      T.msg("returned " + iPts.size() + " candidates");

    // filter out false intersections
    for (int i = 0; i < iPts.size(); i++) {
      FPoint2 pt = iPts.getFPoint2(i);
      if (applyFilter) {
        if (db && T.update())
          T.msg("candidate #" + i + T.show(pt));

        if (!b1.onBisector(pt) || !b2.onBisector(pt)) {
          if (db && T.update())
            T.msg("not on one of the bisectors" + T.show(pt));
          continue;
        }
      }
      jPts.addDouble(b1.parameterFor(pt));
      jPts.addDouble(b2.parameterFor(pt));
    }
    return jPts;
    //
    //    
    //    
    //    DArray ip = findIntersections(b1, b2);
    //
    //    DArray ret = new DArray();
    //
    //    for (int i = 0; i < ip.size(); i++) {
    //      FPoint2 pt = ip.getFPoint2(i);
    //
    //      FPoint2 data = b1.calcParameterAndDistance(pt);
    //      ret.addDouble(data.x);
    //
    //      data = b2.calcParameterAndDistance(pt);
    //      ret.addDouble(data.x);
    //    }
    //    return ret;
  }

  /**
   * Construct array of FPoint2s corresponding to a list of intersection
   * parameters returned by intersect(...)
   * @param b1 curve
   * @param tPairs array of t parameter pairs; only first in each pair is
   *   assumed to correspond to b1
   * @return array of FPoint2s 
   */
  public static DArray calcPoints(IPlaneCurve b1, DArray tPairs) {
    DArray r = new DArray();
    for (int i = 0; i < tPairs.size(); i += 2)
      r.add(b1.pt(tPairs.getDouble(i)));

    return r;
  }

  private static DArray intersectCurves(IPlaneCurve b1, IPlaneCurve b2) {
    final boolean db = false;
    DArray ipts = new DArray();

    if (db && T.update())
      T.msg("findIntersect for\nb1=" + b1 + " deg=" + b1.degree() + " b2=" + b2
          + " deg=" + b2.degree());

    // are both curves conics?
    do {
      if (b1.degree() == 2) {
        if (b2.degree() == 2) {
          intersectQuadratics(b1, b2, ipts);
          if (db && T.update())
            T.msg("found ipts" + T.show(ipts));

          break;
        }
        intersectConicWithLine(b1, (LineEqn) b2, ipts);
        break;
      }
      if (b2.degree() == 2) {
        intersectConicWithLine(b2, (LineEqn) b1, ipts);
        break;
      }

      {
        FPoint2 ipt = LineEqn.intersection((LineEqn) b1, (LineEqn) b2);
        if (ipt != null)
          ipts.add(ipt);
        break;
      }
    } while (false);
    return ipts;
  }

  private static DArray intersectConicWithLine(IPlaneCurve conic, LineEqn line,
      DArray ipts) {

    final boolean db = false;

    if (db && T.update())
      T.msg("intersectConicWithLine");
    if (ipts == null)
      ipts = new DArray();
    ipts.clear();

    double L = line.coeff(2);
    double M = line.coeff(1);
    double P = line.coeff(0);

    double A = conic.coeff(5), B = conic.coeff(4), C = conic.coeff(3), D = conic
        .coeff(2), E = conic.coeff(1), F = conic.coeff(0);

    if (db && T.update()) {
      FPoint2 tp = new FPoint2(56, 35);
      double t = conic.parameterFor(tp);
      FPoint2 tp2 = conic.pt(t);

      double t3 = line.parameterFor(tp);
      FPoint2 tp3 = line.pt(t3);

      double x = 56;
      double y = 35;
      double tt = A * x * x + B * x * y + C * y * y + D * x + E * y + F;

      T.msg("line eqn L=" + L + "\n M=" + M + "\n P=" + P + "\n\nconic:\nA="
          + A + " B=" + B + " C=" + C + " D=" + D + " E=" + E + " F=" + F
          + "\n tp=" + tp + "\n t=" + t + "\n tp2=" + tp2 + "\n t3=" + t3
          + "\n tp3=" + tp3 + "\ntt=" + tt);

    }

    if (Math.abs(M) >= Math.abs(L)) {
      double Q = -L / M;
      double R = -P / M;

      double a = B * Q + A + C * Q * Q;
      double b = 2 * C * R * Q + D + E * Q + B * R;
      double c = F + C * R * R + E * R;

      Polyn q = new Polyn(a, b, c);
      DArray roots = new DArray();
      q.solve(roots, 1e-8);
      for (int i = 0; i < roots.size(); i++) {
        double x = roots.getDouble(i);
        FPoint2 isect = new FPoint2(x, Q * x + R);
        if (db && T.update())
          T.msg("isect pt" + T.show(isect));
        ipts.add(isect);
      }
    } else {

      double Q = -M / L;
      double R = -P / L;

      double a = A * Q * Q + B * Q + C;
      double b = B * R + E + 2 * A * R * Q + D * Q;
      double c = A * R * R + F + D * R;

      Polyn q = new Polyn(a, b, c);
      if (db && T.update()) {
        double rr = b * b - 4 * a * c;

        T.msg("solving quadratic to get y:\n" + q + "\n b^2-4ac=" + rr);

      }

      DArray roots = new DArray();
      q.solve(roots, 1e-8);
      for (int i = 0; i < roots.size(); i++) {
        double y = roots.getDouble(i);
        FPoint2 isect = new FPoint2(Q * y + R, y);
        if (db && T.update())
          T.msg("isect pt" + T.show(isect));
        ipts.add(isect);
      }

    }
    if (db && T.update())
      T.msg("returning " + ipts.size() + " points" + T.show(ipts));
    return ipts;
  }

  //  private static IPlaneCurve exchangeVars(IPlaneCurve c) {
  //    if (c.degree() < 2)
  //      throw new IllegalArgumentException();
  //
  //    Conic ch = (Conic) c;
  //    return new Conic(ch.coeff(3), ch.coeff(4), ch.coeff(5), ch.coeff(1), ch
  //        .coeff(2), ch.coeff(0));
  //
  //  }

  //  /**
  //   * Find intersection points of two curves.
  //   * @param c1 : first curve
  //   * @param c2 : second curve
  //   * @param ipts : FPoint2 objects are returned in this array
  //   */
  //  public static void findIntersect(IPlaneCurve c1, IPlaneCurve c2, DArray ipts) {
  //    final boolean db = true;
  //    if (db) {
  //      System.out.println("findIntersect for\nc1=" + c1 + "\nc2=" + c2);
  //    }
  //    ipts.clear();
  //
  //    for (int pass = 0; pass < 3; pass++) {
  //      try {
  //        boolean swapResults = false;
  //        do {
  //
  //          // both y^2 not zero:
  //
  //          // if ysquared term is not near zero in either curve,
  //          // we can use the long method above.
  //
  //          if (!(Polyn.nearZero(c1.coeff(3)) || Polyn.nearZero(c2.coeff(3)))) {
  //            if (db) {
  //              System.out.println(" Neither C1 nor C2 is near zero;");
  //            }
  //            intersectQuadratics(c1, c2, ipts);
  //            break;
  //          }
  //
  //          if (db) {
  //            System.out.println("c1=\n" + c1 + "\nc2=\n" + c2);
  //          }
  //
  //          // both x^2 not zero:
  //
  //          // try xsquared term.  If nonzero, we can exchange x & y.
  //          if (!(Polyn.nearZero(c1.coeff(5)) || Polyn.nearZero(c2.coeff(5)))) {
  //            if (db) {
  //              System.out.println("... long method, exchanging vars");
  //            }
  //            IPlaneCurve d1 = exchangeVars(c1);
  //            IPlaneCurve d2 = exchangeVars(c2);
  //            intersectQuadratics(d1, d2, ipts);
  //            swapResults = true;
  //            break;
  //          }
  //          if (db) {
  //            System.out.println("c1 degree=" + c1.degree() + ", c2="
  //                + c2.degree());
  //          }
  //
  //          // all x^2, xy, y^2 zero:
  //          if (c1.degree() < 2 && c2.degree() < 2) {
  //            double D = c1.coeff(2), E = c1.coeff(1), F = c1.coeff(0);
  //            double J = c2.coeff(2), K = c2.coeff(1), L = c2.coeff(0);
  //
  //            double numer = J * F - D * L, denom = D * K - J * E;
  //            if (db) {
  //              System.out.println("numer=" + numer + "\ndenom=" + denom);
  //            }
  //            if (!Polyn.nearZero(denom)) {
  //              double y = numer / denom;
  //              double x = (-E * y - F) / D;
  //              if (db) {
  //                System.out.println(" x=" + x + "\n y=" + y);
  //              }
  //              ipts.add(new FPoint2(x, y));
  //            }
  //            break;
  //          }
  //
  //          if (c1.degree() < 2) {
  //            // solve c1 for x or y
  //            if (!(Polyn.nearZero(c1.coeff(2)))) {
  //              // get expression for x in first eqn
  //              double D = c1.coeff(2), E = c1.coeff(1), F = c1.coeff(0);
  //              double j = -F / D, k = -E / D;
  //              // get quadratic for y from second eqn
  //              double a = c2.coeff(5), b = c2.coeff(4), c = c2.coeff(3), d = c2
  //                  .coeff(2), e = c2.coeff(1), f = c2.coeff(0);
  //              Polyn q = new Polyn(a * k * k + b * k + c, 2 * a * j * k + b * j
  //                  + d * k + e, a * j * j + d * j + f);
  //              DArray rts = new DArray();
  //              q.solve(rts);
  //              for (int i = 0; i < rts.size(); i++) {
  //                double y = rts.getDouble(i);
  //                double x = j + k * y;
  //                ipts.add(new FPoint2(x, y));
  //              }
  //            }
  //            break;
  //          } else if (c2.degree() < 2) {
  //            IPlaneCurve d1 = exchangeVars(c1);
  //            IPlaneCurve d2 = exchangeVars(c2);
  //            intersectQuadratics(d1, d2, ipts);
  //            swapResults = true;
  //            break;
  //          }
  //
  //          throw new FPError("short method not supported yet for\n" + c1
  //              + "\nand\n" + c2);
  //        } while (false);
  //
  //        if (swapResults) {
  //          // swap x,y coordinates of answers
  //          for (int i = 0; i < ipts.size(); i++) {
  //            FPoint2 pt = (FPoint2) ipts.get(i);
  //            pt.setLocation(pt.y, pt.x);
  //          }
  //        }
  //      } catch (FPError err) {
  //        System.out.println("failed on pass:"+pass);
  //        if (pass == 0) {
  //          IPlaneCurve tmp = c1;
  //          c1 = c2;
  //          c2 = tmp;
  //          continue;
  //        }
  //        //        else if (pass == 1) {
  //        //
  //        //          // use new solver
  //        //          PolynOfPolyn p1 = c1.ppolyn(), p2 = c2.ppolyn();
  //        //          PolynOfPolyn.solve(p1, p2, ipts);
  //        //          break;
  //        //        }
  //        {
  //          System.out.println("c1=" + c1 + "\nc2=" + c2);
  //          throw new FPError("Failed on pass two:\n" + c1 + "\n" + c2 + "\n"
  //              + err);
  //        }
  //      }
  //      break;
  //    }
  //  }

  /**
   * Find intersection points of two curves
   * @param p1 : first curve
   * @param p2 : second curve
   * @param ipts : (must be empty!): intersection points
   *  are returned here
   */
  private static DArray intersectQuadratics(IPlaneCurve p1, IPlaneCurve p2,
      DArray ipts) {

    final boolean db = false;

    if (ipts == null)
      ipts = new DArray();
    ipts.clear();

    if (db && T.update())
      T.msg("intersectQuadratics:\np1=" + p1 + "\np2=" + p2);

    // give the coefficients symbolic names
    double A1 = p1.coeff(5), B1 = p1.coeff(4), C1 = p1.coeff(3), D1 = p1
        .coeff(2), E1 = p1.coeff(1), F1 = p1.coeff(0);
    double A2 = p2.coeff(5), B2 = p2.coeff(4), C2 = p2.coeff(3), D2 = p2
        .coeff(2), E2 = p2.coeff(1), F2 = p2.coeff(0);

    // Determine symbolic coefficients for change of variables
    // to express y in terms of y'.

    // ************************************************
    // C1 must not be near zero for this to be accurate.
    // ************************************************

    if (db && T.update())
      T.msg("C1=" + C1);

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

    double R = B1 * B1 - 4 * A1 * C1, S = 2 * B1 * E1 - 4 * C1 * D1, TT = E1
        * E1 - 4 * C1 * F1;

    // determine symbol coefficients to express y' in terms of x.
    double U = -(K + M * R), V = -(M * S + N), W = -(M * TT + Q);

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

    if (db && T.update())
      T.msg("quartic polynomial for x =\n" + xp.toString(true));

    // find roots of the polynomial for x.  These become
    // candidates for the intersection points

    DArray xRoots = new DArray();
    xp.solve(xRoots);
    if (db && T.update())
      T.msg(" # candidate roots=" + xRoots.size());

    //    // candidate x,y points to be filtered later
    //      DArray cdPoints = new DArray();
    //    cdPoints.clear();

    for (int ind = 0; ind < xRoots.size(); ind++) {
      double x = xRoots.getDouble(ind);

      if (db && T.update())
        T.msg(" candidate x = " + x);

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
        if (db && T.update())
          T.msg("denominator is zero; e=" + e + ", f=" + f + ", x=" + x);

        // construct a quadratic in y.
        Polyn qy = new Polyn(C1, B1 * x + E1, (A1 * x + D1) * x + F1);

        if (db) {
          System.out.println("quadratic in y =\n" + qy.toString(true));
        }
        final DArray qSoln = new DArray();
        qy.solve(qSoln);
        ptLoop: for (int i = 0; i < qSoln.size(); i++) {
          FPoint2 pt = new FPoint2(x, qSoln.getDouble(i));
          if (db && T.update())
            T.msg("  root #" + i + " pt=" + pt);

          // don't add the point if it's essentially the same as another.
          for (int j = 0; j < ipts.size(); j++) {
            FPoint2 tp = (FPoint2) ipts.get(j);
            if (db) {
              System.out.println("  comparing to previous " + tp);
            }
            if (Math.abs(pt.x - tp.x) + Math.abs(pt.y - tp.y) < .002) {
              if (db && T.update())
                T.msg(" same as root #" + j + ", skipping");

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

        ipts.add(new FPoint2(x, y));
        if (db && T.update())
          T.msg("y = " + y + T.show(ipts.last()));
      }
    }
    return ipts;
  }

  public static void render(IPlaneCurve curve, double t0, double t1) {
    curve.render(t0, t1);
    //    
    //    if (curve instanceof LineEqn) {
    //      renderLine(curve.pt(t0), curve.pt(t1));
    //    } else {
    //      if (curve instanceof Hyperb)
    //        ((Hyperb) curve).render(t0, t1);
    //      else
    //        ((Parabola) curve).render(t0, t1);
    //    }
  }
  public static void renderLine(FPoint2 p0, FPoint2 p1) {
    //  V.pushStroke(stroke, Globals.STRK_NORMAL);
    //  V.pushColor(c, MyColor.cDARKGREEN);

    //    FPoint2 p0 = pointAt(t0);
    //    FPoint2 p1 = pointAt(t1);

    {
      if (MyMath.clipSegmentToRect(p0, p1, V.viewRect)) {
        V.drawLine(p0, p1);
        //        EdSegment.plotDirectedLine(p0, p1);

        if (false) {
          FPoint2 np = MyMath.ptOnCircle(FPoint2.midPoint(p0, p1), MyMath
              .polarAngle(p0, p1)
              + Math.PI / 2, 4);
          V.pushScale(.5);
          //   V.draw(Tools.f(t0) + "..." + Tools.f(t1), np);
          V.pop();
        }

      }
    }

    //  V.pop(2);
  }

//  public static void render(IPlaneCurve curve, Color color, int stroke,
//      int markType) {
//    // final boolean db = false;
//
//    V.pushStroke(stroke, Globals.STRK_NORMAL);
//    V.pushColor(color, MyColor.cBLUE);
//
//    //    double step = Math.max(1, ((c - a) * (c - a) / c) * .3);
//    //
//    //    if (db)
//    //      Streams.out.println(" step=" + step);
//    //
//    // plot each visible segment
//
//    for (int seg = 0; seg < nComponents(); seg++) {
//      curve.render(cStart(seg), cEnd(seg));
//
//      //curve.render(cStart(seg), cEnd(seg));
//      //      double t0 = this.cStart(seg), t1 = this.cEnd(seg);
//      //      t0 = MyMath.clamp(t0, -500.0, 500.0);
//      //      t1 = MyMath.clamp(t1, -500.0, 500.0);
//      //
//      //      double s0, s1;
//      //
//      //      s0 = t0;
//      //      s1 = t1;
//      //
//      //      if (s0 > s1) {
//      //        double tmp = s0;
//      //        s0 = s1;
//      //        s1 = tmp;
//      //      }
//      //      FPoint2 p0 = this.pointAt(s0), p1 = pointAt(s1);
//      //      if (db)
//      //        Streams.out.println(" p0=" + p0 + ", p1=" + p1);
//      //      {
//      //        int count = 0;
//      //        boolean first = true;
//      //        for (double t = t0;; t += step, count++) {
//      //          boolean last = (t >= t1);
//      //          if (last)
//      //            t = t1;
//      //          pointAt(t, p1);
//      //
//      //          if (!p1.isValid()) {
//      //            if (last) {
//      //              break;
//      //            }
//      //            continue;
//      //          }
//      //
//      //          if (db) {
//      //            System.out.println(" calcPt " + Tools.f(t) + " = " + p1.x + ","
//      //                + p1.y);
//      //          }
//      //          if (!first)
//      //            V.drawLine(p0, p1);
//      //          if (last)
//      //            break;
//      //          p0.setLocation(p1);
//      //          first = false;
//      //        }
//    }
//    V.pop(2);
//
//  }

}
