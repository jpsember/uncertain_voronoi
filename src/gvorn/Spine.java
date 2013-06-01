package gvorn;

import java.util.*;
import testbed.*;
import base.*;

public class Spine {

  //  public static void main(String[] args) {
  //    
  //    
  //    int A = 1, B = 2;
  //    
  //    FPoint2 p2 = new FPoint2(1,6);
  //    
  //    for (double t = 2.45; t < 2.46; t += .0001) {
  //      double xt = Math.sqrt(A+B*t*t);
  //      double yt = t;
  //      FPoint2 p1 = new FPoint2(xt,yt);
  //      
  //      double dist = p1.distance(p2);
  //      Streams.out.println("t="+Tools.f(t*100)+" dist="+dist);
  //    }
  //  }
  /**
   * Generate guaranteed Delaunay edges for a set of discs
   * @param d : discs
   * @return array of 2n integers, each two a pair of indexes into discs
   */
  public static DArray buildGuaranteedDelaunay(EdDisc[] d) {
    final boolean db = false;

    DArray ret = new DArray();

    for (int i = 0; i < d.length; i++) {
      nextEdge: for (int j = i + 1; j < d.length; j++) {
        EdDisc di = d[i];
        EdDisc dj = d[j];
        if (di.getRadius() < dj.getRadius()) {
          EdDisc t = di;
          di = dj;
          dj = t;
        }

        if (EdDisc.contains(di, dj)) {
          for (int k = 0; k < d.length; k++) {
            if (k == i || k == j)
              continue;
            EdDisc dk = d[k];
            if (EdDisc.overlap(di, dk))
              continue nextEdge;
          }
        } else {

          Hyperbola h = null;
          final double MAX = 20000;
          h = DiscUtil.supportingHyperbola(di, dj);
          double t0 = -MAX + 1;
          double t1 = MAX - 1;
          //          if (db)
          //            hTrace = h;

          if (db && T.update())
            T.msg("testing edge " + di.getLabel() + " to " + dj.getLabel());

          for (int k = 0; k < d.length; k++) {
            if (k == i || k == j)
              continue;
            EdDisc dk = d[k];
            if (db && T.update())
              T.msg("testing supported disc penetrated by " + dk.getLabel()
                  + "\n t0=" + Tools.f(t0) + " t1=" + Tools.f(t1));
            final double EPS = 1e-5;

            // determine range of t that doesn't cause dk to penetrate supported disc

            // Perform two passes; pass 0 finds max support t, pass 1 finds min

            for (int pass = 0; pass < 2; pass++) {

              boolean penetrated = true;

              double tLow = -MAX + 1, tHigh = MAX - 1;
              double tMid = 0;

              while (tLow < tHigh - EPS) {
                tMid = (tLow + tHigh) * .5;

                EdDisc ds = getSupportedDisc(h.calcPoint(tMid), di);
                boolean pen = EdDisc.overlap(ds, dk);

                //                boolean pen = suppDiscPenetrated(h, di, dj, tMid, dk, db);
                int sign = penetratingDiscPosition(ds, di, dj, EdDisc
                    .nearestPointToOrigin(ds, dk));
                //penetratingDiscPosition(ds, di, dj, dk);

                if (db && T.update())
                  T.msg("  pass=" + pass + " tm=" + Tools.f(tMid) + " pen="
                      + pen + " tlow=" + Tools.f(tLow) + " thigh="
                      + Tools.f(tHigh) + "\n sign=" + sign);

                if (pen) {
                  if (sign > 0)
                    tLow = tMid;
                  else
                    tHigh = tMid;
                } else {
                  penetrated = false;
                  if (pass == 0)
                    tLow = tMid;
                  else
                    tHigh = tMid;
                }
              }
              if (penetrated)
                continue nextEdge;

              if (db && T.update())
                T.msg(" pass=" + pass + " tm=" + Tools.f(tMid) + " t0="
                    + Tools.f(t0) + " t1=" + Tools.f(t1));

              // clip bounds on t; abort if it becomes empty
              if (pass == 0) {
                if (tMid < t0)
                  continue nextEdge;
                t1 = Math.min(t1, tMid);
              } else {
                if (tMid > t1)
                  continue nextEdge;
                t0 = Math.max(t0, tMid);
              }
            }
          }
        }
        if (db && T.update())
          T.msg("adding edge between " + di.getLabel() + " and "
              + dj.getLabel());
        ret.addInt(i);
        ret.addInt(j);
      }
    }
    return ret;
  }

  /**
   * Construct disc that has another inside-tangent to it
   * @param p center of disc
   * @param da disc that is to be inside-tangent to this one
   * @return
   */
  public static EdDisc getSupportedDisc(FPoint2 p, EdDisc da) {
    double pr = p.distance(da.getOrigin()) + da.getRadius();
    return new EdDisc(p, pr);
  }

  /**
   * Determine where a disc is with respect to a supported disc
   * @param dS : supported disc
   * @param da, dB : discs inside tangent to dS; da is disc on right
   * @param dk : disc to test
   * @return -1 if support disc should move lower on hyperbolic arm to
   *   move farther away from dS; 1 if it should move higher
   */
  public static int penetratingDiscPosition(EdDisc dS, EdDisc da, EdDisc dB,
      FPoint2 q) {
    final boolean db = false;

    FPoint2 p = dS.getOrigin();
    double pr = dS.getRadius();

    double aAng = MyMath.polarAngle(p, da.getOrigin());
    double bAng = MyMath.polarAngle(p, dB.getOrigin());

    int ret = 0;

    // If test disc is enclosed by supported disc, test which side of the
    // line between the support points the test discs's center lies on

    if (EdDisc.encloses(dS, q)
    //  FPoint2.distance(p, dk.getOrigin()) < pr 
    ) {
      FPoint2 ta = MyMath.ptOnCircle(p, aAng, pr);
      FPoint2 tb = MyMath.ptOnCircle(p, bAng, pr);
      ret = (int) Math.signum(MyMath.sideOfLine(ta, tb, q));
    } else {
      // Otherwise, calculate the angle of a point of intersection of the
      // boundaries (the ray from the supported discs's center through the
      // test disc's center will contain such a point); see where this lies
      // relative to the tangent points on the boundary

      double kAng = MyMath.polarAngle(p, q);

      if (angBetween(aAng, bAng, kAng))
        ret = -1;
      else
        ret = 1;

      if (db && T.update())
        T.msg("discPos a=" + da.getLabel() + " b=" + dB.getLabel() + " q=" + q
            + "\n aAng=" + Tools.fa(aAng) + " bAng=" + Tools.fa(bAng)
            + " kAng=" + Tools.fa(kAng) + " ret=" + ret);
    }
    return ret;
  }

  /**
   * Determine where a disc is with respect to a supported disc
   * @param dS : supported disc
   * @param da, dB : discs inside tangent to dS; da is disc on right
   * @param dk : disc to test
   * @return -1 if support disc should move lower on hyperbolic arm to
   *   move farther away from dS; 1 if it should move higher
   */
  public static int penetratingDiscPosition(EdDisc dS, EdDisc da, EdDisc dB,
      EdDisc dk) {
    final boolean db = false;

    FPoint2 p = dS.getOrigin();
    double pr = dS.getRadius();

    double aAng = MyMath.polarAngle(p, da.getOrigin());
    double bAng = MyMath.polarAngle(p, dB.getOrigin());

    int ret = 0;

    // If test disc is enclosed by supported disc, test which side of the
    // line between the support points the test discs's center lies on

    if (EdDisc.encloses(dS, dk)
    //  FPoint2.distance(p, dk.getOrigin()) < pr 
    ) {
      FPoint2 ta = MyMath.ptOnCircle(p, aAng, pr);
      FPoint2 tb = MyMath.ptOnCircle(p, bAng, pr);
      ret = (int) Math.signum(MyMath.sideOfLine(ta, tb, dk.getOrigin()));
    } else {
      // Otherwise, calculate the angle of a point of intersection of the
      // boundaries (the ray from the supported discs's center through the
      // test disc's center will contain such a point); see where this lies
      // relative to the tangent points on the boundary

      double kAng = MyMath.polarAngle(p, dk.getOrigin());

      if (angBetween(aAng, bAng, kAng))
        ret = -1;
      else
        ret = 1;

      if (db && T.update())
        T.msg("discPos a=" + da.getLabel() + " b=" + dB.getLabel() + " dk="
            + dk.getLabel() + "\n aAng=" + Tools.fa(aAng) + " bAng="
            + Tools.fa(bAng) + " kAng=" + Tools.fa(kAng) + " ret=" + ret);
    }
    return ret;
  }

  private static boolean angBetween(double a, double b, double t) {
    boolean ret = false;
    if (b >= a) {
      ret = t >= a && t <= b;
    } else {
      ret = t >= a || t <= b;
    }
    return ret;
  }

  private static Hyperbola buildSpine(Graph g, int iId, int jId) {
    EdDisc iCircle = (EdDisc) g.nodeData(iId);
    EdDisc jCircle = (EdDisc) g.nodeData(jId);

    Hyperbola hSpine = null;

    // if one disc contains the other, use bisector (containment method #1)
    if (EdDisc.contains(iCircle, jCircle) || EdDisc.contains(jCircle, iCircle)) {
      hSpine = new Hyperbola(iCircle.getOrigin(), jCircle.getOrigin());
    } else {
      double distAB = FPoint2
          .distance(iCircle.getOrigin(), jCircle.getOrigin());
      double distPA = (distAB + jCircle.getRadius() - iCircle.getRadius()) * .5;
      FPoint2 ptOnArm = FPoint2.interpolate(iCircle.getOrigin(), jCircle
          .getOrigin(), distPA / distAB);
      hSpine = new Hyperbola(iCircle.getOrigin(), jCircle.getOrigin(), ptOnArm);
    }

    hSpine.setLabel(iCircle.getLabel() + "|" + jCircle.getLabel());
    hSpine.setData(iCircle, jCircle);

    double tMin = hSpine.minParameter();
    double tMax = hSpine.maxParameter();

    // clip against every <i,z>, <j,z> arm, where z is a neighbor node
    // in dual graph
    for (int pass = 0; pass < 2; pass++) {
      int id = pass == 0 ? iId : jId;
      for (int ni = 0; ni < g.nCount(id); ni++) {
        int id2 = g.neighbor(id, ni);
        if (id2 == iId || id2 == jId)
          continue;
        EdDisc zCircle = (EdDisc) g.nodeData(id2);

        for (int pass2 = 0; pass2 < 2; pass2++) {
          EdDisc cSrc = (pass2 == 0) ? iCircle : jCircle;
          if (EdDisc.overlap(cSrc, zCircle)) {
            hSpine.clipAll();
            break;
          }

          Hyperbola hArm = GuaranteedDiscBisector.S.getBisector(cSrc, zCircle); //VornUtil.constructWCBisector(cSrc, zCircle);

          FPoint2[] iPts = null;
          iPts = (FPoint2[]) Hyperbola.findIntersections(hSpine, hArm, null)
              .toArray(FPoint2.class);

          switch (iPts.length) {
          case 0:
            {
              FPoint2 spinePt = hSpine.calcPoint(0);
              if (hArm.testPoint(spinePt) < 0)
                hSpine.clipAll();
            }
            break;
          case 1:
            {
              FPoint2 ipt = iPts[0];

              double ti = hSpine.calcParameter(ipt);
              // move a little further along, and compare resulting point
              FPoint2 i2 = hSpine.calcPoint(ti + .5);
              if (hArm.testPoint(i2) < 0) {
                hSpine.clip(ti, Hyperbola.CLIP_MAX);
              } else {
                hSpine.clip(Hyperbola.CLIP_MIN, ti);
              }
            }
            break;

          case 2:
            {
              FPoint2 ipt0 = iPts[0];
              FPoint2 ipt1 = iPts[1];
              Hyperbola ha = hSpine;
              Hyperbola hb = hArm;
              double t0 = ha.calcParameter(ipt0);
              double t1 = ha.calcParameter(ipt1);
              FPoint2 i2 = ha.calcPoint((t0 + t1) * .5);
              if (hb.testPoint(i2) < 0) {
                ha.clip(t0, t1);
              } else {
                ha.clipAllBut(t0, t1);
              }
            }
            break;
          }

          if (!hSpine.isEmpty()) {
            double new0 = hSpine.minParameter();

            if (new0 > tMin) {
              hSpine.setData(Hyperbola.MINCLIPPED, hArm.toString()); //zCircle);
              tMin = new0;
            }

            double new1 = hSpine.maxParameter();

            if (new1 < tMax) {
              hSpine.setData(Hyperbola.MAXCLIPPED, hArm.toString());
              tMax = new1;
            }
          }
        }
      }
    }
    return hSpine;
  }
  /**
   * Construct spines for a set of discs, using dual of additively-weighted
   * Voronoi diagram to restrict number of sites to examine
   * @param edCircles : DArray of EdCircles
   * @return array of Hyperbolas
   */
  public static Hyperbola[] construct(EdDisc[] edCircles) {
    final boolean db = false;

    // build additively-weighted V. diagram for discs
    Hyperbola[] vdiag = VornUtil.build(edCircles, 1, StandardDiscBisector.S);

    //      VornUtil.build(edCircles);
    Graph g = VornUtil.buildDualGraph(vdiag);

    DArray edgeList = new DArray();

    //    DArray cl  =new DArray();
    /*
     * The Fortune sweep line algorithm can detect if a disc is contained
     * within another.  When the line encounters the weighted disc center, 
     * if the actual disc center is behind the wavefront, then the disc is
     * contained by the disc associated with the parabola which contains that
     * point of the wavefront.
     * 
     * We augment the graph to contain hyperbolic arcs that are simply 
     * the bisectors of the two site centers, and clip accordingly.
     */
    {
      // get mapping of discs to graph nodes
      Map m = new HashMap();
      for (Iterator it = g.getNodeList().iterator(); it.hasNext();) {
        int k = ((Integer) it.next()).intValue();
        EdDisc c = (EdDisc) g.nodeData(k);
        m.put(c, new Integer(k));
      }

      // detect contained discs using slow but sure method (we could use
      // Fortune here)

      for (int i = 0; i < edCircles.length; i++)
        edCircles[i].clearFlags(DiscUtil.DISC_CONTAINED);

      for (int i = 0; i < edCircles.length; i++) {
        EdDisc ca = edCircles[i];
        if (DiscUtil.contained(ca))
          continue;

        int ccount = 0;
        EdDisc cb = null;

        for (int j = 0; j < edCircles.length; j++) {
          if (j == i)
            continue;
          if (db)
            Streams.out.println("seeing if " + ca.getLabel() + " contains "
                + edCircles[j].getLabel());

          if (DiscUtil.contained(edCircles[j])) {
            if (db)
              Streams.out.println("already contained");
            continue;
          }

          if (EdDisc.contains(ca, edCircles[j])) {
            ccount++;
            cb = edCircles[j];
            cb.addFlags(DiscUtil.DISC_CONTAINED); //)cb.setContained(ca);
          }
        }

        if (ccount == 1) {
          if (db)
            Streams.out.println("ccount=1 for " + ca.getLabel() + " and "
                + cb.getLabel());

          Integer caId = (Integer) m.get(ca);
          if (caId == null) {
            caId = new Integer(g.newNode(ca));
            m.put(ca, caId);
          }

          Integer ccId = (Integer) m.get(cb);
          if (ccId == null) {
            ccId = new Integer(g.newNode(cb));
            m.put(cb, ccId);
          }

          //            Streams.out.println("containment spine for "+ca+" and "+cb);
          Hyperbola hSpine = new Hyperbola(ca.getOrigin(), cb.getOrigin());

          hSpine.setLabel(ca.getLabel() + "|" + cb.getLabel());
          hSpine.setData(ca, cb);
          if (db)
            Streams.out.println("adding edge between " + caId + " and " + ccId);

          //          cl.add(hSpine);
          g.addEdgesBetween(caId.intValue(), ccId.intValue(), hSpine, hSpine);
        }
      }
    }

    if (db)
      Streams.out.println("iterating through graph:\n" + g);

    for (Iterator qi = g.getNodeList().iterator(); qi.hasNext();) {

      int iId = ((Integer) qi.next()).intValue();
      EdDisc a = (EdDisc) g.nodeData(iId);
      if (db)
        Streams.out.println(" node " + iId + " is " + g.nodeData(iId));

      for (int jn = 0; jn < g.nCount(iId); jn++) {
        int jId = g.neighbor(iId, jn);
        if (jId < iId)
          continue;
        EdDisc c = (EdDisc) g.nodeData(jId);

        if (db)
          Streams.out.println("  neighbor " + c);

        if (DiscUtil.contained(c) && DiscUtil.contained(a)) {
          //)c.contained() && a.contained()) {
          if (db)
            Streams.out.println(" contained");
          continue;
        }

        Hyperbola hSpine = null;
        hSpine = buildSpine(g, iId, jId);

        if (!hSpine.isEmpty())
          edgeList.add(hSpine);
      }

    }

    return (Hyperbola[]) edgeList.toArray(Hyperbola.class);
  }
  //  /**
  //   * Construct spines for a set of discs
  //   * 
  //   * @param edCircles : DArray of EdCircles
  //   * @param vdiag : Voronoi diagram of disc centerpoints
  //   * @return array of Hyperbolas
  //   * @deprecated
  //   */
  //  public static Hyperbola[] construct0(EdCircle[] edCircles, Hyperbola[] vdiag) {
  //
  //    DArray edgeList = new DArray();
  //    outer: for (int qq = 0; qq < vdiag.length; qq++) {
  //      Hyperbola ve = vdiag[qq];
  //      EdCircle cj = (EdCircle) ve.getData(0);
  //      EdCircle ck = (EdCircle) ve.getData(1);
  //
  //      if (cj.getLabel().compareTo(ck.getLabel()) > 0) {
  //        continue;
  //      }
  //
  //      Hyperbola hSpine = null;
  //
  //      // If one disc contains the other, construct special type of 
  //      // hyperbola, since there is no spine.
  //      for (int p2 = 0; p2 < 2; p2++) {
  //        EdCircle ca = p2 == 0 ? cj : ck;
  //        EdCircle cb = p2 == 0 ? ck : cj;
  //
  //        if (EdCircle.contains(ca, cb)) {
  //          for (int z = 0; z < edCircles.length; z++) {
  //            EdCircle cz = edCircles[z];
  //            if (cz == cj || cz == ck)
  //              continue;
  //            if (EdCircle.overlap(cz, ca)) {
  //              continue outer;
  //            }
  //          }
  //          FPoint2 or = ca.getOrigin();
  //
  //          hSpine = new Hyperbola(or, new FPoint2(or.x + 10, or.y), 3);
  //          hSpine.setLabel(ca.getLabel() + "|" + cb.getLabel());
  //          hSpine.setData(ca, cb);
  //          hSpine.setUserFlag(true);
  //
  //          edgeList.add(hSpine);
  //          continue outer;
  //        }
  //      }
  //
  //      // construct spine(j,k)
  //      double distAB = FPoint2.distance(cj.getOrigin(), ck.getOrigin());
  //      double distPA = (distAB + ck.radius() - cj.radius()) * .5;
  //      FPoint2 ptOnArm = FPoint2.interpolate(cj.getOrigin(), ck.getOrigin(), distPA
  //          / distAB);
  //      hSpine = new Hyperbola(cj.getOrigin(), ck.getOrigin(), ptOnArm);
  //      hSpine.setLabel(cj.getLabel() + "|" + ck.getLabel());
  //      hSpine.setData(cj, ck);
  //
  //      // clip against every <j,z>, <k,z> arm
  //      for (int z = 0; z < edCircles.length; z++) {
  //        EdCircle cz = edCircles[z];
  //        if (cz == cj || cz == ck)
  //          continue;
  //
  //        for (int pass = 0; pass < 2; pass++) {
  //          EdCircle cSrc = (pass == 0) ? cj : ck;
  //          if (EdCircle.overlap(cSrc, cz)) {
  //            hSpine.clipAll();
  //            break;
  //          }
  //
  //          Hyperbola hArm = DiscVornOper.constructWCBisector(cSrc, cz);
  //
  //          FPoint2[] iPts = (FPoint2[]) Hyperbola.findIntersections(hSpine,
  //              hArm, null).toArray(FPoint2.class);
  //
  //          switch (iPts.length) {
  //          case 0:
  //            {
  //              FPoint2 spinePt = hSpine.calcPoint(0);
  //              if (hArm.testPoint(spinePt) < 0)
  //                hSpine.clipAll();
  //            }
  //            break;
  //          case 1:
  //            {
  //              FPoint2 ipt = iPts[0];
  //
  //              double ti = hSpine.calcParameter(ipt);
  //              // move a little further along, and compare resulting point
  //              FPoint2 i2 = hSpine.calcPoint(ti + .5);
  //              if (hArm.testPoint(i2) < 0) {
  //                hSpine.clip(ti, Hyperbola.CLIP_MAX);
  //              } else {
  //                hSpine.clip(Hyperbola.CLIP_MIN, ti);
  //              }
  //            }
  //            break;
  //
  //          case 2:
  //            {
  //              FPoint2 ipt0 = iPts[0];
  //              FPoint2 ipt1 = iPts[1];
  //              Hyperbola ha = hSpine;
  //              Hyperbola hb = hArm;
  //              double t0 = ha.calcParameter(ipt0);
  //              double t1 = ha.calcParameter(ipt1);
  //              FPoint2 i2 = ha.calcPoint((t0 + t1) * .5);
  //              if (hb.testPoint(i2) < 0) {
  //                ha.clip(t0, t1);
  //              } else {
  //                ha.clipAllBut(t0, t1);
  //              }
  //            }
  //            break;
  //          }
  //          if (hSpine.isEmpty())
  //            break;
  //        }
  //      }
  //
  //      if (!hSpine.isEmpty()) {
  //        edgeList.add(hSpine);
  //      }
  //    }
  //
  //    //    edgeList.sort(hComparator);
  //    return (Hyperbola[]) edgeList.toArray(Hyperbola.class);
  //  }

  //  private static Comparator hComparator = new Comparator() {
  //
  //    public int compare(Object a, Object b) {
  //      Hyperbola ha = (Hyperbola) a;
  //      Hyperbola hb = (Hyperbola) b;
  //      return ha.toString().compareTo(hb.toString());
  //    }
  //  };

}
