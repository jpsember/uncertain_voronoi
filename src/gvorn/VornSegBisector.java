package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class VornSegBisector extends VornBisector {

  /**
   * Subcurve that makes up segment bisector
   */
  private static class Piece extends VornBisector {

    public IPlaneCurve curve(int curveIndex) {
      return curve;
    }

    /**
     * Do prepare stage 1 for piece.
     * @param startPt point representing first endpoint encountered on 
     * larger bisector
     */
    public void prepare1(FPoint2 startPt) {
      double t = this.parameterFor(startPt);

      t0 = cStart(0);
      t1 = cEnd(0);
      if (Math.abs(t0 - t1) < .5) {
        t1 = t0;
      }
      boolean flipped = Math.abs(t - t0) > Math.abs(t - t1);
      if (flipped) {
        double tmp = t0;
        t0 = t1;
        t1 = tmp;
      }
      //  inPath = true;
      tRange = Math.abs(t1 - t0);
    }

    public void prepare2(double s0, double s1) {
      // determine group->individual t conversion
      double t0 = cStart(0);
      double t1 = cEnd(0);
      a = (t0 - t1) / (s0 - s1);
      b = t0 - a * s0;
    }

    /**
     * Determine if curve has been incorporated into global path
     * @return
     */
    public boolean inPath() {
      return t0 != t1;
    }

    public double tToExternal(double t) {
      return (t - b) / a;
    }

    public double sToInternal(double s) {
      return a * s + b;
    }

    public double tRange() {
      return tRange;
    }
    public String toString() {
      StringBuilder sb = new StringBuilder();
      sb.append("Piece[");
      sb.append(pointAt(this.cStart(0)));
      sb.append(" ... ");
      sb.append(pointAt(this.cEnd(0)));
      sb.append("]");
      return sb.toString();
    }

    /**
     * Construct linear bisector between two endpoints
     * @param sa
     * @param sb
     * @param pa
     * @param pb
     */
    public Piece(VornSegSite sa, VornSegSite sb, int pa, int pb) {
      super(sa, sb);

      FPoint2 pta = sa.pt(pa);
      FPoint2 ptb = sb.pt(pb);

      FPoint2 mid = FPoint2.midPoint(pta, ptb);
      double theta = MyMath.polarAngle(pta, ptb) + Math.PI / 2;

      curve = new LineCurve(mid, theta);
    }

    /**
     * Construct segment bisector that is parabola
     * @param sa first site (represents line of parabola)
     * @param sb second site 
     * @param pb endpoint of second site that is focus
     */
    public Piece(VornSegSite sa, VornSegSite sb, int pb) {
      super(sa, sb);

      FPoint2 a0 = sa.pt(0);
      FPoint2 a1 = sa.pt(1);
      FPoint2 ptb = sb.pt(pb);

      curve = new Parabola(a0, a1, ptb);
    }

    /**
     * Construct linear bisector
     * @param sa
     * @param sb sites
     * @param mPt point on line
     * @param mTheta angle of line
     */
    public Piece(VornSegSite sa, VornSegSite sb, FPoint2 mPt, double mTheta) {
      super(sa, sb);

      curve = new LineCurve(mPt, mTheta);
    }

    @Override
    public double parameterFor(FPoint2 pt) {
      return curve.parameterFor(pt);
    }

    @Override
    public FPoint2 pointAt(double t) {
      return curve.pt(t);
    }

    @Override
    public void render(Color c, int stroke, int markType) {
      render(curve, c, stroke, markType);
    }

    public void markInvalid() {
      invalid = true;
    }
    public boolean invalid() {
      return invalid;
    }

    // local start/end parameter values
    private double t0, t1;
    // constants for converting global <=> local parameters
    private double a, b;

    //    private boolean inPath;
    private double tRange;

    private IPlaneCurve curve;

    private boolean invalid;

    public FPoint2 oppositePoint(FPoint2 pt) {
      double t = parameterFor(pt);

      // don't use t values, since maybe not yet defined
      double c0 = cStart(0);
      double c1 = cEnd(0);

      double tOpp = (Math.abs(t - c0) < Math.abs(t - c1)) ? c1 : c0;
      return pointAt(tOpp);
    }

  }

  public double parameterFor(FPoint2 pt) {
    // return the parameter that is the closest match of the pieces
    double tBest = -1;
    boolean bestFound = false;
    double bestDist = 0;

    for (int i = 0; i < nPieces(); i++) {
      Piece p = piece(i);
      double s = p.parameterFor(pt);
      FPoint2 tPt = p.pointAt(s);
      double dist = tPt.distance(pt);
      if (!bestFound || dist < bestDist) {
        bestFound = true;
        bestDist = dist;

        // convert t from internal to externa
        tBest = p.tToExternal(s);
      }
    }
    return tBest;
  }

  @Override
  public FPoint2 pointAt(double s) {

    // find piece that contains s
    int i = 0;
    Piece p = piece(i);
    double pEnd = sStart + p.tRange();
    while (++i < nPieces() && s > pEnd) {
      p = piece(i);
      pEnd += p.tRange();
    }
    return p.pointAt(p.sToInternal(s));
  }

  @Override
  public void render(Color c, int stroke, int markType) {
    V.pushStroke(stroke, Globals.STRK_RUBBERBAND);
    V.pushColor(c, MyColor.cDARKGRAY);
    for (int i = 0; i < this.nComponents(); i++) {
      double s0 = this.cStart(i);
      double s1 = this.cEnd(i);

      //double tp0 = 0;
      for (int j = 0; j < nPieces(); j++) {
        Piece p = piece(j);

        double t0 = p.sToInternal(s0);
        double t1 = p.sToInternal(s1);
        //        double tp1 = p.tRange() + tp0;

        double tMin = p.cStart(0);
        double tMax = p.cEnd(0);

        t0 = MyMath.clamp(t0, tMin, tMax);
        t1 = MyMath.clamp(t1, tMin, tMax);
        if (t1 > t0) {
          p.curve.render(t0, t1);
        }
      }
    }
    V.pop(2);
  }

  public VornSegBisector(VornSegSite sa, VornSegSite sb) {
    super(sa, sb);

    final boolean db = false;

    if (db && T.update())
      T.msg("VornSegBisector constructor " + T.show(sa, MyColor.cRED)
          + T.show(sb, MyColor.cPURPLE));

    DArray r = new DArray();

    for (int i = 0; i < 2; i++) {
      VornSegSite s = (i == 0) ? sa : sb;
      VornSegSite sopp = (i == 0) ? sb : sa;

      for (int j = 0; j < 2; j++) {

        {
          r.add(new Piece(sa, sb, i, j));
          if (db && T.update())
            T.msg(" pt/pt edge" + T.show(r.last()) + T.show(sa.pt(i))
                + T.show(sb.pt(j)));
        }

        {
          FPoint2 pt = sopp.pt(j);
          if (MyMath.ptDistanceToLine(pt, s.pt(0), s.pt(1), null) < 1e-3)
            continue;
          r.add(new Piece(s, sopp, j));
          if (db && T.update())
            T.msg(" pt/seg edge" + T.show(r.last()) + T.show(s, MyColor.cRED)
                + T.show(sopp.pt(j)));
        }
      }

    }

    // add bisector for medial axis
    {
      FPoint2 a0 = sa.pt(0), a1 = sa.pt(1), b0 = sb.pt(0), b1 = sb.pt(1);

      double th0 = MyMath.polarAngle(a0, a1);
      double th1 = MyMath.polarAngle(b0, b1);
      if (th0 < 0)
        th0 += Math.PI * 2;
      if (th1 < 0)
        th1 += Math.PI * 2;

      // point on medial axis is intersection point, if it exists
      FPoint2 mPt = MyMath.linesIntersection(a0, a1, b0, b1, null);

      if (mPt == null) {
        mPt = FPoint2.midPoint(a0, b0);
      }
      for (int pass = 0; pass < 2; pass++) {
        double m0 = (th0 + th1) / 2;
        if (pass == 1)
          m0 += Math.PI / 2;
        if (db && T.update())
          T.msg("medial axis, th0=" + Tools.fa(th0) + " th1=" + Tools.fa(th1)
              + " m0=" + Tools.fa(m0));
        r.add(new Piece(sa, sb, mPt, m0));
        if (db && T.update())
          T.msg(" seg/seg edge" + T.show(r.last()) + T.show(sa, MyColor.cRED)
              + T.show(sb, MyColor.cPURPLE));
      }
    }

    // construct arrangement of bisectors, clip dishonest edges
    r = clipArrangement(sa, sb, r);

    // stitch together the pieces into a path
    stitch(r);
  }

  /**
   * For stitching together segments:  graph vertices
   */
  private static class PathVertex implements Renderable {
    final boolean db = true;
    /**
     * Constructor
     * @param pt location of vertex
     */
    public PathVertex(FPoint2 pt) {
      this.pt = pt;
    }
    /**
     * Get # edges incident with vertex
     * @return
     */
    public int degree() {
      return edges.size();
    }

    /**
     * Add edge to vertex
     * @param edge
     */
    public void addEdge(Piece edge) {
      edges.add(edge);
      if (db && T.update())
        T.msg("adding edge, degree now " + degree()
            + T.show(edges, MyColor.cDARKGRAY, Globals.STRK_RUBBERBAND, -1)
            + T.show(edge, MyColor.cRED) + T.show(pt));
    }

    /**
     * Get edge
     * @param i edge number
     * @return edge
     */
    public Piece getEdge(int i) {
      return (Piece) edges.get(i);
    }
    public void render(Color c, int stroke, int markType) {
      V.mark(pt, Globals.MARK_X);
      V.pushColor(MyColor.cLIGHTGRAY);
      V.pushStroke(Globals.STRK_RUBBERBAND);
      for (int i = 0; i < degree(); i++) {
        Piece p = getEdge(i);
        V.drawLine(p.pointAt(p.cStart(0)), p.pointAt(p.cEnd(0)));
      }
      V.pop(2);
    }
    public String toString() {
      StringBuilder sb = new StringBuilder();
      sb.append("EndPt ");
      sb.append(pt);
      sb.append(" degree=" + degree());
      return sb.toString();
    }
    private FPoint2 pt;
    private DArray edges = new DArray();
  }

  private PathVertex findEndPt(FPoint2 pt) {
    PathVertex ret = null;
    for (int k = 0; k < epts(); k++) {
      PathVertex ep = ept(k);
      if (ep.pt.distance(pt) < .5) {
        ret = ep;
        break;
      }
    }
    return ret;
  }
  private int epts() {
    return endPts.size();
  }
  private PathVertex ept(int i) {
    return (PathVertex) endPts.get(i);
  }

  private Piece piece(int i) {
    return (Piece) pieces.get(i);
  }

  //  private static double[] steps = { .5,
  //  //.2, .7,
  //  };

  private static DArray clipArrangement(VornSegSite sa, VornSegSite sb,
      DArray bl) {

    final boolean db = false;

    if (db && T.update())
      T.msg("clipArrangement" + T.showAll(bl, MyColor.cRED, -1, -1));

    DArray[] crossLists = new DArray[bl.size()];
    for (int i = 0; i < crossLists.length; i++)
      crossLists[i] = new DArray();

    for (int i = 0; i < bl.size(); i++) {
      VornBisector bi = (VornBisector) bl.get(i);
      for (int j = i + 1; j < bl.size(); j++) {
        VornBisector bj = (VornBisector) bl.get(j);

        DArray ipts = VUtil2.intersect(bi, bj, false);

        //        if (false && db && T.update()) {
        //          T.msg("intersections between " + T.show(bi) + T.show(bj)
        //              + T.show(VUtil2.calcPoints(bi.curve, ipts)));
        //        }

        for (int k = 0; k < ipts.size(); k += 2) {
          crossLists[i].add(ipts.get(k));
          crossLists[j].add(ipts.get(k + 1));
        }
      }
    }

    DArray ret = new DArray();
    for (int i = 0; i < crossLists.length; i++) {
      Piece bi = (Piece) bl.get(i);
      DArray cl = crossLists[i];

      cl.addDouble(VornBisector.TMIN);
      cl.addDouble(VornBisector.TMAX);

      cl.sort(DArray.COMPARE_DOUBLES);

      if (db && T.update()) {
        DArray ipts = new DArray();
        for (int k = 1; k < cl.size() - 1; k++)
          ipts.add(bi.pointAt(cl.getDouble(k)));
        T.msg("cross points" + T.show(ipts) + T.show(bi) + cl.toString(true));
      }

      double t1 = cl.getDouble(0);
      double t0 = t1;

      for (int j = 1; j < cl.size(); j++, t0 = t1) {
        double t = cl.getDouble(j);

        //        // if essentially the same as the previous point, treat
        //        // as if this point doesn't exist
        //        if (false && t - t0 < 1e-2) {
        //          t1 = t0;
        //          continue;
        //        }
        t1 = t;

        //   passLoop: for (int pass = 0; pass < steps.length; pass++)
        {
          double frac = .5; //steps[pass];
          t = (t0 * frac + t1 * (1 - frac));

          if (false) {
            double f0 = Math.max(Math.min(t1 - 20, 100), t0);
            double f1 = Math.min(Math.max(t0 + 20, -100), t1);
            t = (f0 + f1) * .5;
          }

          FPoint2 tPt = bi.pointAt(t);

          if (false && db && T.update())
            T.msg("testing if midpoint is on bisector; t0=" + t0 + " t1=" + t1
                + T.show(bi.pointAt(t0), null, -1, Globals.MARK_X)
                + T.show(bi.pointAt(t1), null, -1, Globals.MARK_X)
                + T.show(tPt));

          //   double tol = (pass == 0) ? 1e-4 : 1e-2;
          if (!bi.onBisector(tPt))
            bi.clipOut(t0, t1);
          //          break passLoop;
        }
      }

      if (db && T.update())
        T.msg("after clipping bisector, nComp=" + bi.nComponents() + T.show(bi)
            + "\n" + bi.componentsStr());
      if (bi.nComponents() > 0)
        ret.add(bi);
    }
    if (T.update())
      T.msg("clipped Arrangement"
          + T.showAll(ret, MyColor.cRED, -1, Globals.MARK_X));

    return ret;
  }
  @Override
  public int nPieces() {
    return pieces.size();
  }

  public IPlaneCurve curve(int i) {
    return piece(i).curve;
  }

  private String toString(DArray r) {
    StringBuilder sb = new StringBuilder("[\n");
    for (int i = 0; i < r.size(); i++) {
      Piece p = (Piece) r.get(i);
      if (p.invalid())
        continue;
      sb.append("piece#" + i + ": " + p + "\n");
    }
    return sb.toString();
  }

  private String dumpVerts() {
    StringBuilder sb = new StringBuilder("vertices:[\n");
    for (int i = 0; i < endPts.size(); i++) {
      PathVertex ep = ept(i);
      sb.append(ep.pt + " deg:" + ep.degree() + "\n");
    }
    sb.append("]");
    return sb.toString();
  }

  /**
   * Stitch pieces together into a path
   * @param r array of Pieces
   */
  private void stitch(DArray r) {
    final boolean db = true;

    if (db && T.update())
      T.msg("stitch pieces" + T.show(r));

    int nEdges = 0;
    pieces = new DArray();
    endPts = new DArray();

    // build graph from pieces; filter out pieces that are too short
    for (int i = 0; i < r.size(); i++) {
      Piece p = (Piece) r.get(i);
      if (db && T.update())
        T.msg("stitch piece #" + i + T.show(p, MyColor.cRED));

      // if piece is too short, skip it
      if (Math.abs(p.cStart(0) - p.cEnd(0)) < 1e-1) {
        if (db && T.update())
          T.msg("piece length is too short, skipping");
        p.markInvalid();
        continue;
      }
      nEdges++;

      for (int j = 0; j < 2; j++) {
        FPoint2 pt = p.pointAt(j == 0 ? p.cStart(0) : p.cEnd(0));
        PathVertex ptVertex = findEndPt(pt);
        if (ptVertex == null) {
          ptVertex = new PathVertex(pt);
          endPts.add(ptVertex);
          if (db && T.update())
            T.msg("adding new vertex " + T.show(pt) + pt);
        }
        ptVertex.addEdge(p);
      }
    }

    // filter out edges that connect degree 1 vertices; they are
    // probably there due to precision problems or something
    if (db && T.update())
      T.msg("filtering out singleton edges\n" + toString(r) + dumpVerts());
    for (int i = 0; i < endPts.size(); i++) {
      PathVertex ep = ept(i);
      if (ep.degree() != 1)
        continue;
      Piece pc = ep.getEdge(0);
      if (db && T.update())
        T.msg("degree 1 vertex=" + ep.pt + T.show(ep.pt) + " piece=" + pc
            + T.show(pc, MyColor.cRED));

      // probably shouldn't be invalid at this point
      if (pc.invalid())
        continue;

      FPoint2 oppPt = pc.oppositePoint(ep.pt);
      PathVertex oppVert = findEndPt(oppPt);
      if (oppVert == null) {
        if (db && T.update())
          T.msg("can't find opp vert! " + oppPt + T.show(oppPt));
      }

      if (oppVert.degree() == 1) {
        if (db && T.update())
          T.msg("marking singleton edge as invalid" + T.show(pc));
        pc.markInvalid();
        nEdges--;
      }
    }

    if (db && T.update()) {
      T.msg(toString(r) + dumpVerts());
    }

    // start path with an endpt of degree 1
    if (db && T.update())
      T.msg("looking for degree 1 vertex to start path");
    PathVertex startPt = null;
    for (int i = 0; i < epts(); i++) {
      PathVertex endPt = ept(i);
      if (db && T.update())
        T.msg("degree=" + T.show(endPt.pt) + endPt.degree() + "\n" + endPt);
      if (endPt.degree() == 1) {
        startPt = endPt;
        break;
      }
    }
    if (startPt == null)
      throw new FPError("can't find degree 1 endpoint");

    double tTotal = 0;

    FPoint2 currPoint = startPt.pt;
    // PathVertex prevPt = startPt;

    while (pieces.size() < nEdges) {

      //   Piece prev = currPiece;
      //   FPoint2 nextPt = currPiece.pointAt(currPiece.t1);
      //      if (pt.distance(prevPt.pt) < .1) {
      //        pt = prev.pointAt(prev.cStart(0));
      //      }
      if (db && T.update())
        T.msg("finding next edge from " + currPoint + T.show(currPoint)
            + T.show(currPoint));

      //      FPoint2 pt = prev
      //          .pointAt((!prev.flipped) ? prev.cEnd(0) : prev.cStart(0));
      PathVertex found = findEndPt(currPoint);
      if (found == null) {
        if (db && T.update())
          T.msg("can't find endpoint " + currPoint + T.show(currPoint));
        //        if (true) {
        //          Tools.warn("ignoring suspected problem");
        //          break;
        //        }
        throw new FPError("can't find endpoint");
      }
      // prevPt = found;

      Piece nextPiece = null;
      for (int i = 0; i < found.degree(); i++) {
        Piece p2 = found.getEdge(i);
        if (db && T.update())
          T.msg("looking for edge incident with " + currPoint
              + " not already in path" + T.show(p2, MyColor.cRED) + "\n"
              + endPts.toString(true));
        if (!p2.inPath()) {
          nextPiece = p2;
          break;
        }
      }
      if (nextPiece == null) {
        //        if (true) {
        //          Tools.warn("ignoring susp. prob.");
        //          break;
        //        }
        throw new FPError("can't find next edge");
      }

      nextPiece.prepare1(currPoint);
      ////      Piece currPiece = startPt.getEdge(0);
      //      FPoint2 currPoint = startPt.pt;
      //      //    Piece p = startPt.getEdge(0);
      //      currPiece.prepare1(currPoint);
      tTotal += nextPiece.tRange();
      pieces.add(nextPiece);

      currPoint = nextPiece.pointAt(nextPiece.t1);

      //      nextPiece.prepare1(nextPt);

      // p.prepare1(pt);

      //      tTotal += nextPiece.tRange(); //Math.abs(p.cStart(0) - p.cEnd(0)); //))p.tRange();
      //      pieces.add(nextPiece);

      //      currPiece = nextPiece;
    }
    // throw out now that not needed
    endPts = null;

    // determine group -> individual t mapping
    sStart = -tTotal / 2;
    double s0 = sStart;

    for (int i = 0; i < nPieces(); i++) {
      Piece p = piece(i);
      double s1 = s0 + p.tRange();
      p.prepare2(s0, s1);
      s0 = s1;
    }
    this.setClipRange(sStart, s0);
  }
  private double sStart;
  // array of pieces
  private DArray pieces;
  // temporary use by stitching
  private DArray endPts;

}
