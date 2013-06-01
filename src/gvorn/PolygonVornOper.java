package gvorn;

import base.*;
import testbed.*;
import static gvorn.Main.FULL;

public class PolygonVornOper implements TestBedOperation, Globals {
  /*! .enum  .private 1800
      _ _ ISOLATECELL _ _ _ PLOTFARTHEST
      _ _ _  retain clear plotstd plotguar showstd showbisector showleading showleading2
      _ _
  */

  private static final int ISOLATECELL = 1802;//!
  private static final int PLOTFARTHEST = 1806;//!
  private static final int RETAIN = 1810;//!
  private static final int CLEAR = 1811;//!
  private static final int PLOTSTD = 1812;//!
  private static final int PLOTGUAR = 1813;//!
  private static final int SHOWSTD = 1814;//!
  private static final int SHOWBISECTOR = 1815;//!
  private static final int SHOWLEADING = 1816;//!
  private static final int SHOWLEADING2 = 1817;//!
  /*!*/
 // private static final boolean FULL = false;

  public void addControls() {
    C.sOpenTab("Polygons");
    C.sStaticText("Plots guaranteed Voronoi cell for uncertain polygons.  "
        + "Hover over guaranteed cell with shift key pressed to show induced disc.");

    {
      {
        C.sOpen();
        C.sIntSpinner(ISOLATECELL, "cell", "Select cell to display", 0, 20, 0,
            1);
        //        C.sHide();
        //        C.sIntSpinner(RES, "res", "sampling resolution", 0, 100, 75, 1);
        if (FULL) {
          C.sNewColumn();
          C.sButton(RETAIN, "Retain", "Add highlighted sample to render list");
          C.sButton(CLEAR, "Clear", "Clear render list");
        }
        C.sClose();

        C.sOpen();
        C.sCheckBox(SHOWSTD, "std edges", "plot portions of standard edges",
            false);
        C.sCheckBox(PLOTFARTHEST, "f.p. Vorn",
            "Plot farthest-point Voronoi diagram", false);
        C.sNewColumn();
        C.sCheckBox(PLOTSTD, "std disc", "highlight std disc", false);
        C.sCheckBox(PLOTGUAR, "guar disc", "highlight guaranteed disc", true);
        if (FULL) {
          C.sCheckBox(SHOWBISECTOR, "bisector",
              "highlight associated bisectors", false);
          C.sCheckBox(SHOWLEADING, "bisector (t+1)",
              "highlight bisector associated with next src component", false);
          C.sCheckBox(SHOWLEADING2, "bisector (r+1)",
              "highlight bisector associated with next opp component", false);
        }
        C.sClose();
      }
      //      {
      //        C.sOpen();
      //        C.sCheckBox(SCALED, "scaled", "plot scaled version of cell's polygon",
      //            false);
      //        C.sNewColumn();
      //        C.sIntSlider(SCALEFACTOR, "factor", "scale factor", 0, 100, 50, 1);
      //        C.sClose();
      //      }
    }

    C.sCloseTab();
  }
  public static PolygonVornOper singleton = new PolygonVornOper();

  private PolygonVornOper() {
  }

  public void processAction(TBAction a) {
    if (a.code == TBAction.HOVER && a.shiftPressed()) {
      Sample hPrev = hSample;
      hSample = cellInfo.nearestSample(a.loc);
      if (hPrev != hSample)
        V.repaint();
    }
    if (a.code == TBAction.CTRLVALUE) {
      switch (a.ctrlId) {
      case RETAIN:
        if (hSample != null)
          retained.add(hSample);
        break;
      case CLEAR:
        retained.clear();
        break;
      }

    }

    //    if (a.code == TBAction.CTRLVALUE && a.ctrlId == STEP)
    //      C.seti(ISOLATECIRC, C.vi(ISOLATECIRC) + 12);
  }

  private void processHash() {
    DArray poly = Main.getPolygons();
    EdPolygon[] pl = (EdPolygon[]) poly.toArray(EdPolygon.class);
    int icell = C.vi(ISOLATECELL);
    icell = Math.min(icell, pl.length - 1);

    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < pl.length; i++) {
      EdPolygon p = pl[i];
      sb.append(p.getHash());
    }
    sb.append("i:" + icell);

    int sr = 50;// C.vi(RES);
    sb.append("s:" + sr);

    //    boolean sc = C.vb(SCALED);
    //    sb.append(" sc:" + sc);
    //    double scaleFactor = C.vi(SCALEFACTOR) / 100.0;
    //    sb.append(scaleFactor);

    String newHash = sb.toString();
    if (!newHash.equals(hash)) {
      hash = newHash;
      polys = pl;
      hSample = null;
      isoCell = icell;
      furthestVDiag = null;
      cellInfo = null;
      sampRes = sr;
      //      scaled = sc;
      // double scaleFactor = C.vi(SCALEFACTOR);

      //      scaledPoly = null;
      //      if (scaled) {
      //        scaledPoly = buildScaledVersion(polys[isoCell], scaleFactor);
      //      }
      retained.clear();
    }
  }

  public void runAlgorithm() {

    processHash();
    if (polys.length == 0)
      return;

    cellInfo = new CellInfo();

    //    if (isoCell < 0) {
    //      for (int i = 0; i < polys.length; i++) {
    //        isoCell = i;
    //        sampleCell(polys[i],cellInfo);
    //      }
    //      isoCell = -1;
    //    }
    //    else
    sampleCell(polys[isoCell], cellInfo);

    //    if (scaledPoly != null)
    //      sampleCell(scaledPoly, cellInfo);

    if (C.vb(PLOTFARTHEST)) {
      DArray poly = Main.getPolygons();
      furthestVDiag = new VornGraph(
          buildFarthestPointVoronoiDiagram(((EdPolygon) poly.get(isoCell))
              .getPoints()));
    }
  }

  //  private EdPolygon buildScaledVersion(EdPolygon p, double t) {
  //    FRect r = p.getBounds();
  //    FPoint2 origin = r.midPoint();
  //    EdPolygon sp = new EdPolygon();
  //    // double t = scaleFactor / 100.0;
  //
  //    for (int i = 0; i < p.nPoints(); i++) {
  //      FPoint2 pt = p.getPoint(i);
  //      pt = FPoint2.interpolate(origin, pt, t);
  //      sp.addPoint(pt);
  //    }
  //    return sp;
  //  }

  public void paintView() {
    Editor.render();
    //    T.render(scaledPoly, MyColor.cPURPLE);

    if (C.vb(PLOTFARTHEST) && furthestVDiag != null)
      furthestVDiag.plot(MyColor.get(MyColor.PURPLE, .7), -1,
          DiscVornOper.hlVert(), DiscVornOper.labelEdges());

    plotVorn();

    if (hSample != null) {
      hSample.highlight();
    }

    for (int i = 0; i < retained.size(); i++) {
      ((Sample) retained.get(i)).highlight();
    }

  }
  private static final double MIN_RADIUS_DIFF = 1e-3;

  /**
   * Sample data for cell
   * @param V
   * @param isoCell
   */
  private void sampleCell(EdPolygon cellPoly, CellInfo ci) {

    double step = MyMath.radians(((sampRes + 10) * .5) / 100);

    // EdPolygon cellPoly = polys[isoCell];

    // test each vertex of polygon for the cell
    for (int cvi = 0; cvi < cellPoly.nPoints(); cvi++) {

      FPoint2 vpt = cellPoly.getPoint(cvi);

      double angMin, angMax;
      {
        FPoint2 v0 = cellPoly.sideStart(cvi - 1);
        FPoint2 v1 = cellPoly.sideEnd(cvi);

        angMin = MyMath
            .normalizeAngle(MyMath.polarAngle(v0, vpt) + Math.PI / 2);
        angMax = MyMath
            .normalizeAngle(MyMath.polarAngle(vpt, v1) + Math.PI / 2);
      }

      // if not enough space between this range of angles, skip this vertex
      if (MyMath.normalizeAngle(angMax - angMin) < step)
        continue;

      // sample from min...max angle for this vertex
      double angSample = angMin + step * .5;

      int tComponent = 1 + cvi;

      while (true) {

        Sample best = null;
        {
          double dx = Math.cos(angSample);
          double dy = Math.sin(angSample);

          // set tangent point back a bit so polygon pt is always inside circle
          final double SETBACK = .0001;
          FPoint2 tpt = new FPoint2(vpt.x - dx * SETBACK, vpt.y - dy * SETBACK);

          // find smallest radius that contains polygon and does not intersect
          // other polygons

          double r0 = .01;
          double r1 = 3000;
          double r = (r0 + r1) * .5;

          // keep track of last neighboring polygon and part that
          // was intersected
          EdPolygon nbrPoly = null;

          // use binary search to find largest disc that
          // is supported by the vertex, and is at that angle
          int rComponent = 0;

          while (true) {
            FPoint2 cpt = new FPoint2(tpt.x + r * dx, tpt.y + r * dy);
            boolean inCirc = cellPoly.withinCircle(cpt, r, null); // extrPt);

            // determine if other polygons intersect this circle
            boolean isectsNeighbor = false;
            for (int neighborPolyIndex = 0; neighborPolyIndex < polys.length; neighborPolyIndex++) {
              if (neighborPolyIndex == isoCell)
                continue;
              EdPolygon neighborPoly = polys[neighborPolyIndex];

              int testIsectPart = neighborPoly.intersectsCircle(cpt, r);
              if (testIsectPart != 0) {
                isectsNeighbor = true;
                nbrPoly = neighborPoly;
                rComponent = testIsectPart;
                break;
              }
            }

            if (!inCirc && isectsNeighbor)
              break;

            if (!inCirc) {
              r0 = r;
            } else if (isectsNeighbor) {
              r1 = r;
            } else {
              if (nbrPoly != null)
                best = new Sample(cellPoly, nbrPoly, cpt, r, tComponent,
                    rComponent);
              r0 = r;
            }

            if (r1 - r0 < MIN_RADIUS_DIFF) {
              break;
            }
            r = (r0 + r1) * .5;
          }
        }

        // add sample, if one was found
        if (best != null) {
          // determine standard voronoi disc, by shrinking toward other poly
          // until it excludes main poly

          // find smallest radius that contains polygon and does not intersect
          // other polygons

          double t0 = 0;
          double t1 = 1.0;
          double t = .5;

          FPoint2 rPt = best.nbrPoly.closestBoundaryPointTo(best.guarPt);

          // use binary search to find largest disc that excludes main polygon
          while (t1 - t0 > .005) {

            t = (t0 + t1) * .5;
            FPoint2 cpt = FPoint2.interpolate(rPt, best.guarPt, t);
            double tRad = cpt.distance(rPt);

            FPoint2 critPt = new FPoint2();
            int comp = cellPoly.closestBoundaryPointTo(cpt, critPt);
            if (critPt.distance(cpt) >= tRad) {
              best.setStdDisc(cpt, comp);

              t0 = t;
            } else
              t1 = t;
          }
          ci.add(best);
        }

        angSample += step;
        double diff = MyMath.normalizeAngle(angSample - angMax);
        if (diff > 0)
          break;
      }
    }
  }

  private static EdPoint[] sitesFor(FPoint2[] pts) {
    EdPoint[] a = new EdPoint[pts.length];
    for (int i = 0; i < pts.length; i++) {
      a[i] = new EdPoint(pts[i]);
    }
    return a;
  }

  /**
   * Build a furthest-point standard Voronoi diagram from a set of FPoint2's
   * @param points
   * @return
   */
  private static Hyperbola[] buildFarthestPointVoronoiDiagram(FPoint2[] pts) {
    return VornUtil.build(sitesFor(pts), pts.length - 1, PointBisector.S);
  }

  private void plotVorn() {
    V.pushColor(MyColor.cDARKGRAY);
    if (cellInfo != null)
      cellInfo.render();
    V.pop();

  }

  private static class Sample {

    /**
     * Constructor
     * @param nbrPolygon
     * @param guarPt origin of disc
     * @param guarRadius
     */
    public Sample(EdPolygon srcPolygon, EdPolygon nbrPolygon, FPoint2 guarPt,
        double guarRadius, int sourceComponent, int nbrComponent) {
      this.srcPoly = srcPolygon;
      this.nbrPoly = nbrPolygon;
      this.guarPt = new FPoint2(guarPt);
      this.radius = guarRadius;
      this.guarSourceComponent = sourceComponent;
      this.nbrComponent = nbrComponent;
    }

    public void setStdDisc(FPoint2 stdPt, int stdSourceComponent) {
      this.stdPoint = stdPt;
      this.stdSourceComponent = stdSourceComponent;
    }

    private static int nextInsideComponent(int c, EdPolygon poly, int adj) {
      if (c > 0) {
        c--;
        c = (c + adj) % poly.nPoints();
        c++;
      }
      return c;
    }

    private static int nextOutsideComponent(int c, EdPolygon poly, int adj) {
      final boolean db = false;

      if (db)
        Streams.out.println("nextOutsideComponent c=" + c + " #npts="
            + poly.nPoints() + " adj=" + adj + "\n"
            + Editor.write(poly, new StringBuilder()));

      if (c > 0) {
        c--;
        c = (c + adj) % poly.nPoints();
        c = -c - 1;
      } else {
        c = -c - 1;
        c = (c + adj + 1) % poly.nPoints();
        c = c + 1;
      }
      if (db)
        Streams.out.println(" returning " + c);
      return c;
    }

    public void highlight() {
      if (C.vb(PLOTGUAR)) {
        V.pushColor(MyColor.cRED);
        V.mark(guarPt, MARK_X);
        V.pushStroke(STRK_RUBBERBAND);
        V.drawCircle(guarPt, radius);

        if (FULL) {
          V.pushColor(MyColor.cBLUE);
          if (C.vb(SHOWLEADING)) {
            plotBisector(srcPoly,
                nextInsideComponent(guarSourceComponent, srcPoly, 1),
                this.nbrPoly, this.nbrComponent, stdPoint, true);
          }

          if (C.vb(SHOWLEADING2)) {

            plotBisector(srcPoly, guarSourceComponent, this.nbrPoly,
                nextOutsideComponent(nbrComponent, nbrPoly, -1), stdPoint, true);

          }
          V.pop();
        }
        plotBisector(srcPoly, guarSourceComponent, this.nbrPoly,
            this.nbrComponent, guarPt, FULL && C.vb(SHOWBISECTOR));

        V.pop(2);
      }

      if (C.vb(PLOTSTD)) {
        if (stdPoint != null) {
          V.pushColor(MyColor.cPURPLE);
          V.mark(stdPoint, MARK_X);
          V.pushStroke(STRK_RUBBERBAND);
          double auxRad = radius - stdPoint.distance(guarPt);
          V.drawCircle(stdPoint, auxRad);
          //          if (C.vb(SHOWBISECTOR)) 
          {
            plotBisector(srcPoly, this.stdSourceComponent, this.nbrPoly,
                this.nbrComponent, stdPoint, FULL && C.vb(SHOWBISECTOR));

          }

          V.pop(2);
        }
      }
    }
    private void plotBisector(EdPolygon poly1, int comp1, EdPolygon poly2,
        int comp2, FPoint2 target, boolean plotCurve) {

      V.pushStroke(STRK_THIN);
      IPlaneCurve curve = null;
      if (comp1 > 0) {
        FPoint2 p1 = poly1.getPoint(comp1 - 1);
        if (comp2 > 0) {
          FPoint2 p2 = poly2.getPoint(comp2 - 1);
          double theta = MyMath.polarAngle(p1, p2) + Math.PI / 2;
          curve = new LineCurve(target, theta);
          V.mark(p1, MARK_DISC);
          V.mark(p2, MARK_DISC);
        } else {
          int edge2 = (-comp2) - 1;

          FPoint2 p2a = poly2.sideStart(edge2);
          FPoint2 p2b = poly2.sideEnd(edge2);
          V.pushStroke(STRK_THICK);
          V.drawLine(p2a, p2b);
          V.pop();
          curve = new Parabola(p2a, p2b, p1);
          V.mark(p1, MARK_DISC);
        }
      } else {
        int edge1 = (-comp1) - 1;

        FPoint2 p1a = poly1.sideStart(edge1);
        FPoint2 p1b = poly1.sideEnd(edge1);
        V.pushStroke(STRK_THICK);
        V.drawLine(p1a, p1b);
        V.pop();

        if (comp2 > 0) {
          FPoint2 p2 = poly2.getPoint(comp2 - 1);
          V.mark(p2, MARK_DISC);
          curve = new Parabola(p1a, p1b, p2);
        } else {

          int edge2 = (-comp2) - 1;

          FPoint2 p2a = poly2.sideStart(edge2);
          FPoint2 p2b = poly2.sideEnd(edge2);
          V.pushStroke(STRK_THICK);
          V.drawLine(p2a, p2b);
          V.pop();

          double theta;
          FPoint2 iPt = MyMath.linesIntersection(p1a, p1b, p2a, p2b, null);
          if (iPt == null)
            theta = MyMath.polarAngle(p1a, p1b);
          else {
            //            double th1 = MyMath.polarAngle(p1a, p1b);
            //            double th2 = MyMath.polarAngle(p2a, p2b);
            theta = MyMath.polarAngle(iPt, target);
          }
          curve = new LineCurve(target, theta);
        }

      }
      if (curve != null && plotCurve) {
        double t = curve.parameterFor(target);
        curve.render(t - 50, t + 50);
      }
      V.pop();
    }

    public void render(Sample next) {
      final double MAX_RAD = 1000;
      if (next != null) {
        if (radius < MAX_RAD || next.radius < MAX_RAD) {
          V.drawLine(guarPt, next.guarPt);

          if (nbrComponent != next.nbrComponent
              || guarSourceComponent != next.guarSourceComponent)
            V.mark(guarPt, MARK_DISC, .65);

          if (C.vb(SHOWSTD)) {
            if (nbrComponent != next.nbrComponent
                || stdSourceComponent != next.stdSourceComponent)
              V.mark(stdPoint, MARK_DISC, .65);
            if (stdPoint != null && next.stdPoint != null
                && nbrPoly == next.nbrPoly)
              V.drawLine(stdPoint, next.stdPoint);
          }
        }
      }
    }

    double radius;

    // FPoint2 tMain;
    FPoint2 guarPt;
    FPoint2 stdPoint;

    EdPolygon srcPoly;
    EdPolygon nbrPoly;

    int guarSourceComponent, nbrComponent;
    int stdSourceComponent;
  }

  private static class CellInfo {

    /**
     * Find sample corresponding to location
     * @param loc location
     * @return sample, or null if not found
     */
    public Sample nearestSample(FPoint2 loc) {
      //      System.out.println("finding nearest sample to "+loc);
      Sample ret = null;
      double bestDist = 5.0;

      for (int i = 0; i < nSamp(); i++) {
        Sample s = sample(i);
        double dist = s.guarPt.distanceSq(loc);
        if (s.stdPoint != null)
          dist = Math.min(dist, s.stdPoint.distanceSq(loc));

        if (dist < bestDist) {
          bestDist = dist;
          ret = s;
        }
      }
      //System.out.println("nearest sample="+ret);
      return ret;
    }

    private DArray samples = new DArray();
    public int nSamp() {
      return samples.size();
    }
    public Sample sample(int i) {
      return (Sample) samples.get(i);
    }
    public CellInfo() {
    }
    public void add(Sample b) {
      samples.add(b);
    }

    public void render() {
      if (samples.isEmpty())
        return;

      Sample sPrev = (Sample) samples.last();
      for (int i = 0; i < samples.size(); i++) {
        Sample s = (Sample) samples.get(i);
        //        if (sPrev != null)
        sPrev.render(s);
        sPrev = s;
        //        Sample sNext = (Sample) samples.get((i + 1) % samples.size());
        //        s.render(sNext);
      }
    }

  }

  private VornGraph furthestVDiag;
  private CellInfo cellInfo;
  private int isoCell;
  private EdPolygon[] polys;
  private String hash;
  private int sampRes;
  private Sample hSample;
  private DArray retained = new DArray();
  //  private boolean scaled;
  // private int scaleFactor;
  // private EdPolygon scaledPoly;
}
