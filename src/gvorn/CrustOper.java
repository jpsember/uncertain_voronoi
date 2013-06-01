package gvorn;

import base.*;
import testbed.*;

public class CrustOper implements   TestBedOperation, Globals {
  /*! .enum  .private 2300
       crust
       delaunay
       betaskel
       beta
       oldgdel
       plotr3
  */

    private static final int CRUST            = 2300;//!
    private static final int DELAUNAY         = 2301;//!
    private static final int BETASKEL         = 2302;//!
    private static final int BETA             = 2303;//!
    private static final int OLDGDEL          = 2304;//!
    private static final int PLOTR3           = 2305;//!
/*!*/

  public void addControls() {
    C.sOpenTab("Crust");
    //   C.sStaticText("Crust and Beta-skeleton");
    {
      {
        C.sOpen();
        C.sCheckBox(CRUST, "Plot crust", null, true);
        C.sNewColumn();
        C
            .sCheckBox(PLOTR3, "Plot |S|>=3 regions",
                "Plots regions of subset diagram for sets with >= 3 elements",
                true);
        C.sClose();
      }
      {
        C.sOpen();
        C.sCheckBox(DELAUNAY, "Plot Delaunay",
            "Plots guaranteed Delaunay edges", false);

        C.sNewColumn();
        C
            .sCheckBox(OLDGDEL, "Old method",
                "Plot edges constructed using old method, for test purposes",
                false);
        C.sClose();
      }
      if (false) {
        C.sOpen();
        C.sCheckBox(BETASKEL, "Plot Beta skeleton", null, false);
        C.sNewColumn();
        C.sDoubleSpinner(BETA, "Beta:", "Value of beta", 1.0, 20, 1.7, .1);
        C.sClose();
      }
    }
    C.sCloseTab();
  }

  public static CrustOper singleton = new CrustOper();

  private CrustOper() {
  }

//  public void paintView() {
//    // build various items in runAlgorithm() method, so we can use
//    // algorithm tracing features if we want to
//    T.runAlgorithm(this);
//
//    // use plotResults() method to display results of algorithm
//  }

  public void runAlgorithm() {


    // clear any objects constructed during last paint cycle
    crust = null;
    delaunay = null;
    hTrace = null;
    dTrace = null;
    r3Regions = null;
    pointSites = null;


    
    //    try {
    if (C.vb(CRUST)) {
      DArray a = build(Main.getDiscs());
      crust = a.getDArray(0);
      r3Regions = (VornGraph) a.get(1);
      pointSites = (VornGraph) a.get(2);
    }
    if (C.vb(DELAUNAY)) {
      delaunay = Spine.buildGuaranteedDelaunay(Main.getDiscs());
    }
    //    } catch (FPError e) {
    //      Tools.warn("caught FPError");
    //    }

  }

  private static String ts(Hyperbola[] ha) {
    StringBuilder sb = new StringBuilder();
    sb.append(" " + ha.length + ":\n");
    for (int i = 0; i < ha.length; i++) {
      Hyperbola h = ha[i];
      sb.append(" " + h);
      sb.append(" " + h.calcPoint(h.minParameter()) + " .. "
          + h.calcPoint(h.maxParameter()));
      sb.append(" " + h.visString());
      sb.append("\n");
    }
    return sb.toString();
  }

  private static boolean overlap(EdDisc d, Hyperbola h, FPoint2 pRet) {

    final boolean db = true;
    if (db && T.update())
      T.msg("overlap, " + d + ", h=" + h);
    double t = h.closestPointTo(d.getOrigin());
    if (db && T.update())
      T.msg(" closestPt = " + t);
    FPoint2 pt = null;
    if (!h.containsPoint(t)) {
      double t0 = h.minParameter();
      pt = h.calcPoint(t0);
      FPoint2 p2 = h.calcPoint(h.maxParameter());
      if (pt.distance(d.getOrigin()) > p2.distance(d.getOrigin()))
        pt = p2;
    } else
      pt = h.calcPoint(t);
    if (pRet != null)
      pRet.setLocation(pt);
    return pt.distance(d.getOrigin()) < d.getRadius();
  }

  public static void plot(EdDisc[] di, DArray edges) {
    // vp vp = TestBed.view();
    V.pushColor(MyColor.get(MyColor.GRAY, .3));
    for (int i = 0; i < edges.size(); i += 2) {
      int p0 = edges.getInt(i), p1 = edges.getInt(i + 1);

      FPoint2 e0 = di[p0].getOrigin();
      FPoint2 e1 = di[p1].getOrigin();
      V.drawLine(e0, e1);
    }
    V.popColor();
  }

  public void paintView() {
    // vp vp = TestBed.view();

    final boolean db = false;

    if (db)
      Streams.out.println("CrustOper paintview");

    EdDisc[] di = Main.getDiscs();

    Editor.render();

    if (delaunay != null) {
      if (C.vb(OLDGDEL)) {
        V.pushColor(MyColor.get(MyColor.GRAY, .3));
        Hyperbola[] spines = Spine.construct(di);
        V.pushStroke(STRK_THICK);
        for (int i = 0; i < spines.length; i++) {
          Hyperbola h = spines[i];
          EdDisc ca = (EdDisc) h.getData(0);
          EdDisc cb = (EdDisc) h.getData(1);
          FPoint2[] ep = EdDisc.lineBetween(ca, cb, false);
          if (ep != null)
            V.drawLine(ep[0], ep[1]);
        }
        V.popStroke();
        V.popColor();
      }

      V.pushColor(MyColor.get(MyColor.GREEN, .4));
      for (int i = 0; i < delaunay.size(); i += 2) {
        EdDisc d0 = di[delaunay.getInt(i)];
        EdDisc d1 = di[delaunay.getInt(i + 1)];
        V.drawLine(d0.getOrigin(), d1.getOrigin());
      }
      V.popColor();
    }

    if (crust != null) {
      for (int i = 0; i < crust.size(); i += 2) {
        int p0 = crust.getInt(i), p1 = crust.getInt(i + 1);

        FPoint2 e0 = di[p0].getOrigin();
        FPoint2 e1 = di[p1].getOrigin();
        V.drawLine(e0, e1);
      }
    }

    Tools.warn("add Trace abiliity");
    T.render(hTrace);
    T.render(dTrace);
    V.pushColor(MyColor.get(MyColor.BROWN, .9));

    if (C.vb(PLOTR3)) {

      if (r3Regions != null) {
        r3Regions.render(null,-1,-1);
        //)render(true, false, false);
      }
      if (pointSites != null) {
        pointSites.render(null, -1, -1); //render(false, true, false);
      }
    }
    V.popColor();

  }

  public void processAction(TBAction a) {
    if (a.code == TBAction.CTRLVALUE) {
      switch (a.ctrlId) {
      }
    }
  }

  /**
   * Build crust for discs
   * 
   * @return DArray, with [0] : array of ints, indexes i,j of guaranteed edges
   *   [1] graph of S3 regions
   *   [2] graph of Voronoi diagram of point sites (for rendering Voronoi vertices)
   */
  public static DArray build(EdDisc [] d) {

    final boolean db = false;

    DArray crust = new DArray();
    VornGraph r3Regions;
    VornGraph pointSites;

    Hyperbola[] h3 = null;

    // get list of region >= 3 edges
    {
      Hyperbola[] h1 = VornUtil.buildSubset(d, true, SubsetDiscBisector.S);
      h1 = Hyperbola.singlify(h1);
      h1 = VornGraph.addCrossings(h1);

      // filter out edges that occur in guaranteed voronoi diagram
      Hyperbola[] hGuarVorn = VornUtil.build(d, 1, GuaranteedDiscBisector.S);
      hGuarVorn = Hyperbola.singlify(hGuarVorn);
      if (db && T.update()) {
        T.msg("subset, with crossings:\n" + ts(h1));
      }
      if (db && T.update()) {
        T.msg("guaranteed:\n" + ts(hGuarVorn));
      }

      DArray f = new DArray();
      outer: for (int i = 0; i < h1.length; i++) {
        Hyperbola h = h1[i];
        if (h.isEmpty())
          continue;
        boolean db2 = false;
        for (int j = 0; j < hGuarVorn.length; j++) {
          Hyperbola b = hGuarVorn[j];
          if (b.isEmpty())
            continue;

          if (db2 && T.update()) {

            T.msg("seeing if " + h + " matches " + b + " range="
                + h.minParameter() + ".." + h.maxParameter() + " range="
                + b.minParameter() + ".." + b.maxParameter());
          }

          if (h.getData(0) == b.getData(0) && h.getData(1) == b.getData(1)
              && b.containsRange(h.minParameter(), h.maxParameter())) {
            if (db2)
              Streams.out.println(h + " matches " + b + ", leaving out");
            continue outer;
          }
        }
        f.add(h);
      }
      h3 = (Hyperbola[]) f.toArray(Hyperbola.class);
      if (db && T.update())
        T.msg("filtered:\n" + ts(h3));
      r3Regions = new VornGraph(h3);
    }

    pointSites = VornUtil.vornForPointSites(d);
    DArray nodeIds = pointSites.getNodeList();

    // crust is a subset of guaranteed Delaunay edges
    DArray ce = Spine.buildGuaranteedDelaunay(d);

    DArray ptSites = new DArray();

    // add disc centerpoints, in case they have zero radius and
    // don't generate any S3 regions
    for (int i = 0; i < d.length; i++)
      ptSites.add(d[i].getOrigin());

    for (int jj = 0; jj < nodeIds.size(); jj++) {
      int id = nodeIds.getInt(jj);

      if (pointSites.isInf(id))
        continue;

      FPoint2 vl = pointSites.getVertLoc(id);
      ptSites.add(vl);
    }

    nextEdge: for (int i = 0; i < ce.size(); i += 2) {

      int ii = ce.getInt(i);
      int ij = ce.getInt(i + 1);
      EdDisc di = d[ii];
      EdDisc dj = d[ij];

      if (di.getRadius() < dj.getRadius()) {
        EdDisc t = di;
        di = dj;
        dj = t;
      }

      boolean contains = EdDisc.contains(di, dj);

      Hyperbola h = null;
      final double MAX = 20000;
      final double EPS = 1e-4;
      double t0 = -MAX + 1;
      double t1 = MAX - 1;

      if (db && T.update())
        T.msg("testing edge " + di.getLabel() + " to " + dj.getLabel());

      if (!contains)
        h = DiscUtil.supportingHyperbola(di, dj);

      // Test against point Voronoi sites
      for (int jj = 0; jj < ptSites.size(); jj++) {
        if (jj == ii || jj == ij)
          continue;

        FPoint2 vl = ptSites.getFPoint2(jj);
        if (db && T.update())
          T.msg(" testing against point site " + vl);
        if (contains) {
          if (EdDisc.encloses(di, vl))
            continue nextEdge;
          continue;
        }

        // Perform two passes; pass 0 finds max support t, pass 1 finds min
        for (int pass = 0; pass < 2; pass++) {

          boolean penetrated = true;

          double tLow = -MAX + 1, tHigh = MAX - 1;
          double tMid = 0;

          while (tLow < tHigh - EPS) {
            tMid = (tLow + tHigh) * .5;

            EdDisc ds = Spine.getSupportedDisc(h.calcPoint(tMid), di);

            boolean pen = EdDisc.encloses(ds, vl);
            int sign = Spine.penetratingDiscPosition(ds, di, dj, vl);
            if (db && T.update())
              T.msg(" supported disc " + ds + "\n penetrates=" + pen + " sign="
                  + sign);
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

      // Test against R3 regions
      for (int k = 0; k < h3.length; k++) {
        Hyperbola ht = h3[k];

        if (db && T.update())
          T.msg("testing supported disc penetrated by " //
              + ht + " " + ht.visString());

        if (contains) {
          if (overlap(di, ht, null))
            continue nextEdge;
          continue;
        }

        // Perform two passes; pass 0 finds max support t, pass 1 finds min

        for (int pass = 0; pass < 2; pass++) {

          boolean penetrated = true;

          double tLow = t0, tHigh = t1;
          double tMid = 0;

          while (tLow < tHigh - EPS) {
            tMid = (tLow + tHigh) * .5;

            EdDisc ds = Spine.getSupportedDisc(h.calcPoint(tMid), di);

            FPoint2 opt = new FPoint2();
            boolean pen = overlap(ds, ht, opt);

            int sign = Spine.penetratingDiscPosition(ds, di, dj, opt);

            if (db && T.update())
              T.msg("  pass=" + pass + " tm=" + Tools.f(tMid) + " pen=" + pen
                  + " tlow=" + Tools.f(tLow) + " thigh=" + Tools.f(tHigh)
                  + "\n sign=" + sign);

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

      if (db && T.update())
        T.msg("adding edge between " + di.getLabel() + " and " + dj.getLabel());
      crust.addInt(ii);
      crust.addInt(ij);
    }
    return DArray.build(crust, r3Regions, pointSites);
  }

  private DArray delaunay;
  private DArray crust;
  private VornGraph r3Regions;
  private VornGraph pointSites;

  private static Hyperbola hTrace;
  private static EdDisc dTrace;
}
