package gvorn;

import base.*;
import testbed.*;
import static gvorn.Main.FULL;

public class DiscVornOper implements TestBedOperation, Globals {
  /*! .enum  .private 1500
  plotguar minmaxposs minmax bendduals
  hledges printedges cliptocircle agstore agframe subsets noclip dualcircles 
  plotstd
  effdiscs labeledges spines gdel  test perturb _ _ _
  stdk guark _ _ _ _ minmaxk
  plotposs psite _ subsetcolor plotcrust plotstored store plotlast
  */

  private static final int PLOTGUAR = 1500;//!
  private static final int MINMAXPOSS = 1501;//!
  private static final int MINMAX = 1502;//!
  private static final int BENDDUALS = 1503;//!
  private static final int HLEDGES = 1504;//!
  private static final int PRINTEDGES = 1505;//!
  private static final int CLIPTOCIRCLE = 1506;//!
  private static final int AGSTORE = 1507;//!
  private static final int AGFRAME = 1508;//!
  private static final int SUBSETS = 1509;//!
  private static final int NOCLIP = 1510;//!
  private static final int DUALCIRCLES = 1511;//!
  private static final int PLOTSTD = 1512;//!
  private static final int EFFDISCS = 1513;//!
  private static final int LABELEDGES = 1514;//!
  private static final int SPINES = 1515;//!
  private static final int GDEL = 1516;//!
  private static final int TEST = 1517;//!
  private static final int PERTURB = 1518;//!
  private static final int STDK = 1522;//!
  private static final int GUARK = 1523;//!
  private static final int PLOTPOSS = 1529;//!
  private static final int PSITE = 1530;//!
  private static final int SUBSETCOLOR = 1532;//!
  private static final int PLOTCRUST = 1533;//!
  private static final int PLOTSTORED = 1534;//!
  private static final int STORE = 1535;//!
  private static final int PLOTLAST = 1536;//!
  /*!*/

  //  public static boolean TRACING;
  private static final int GROWTHMAX = 1000;

  static boolean hlVert() {
    return FULL && C.vb(HLEDGES);
  }

  static boolean labelEdges() {
    return C.vb(LABELEDGES);
  }

  //private static final boolean FULL = false;

  public void addControls() {
    C.sOpenTab("Discs");
    C.sStaticText("Plots guaranteed Voronoi diagram for uncertain discs");
    {
      C.sOpen();
      {
        C.sOpen();
        C.sCheckBox(PLOTSTD, "Standard",
            "Plot standard Voronoi diagram of sites", false);
        if (FULL) {
          C.sCheckBox(PLOTSTORED, "plot stored", null, false);
          C.sCheckBox(PLOTLAST, "editor last", "plot editor after other items",
              false);
          C.sNewColumn();
          C.sIntSpinner(STDK, "k:",
              "Order-k value for standard Voronoi diagram", 1, 5, 1, 1);
          C.sButton(STORE, "Store", null);
        }
        C.sClose();
      }

      {
        C.sOpen();
        C.sCheckBox(PLOTGUAR, "Guaranteed",
            "Plot Guaranteed Voronoi diagram of sites", true);
        if (FULL) {
          C.sNewColumn();
          C.sIntSpinner(GUARK, "k:",
              "Order-k value for guaranteed Voronoi diagram", 1, 5, 1, 1);
        }
        C.sClose();
      }

      {
        C.sOpen();
        C.sCheckBox(MINMAX, "MinMax", "Plots MinMax Voronoi diagram of sites",
            false);
        C.sNewColumn();
        C.sIntSpinner(MINMAXPOSS, "p:", "Disc to adjust possible site for", -1,
            100, -1, 1);

        C.sCheckBox(EFFDISCS, "plot eff", "Plots versions of the discs with "
            + "the radii they were  given to generate the MinMax diagram",
            false);

        C.sClose();
      }

      {
        C.sOpen();
        C.sCheckBox(SUBSETS, "Subset",
            "Plots guaranteed subset Voronoi diagram", false);
        if (FULL) {
          C.sNewColumn();
          C.sCheckBox(SUBSETCOLOR, "Different color",
              "Plots subset diagram in different color", true);
        }
        if (FULL)
          C.sCheckBox(NOCLIP, "No clipping", "Don't clip any bisectors", false);
        C.sClose();
      }
      {
        C.sOpen();
        C.sCheckBox(PLOTPOSS, "Possible cell",
            "Plots possible cell for a site", false);

        C.sNewColumn();
        C.sIntSpinner(PSITE, "site #:", null, 0, 1000, 0, 1);
        C.sClose();
      }
      {
        C.sOpen();

        C.sCheckBox(SPINES, "Spines", "Plots spines", false);
        C.sCheckBox(GDEL, "Delaunay", "Plots guaranteed Delaunay edges", false);
        if (FULL)
          C.sCheckBox(PLOTCRUST, "Guaranteed crust", null, false);

        if (FULL) {
          C.sNewColumn();

          C.sCheckBox(HLEDGES, "Highlight edges",
              "Show interfaces between edges", false);
          C.sCheckBox(LABELEDGES, "Label edges", null, false);

          C.sButton(PERTURB, "Perturb", null);
        }
        C.sClose();
      }
      //      C.sNewColumn();
      if (FULL) {
        C.sOpen("Smooth Growth");
        C.sButton(AGSTORE, "Store final",
            "Set current sites as final growth stage");
        C.sIntSlider(AGFRAME, "Frame", "Position in growth", 0, GROWTHMAX,
            GROWTHMAX, 1);
        C.sClose();
      }

      C.sClose();

      //      C.sHide();
      //      C.sCheckBox(NEWCONTAINMENT, "New containment method", null, true);
      C.sHide();
      C.sCheckBox(PRINTEDGES, "Print edges", "Display edges in text form",
          false);
      C.sHide();
      C.sCheckBox(BENDDUALS, "Bend dual edges",
          "Redirect dual edges to intersect Voronoi edge", false);
      C.sHide();
      C.sCheckBox(CLIPTOCIRCLE, "Clip to boundary",
          "Clip dual edges to exterior of site boundaries", false);
      C.sHide();
      C.sCheckBox(DUALCIRCLES, "Plot dual circles",
          "Plot empty circles associated with dual edges", false);

      C.sHide();
      C.sCheckBox(TEST, "Test", null, false);
    }
    C.sCloseTab();
  }

  public static DiscVornOper singleton = new DiscVornOper();

  private DiscVornOper() {
  }

  private void buildPossibleCell() {
    EdDisc[] circ = Main.getDiscs();
    int s = C.vi(PSITE);
    if (s < 10) {
      if (s < circ.length) {
        possCells.add(PossibleCell.buildPossibleCell(s));
      }
    } else {
      String q = Integer.toString(s);
      for (int i = 0; i < q.length(); i++) {
        int k = q.charAt(i) - '0';
        if (k < circ.length) {
          possCells.add(PossibleCell.buildPossibleCell(k));
        }
      }
    }
  }

  public void runAlgorithm() {
    // clear any objects constructed during last paint cycle
    showables.clear();
    possCells.clear();
    possDisc = null;
    crust = null;
    if (C.vb(PLOTPOSS))
      buildPossibleCell();
    else {
      if (C.vb(MINMAX)) {
        int site = C.vi(MINMAXPOSS);
        if (site >= 0) {
          EdDisc[] circ = Main.getDiscs();
          if (site < circ.length) {
            VornGraph g = PossibleCell.buildPossibleCell(site);
            g.setAppearance(MyColor.cRED, STRK_THICK, -1, false);
            showables.add(g);
            possDisc = circ[site];
            //            showables.add(new EdDisc()
          }
        }
      }
    }

    if (FULL && C.vb(PLOTCRUST)) {
      try {
        DArray a = CrustOper.build(Main.getDiscs());
        crust = a.getDArray(0);
      } catch (FPError e) {
        Tools.warn("caught FPError");
      }
    }

  }
  private DArray crust;
  private Renderable possDisc;
  //  private VornGraph pointSites;

  public void paintView() {

    final boolean db = false;

    if (db)
      Streams.out.println("DiscVornOper paintview");

    EdDisc[] circ = Main.getDiscs();
    //    OneCenterOper.singleton.plotSamples();
    //    GuarOneCenterOper.singleton.plotSamples();

    if (!(FULL && C.vb(PLOTLAST)))
      Editor.render();

    if (FULL && C.vb(PLOTSTORED) && savedGraph != null) {
      savedGraph.setAppearance(null, STRK_THIN, hlVert() ? MARK_DISC : -1,
          labelEdges());
      savedGraph.render(null, -1, -1);
    }

    if (C.vb(PLOTSTD)) {
      int k = 1;
      if (FULL)
        k = C.vi(STDK);
      V.pushStroke(TestBed.STRK_NORMAL);
      V.pushColor(MyColor.get(MyColor.GREEN, .4));
      Hyperbola[] vDiag = VornUtil.build(filterContained(circ), k,
          StandardDiscBisector.S);
      VornGraph g = new VornGraph(vDiag);
      lastPlotted = g;
      g.setAppearance(null, -1, hlVert() ? MARK_DISC : -1, labelEdges());
      g.render(null, -1, -1);
      V.popColor();
      V.popStroke();
    }

    int kGuar = C.vi(GUARK);
    if (C.vb(SUBSETS))
      kGuar = 1;

    if (C.vb(SUBSETS)) {
      if (db)
        Streams.out.println(" plotting subsets");
      Hyperbola[] hypList = VornUtil.buildSubset(circ, !(FULL && C.vb(NOCLIP)),
          SubsetDiscBisector.S);

      VornGraph vg = new VornGraph(hypList);

      boolean cf = FULL && C.vb(SUBSETCOLOR);

      vg.plot(cf ? MyColor.get(MyColor.BROWN, .9) : null, -1, hlVert(),
          labelEdges());
      lastPlotted = vg;

    }

    if (C.vb(PLOTGUAR)) {
      Hyperbola[] edges = VornUtil.build(circ, kGuar, GuaranteedDiscBisector.S);
      VornGraph g = new VornGraph(edges);
      g.plot(MyColor.get(MyColor.DARKGRAY), -1, hlVert(), labelEdges());
      lastPlotted = g;

      if (C.vb(PRINTEDGES))
        g.printEdges();
    }

    spines = null;
    if (C.vb(SPINES)) {
      //      V.pushColor(MyColor.get(MyColor.PURPLE, .7));
      if (db)
        Streams.out.println(" rendering spines:\n"
            + DArray.toString(getSpines()));

      VornGraph g = new VornGraph(getSpines());
      g.plot(MyColor.get(MyColor.PURPLE, .7), -1, false, labelEdges());
      //      g.render(true, false, labelEdges());
      //      V.popColor();
    }

    if (C.vb(GDEL)) {
      V.pushColor(MyColor.get(MyColor.MAGENTA, .7));
      Hyperbola[] spines = getSpines();
      V.pushStroke(STRK_NORMAL);
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
    if (crust != null)
      CrustOper.plot(Main.getDiscs(), crust);

    T.renderAll(possCells, MyColor.cRED, -1, -1);

    if (C.vb(MINMAX)) {
      constructMinMax();
      if (C.vb(EFFDISCS))
        T.renderAll(effDiscs, MyColor.cPURPLE, STRK_RUBBERBAND, -1);
    }
    T.renderAll(showables, null);
    T.render(possDisc, MyColor.cRED, STRK_THICK, -1);
    if (FULL && C.vb(PLOTLAST))
      Editor.render();
  }

  private DArray effDiscs;
  private void constructMinMax() {
    EdDisc[] circ = Main.getDiscs();
    effDiscs = new DArray();

    int possSite = C.vi(MINMAXPOSS);
    // modify radii 
    double maxRad = -1;
    for (int i = 0; i < circ.length; i++)
      maxRad = Math.max(maxRad, circ[i].getRadius());

    for (int i = 0; i < circ.length; i++) {
      EdDisc d = new EdDisc(circ[i]);
      circ[i] = d;
      double rad = maxRad;
      if (i == possSite)
        rad += d.getRadius();
      else
        rad -= d.getRadius();
      d.setRadius(rad);
      effDiscs.add(d);
    }

    int k = 1; // C.vi(MINMAXK);

    Hyperbola[] vDiag = VornUtil.build(filterContained(circ), //circ,
        k, StandardDiscBisector.S);
    VornGraph g = new VornGraph(vDiag);
    g.setAppearance(MyColor.cDARKGREEN, -1, hlVert() ? MARK_DISC : -1,
        labelEdges());
    showables.add(g);

    lastPlotted = g;
  }
  private DArray showables = new DArray();
  private DArray possCells = new DArray();
  private Hyperbola[] spines;

  private Hyperbola[] getSpines() {
    if (spines == null) {
      spines = Spine.construct(Main.getDiscs());
    }
    return spines;
  }

  public void processAction(TBAction a) {
    if (a.code == TBAction.CTRLVALUE) {
      switch (a.ctrlId) {
      case PERTURB:
        Main.perturbDiscs();
        break;
      case AGSTORE:
        sgStoreInitialValues();
        break;
      case AGFRAME:
        sgScaleSelected(C.vi(AGFRAME));
        break;
      case STORE:
        savedGraph = lastPlotted;
        C.setb(PLOTSTORED, true);
        break;
      }
    }
  }

  /**
   * Scale selected discs' radii
   * @param scaleFactor : scale factor 0..GROWTHMAX
   */
  private void sgScaleSelected(int scaleFactor) {

    sgInitIfNec();
    DArray a = Editor.editObjects(EdDisc.FACTORY, true, false);

    double f = C.vi(AGFRAME) / (double) GROWTHMAX;

    for (int i = 0; i < a.size(); i++) {
      EdDisc c = (EdDisc) a.get(i); //))EdDisc.isCompletedDisc(a.obj(i));
      c.setRadius(sgRadii[i] * f);
    }
  }

  private void sgInitIfNec() {
    DArray items = Editor.editObjects(EdDisc.FACTORY, true, false);
    boolean match = sgRadii != null && sgRadii.length == items.size();
    if (!match) {
      sgRadii = new double[items.size()];
      for (int i = 0; i < items.size(); i++) {
        EdDisc c = (EdDisc) items.get(i);
        sgRadii[i] = c.getRadius();
      }
    }
  }

  private void sgStoreInitialValues() {
    sgRadii = null;
    sgInitIfNec();
    C.seti(AGFRAME, GROWTHMAX);
  }

  private double[] sgRadii;

  private VornGraph savedGraph;

  private VornGraph lastPlotted;

  private static EdDisc[] filterContained(EdDisc[] sites) {
    int net = sites.length;
    for (int j = 0; j < sites.length; j++) {
      EdDisc c1 = sites[j];
      for (int k = j + 1; k < sites.length; k++) {
        EdDisc c2 = sites[k];
        {
          // if one circle contains the other, mark smaller 
          if (EdDisc.contains(c1, c2))
            c2.addFlags(DiscUtil.DISC_CONTAINED);
          else if (EdDisc.contains(c2, c1))
            c1.addFlags(DiscUtil.DISC_CONTAINED);
        }
      }
      if (DiscUtil.contained(c1))
        net--;
    }
    EdDisc[] ret = new EdDisc[net];
    int k = 0;
    for (int j = 0; j < sites.length; j++)
      if (!DiscUtil.contained(sites[j]))
        ret[k++] = sites[j];
    return ret;
  }
}
