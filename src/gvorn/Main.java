package gvorn;

import base.*;
import testbed.*;
import java.util.*;

public class Main extends TestBed {
  /*! .enum  .private .prefix G_ 4000
     _  _ _ _ togglediscs _ makesupported maketangent
     _ _ _   random rndtest 
  */

  private static final int G_TOGGLEDISCS = 4004;//!
  private static final int G_MAKESUPPORTED = 4006;//!
  private static final int G_MAKETANGENT = 4007;//!
  private static final int G_RANDOM = 4011;//!
  private static final int G_RNDTEST = 4012;//!
  /* !*/

  public static void main(String[] args) {
    new Main().doMainGUI(args);
  }
public static final boolean FULL = false;

  // -------------------------------------------------------
  // TestBed overrides
  // -------------------------------------------------------

  public void addOperations() {
    addOper(DiscVornOper.singleton);
    addOper(PolygonVornOper.singleton);
    addOper(new GabrielOper());
    addOper(new GeneratorOper());

    //addOper(new VornTestOper());
  }
  public void addControls() {
    C.sOpen();
    C.sButton(G_RANDOM, "Random", "Generate random discs (or polygons)");
    if (FULL)
    C.sCheckBox(G_RNDTEST, "Test", "Repeatedly generate random discs", false);
    C.sClose();
  }
  public void initEditor() {
    Editor.addObjectType(EdPolygon.FACTORY);
    Editor.addObjectType(EdDisc.FACTORY);
    Editor.addObjectType(EdSegment.FACTORY);
    Editor.addObjectType(EdDiameter.FACTORY);
    Editor.addObjectType(EdPoint.FACTORY);

    Editor.openMenu();
    C.sMenuItem(G_TOGGLEDISCS, "Toggle discs/points", "!^t");
    C.sMenuItem(G_MAKETANGENT, "Set disc tangent", "!^3"); //"!^g");
    C.sMenuItem(G_MAKESUPPORTED, "Set disc supported", "!^4"); //"!^u");
    Editor.closeMenu();
  }

  public void processAction(TBAction a) {
    if (a.code == TBAction.CTRLVALUE) {
      switch (a.ctrlId) {
      case G_RANDOM:
        GeneratorOper.generateRandom();
        break;
      case G_RNDTEST:
        if (C.vb(G_RNDTEST))
          GeneratorOper.generateRandom();
        break;

      case G_TOGGLEDISCS:
        {
          for (int i = 0; i < discs2.length; i++) {
            EdDisc c = discs2[i];
            if (!c.isSelected())
              continue;
            c.togglePointMode();
          }
        }
        break;
      case G_MAKESUPPORTED:
        makeSupported();
        break;
      case G_MAKETANGENT:
        makeTangent();
        break;
      }
    }
  }

  /**
   * Make highlighted inactive discs supported by best-fit pair of candidates
   */
  private void makeSupported() {

    DArray a = Editor.editObjects(EdDisc.FACTORY, true, false);
    for (int k = 0; k < a.size(); k++) {
      EdDisc c = (EdDisc) a.get(k);

      {
        if (c.isActive())
          continue;

        // find best two candidates for supporting this disc
        DArray can = getCandidate(c);
        if (can.size() < 2)
          continue;

        EdDisc ca = (EdDisc) can.get(0);
        EdDisc cb = (EdDisc) can.get(1);
        boolean ia = DiscUtil.itan(c, ca);

        c.setPoint(0,
            DiscUtil.supportingHyperbola(c, ca, cb).snap(c.getOrigin()));

        c.setRadius(FPoint2.distance(c.getOrigin(), ca.getOrigin())
            + (ia ? ca.getRadius() : -ca.getRadius()));
      }
    }
  }

  private static DArray getCandidate(final EdDisc c) {
    // find best two candidates for supporting this disc
    DArray can = new DArray(discs);
    can.sort(new Comparator() {
      public int compare(Object arg0, Object arg1) {
        EdDisc c1 = (EdDisc) arg0, c2 = (EdDisc) arg1;

        double d1 = DiscUtil.itan(c, c1) ? DiscUtil.itanDist(c, c1) : DiscUtil
            .otanDist(c, c1);
        double d2 = DiscUtil.itan(c, c2) ? DiscUtil.itanDist(c, c2) : DiscUtil
            .otanDist(c, c2);

        return (int) Math.signum(d1 - d2);
      }
    });
    return can;
  }

  /**
   * Make highlighted inactive discs tangent to best-fit candidate
   */
  private void makeTangent() {
    DArray a = Editor.editObjects(EdDisc.FACTORY, true, false);
    for (int k = 0; k < a.size(); k++) {
      EdDisc c = (EdDisc) a.get(k);

      // find best two candidates for supporting this disc
      DArray can = getCandidate(c);
      if (can.size() < 1)
        continue;

      EdDisc ca = (EdDisc) can.get(0);
      boolean ia = DiscUtil.itan(c, ca);
      c.setRadius(FPoint2.distance(c.getOrigin(), ca.getOrigin())
          + (ia ? ca.getRadius() : -ca.getRadius()));
    }
  }

  public void setParameters() {
    parms.appTitle = "Guaranteed Voronoi Diagrams";
    parms.menuTitle = "GVorn";
    parms.fileExt = "dat";
  }

  public void paintView() {
    discs = null;
    discs2 = null;
    polygons = null;
    points = null;
    segs = null;

    super.paintView();

    //    if (TestBed.oper() != EnvOper.singleton)
    //      EnvOper.singleton.plotSamples();

    if (FULL && C.vb(G_RNDTEST)) {
      if (T.lastEvent() == null) {
        GeneratorOper.generateRandom();
        V.repaint();
      } else
        C.setb(G_RNDTEST, false);
    }
  }

  public static DArray getPolygons() {
    if (polygons == null) {
      polygons = Editor.readObjects(EdPolygon.FACTORY, false, true);
      for (int i = 0; i < polygons.size(); i++) {
        EdPolygon ed = (EdPolygon) polygons.get(i);
        ed = EdPolygon.normalize(ed);
        polygons.set(i, ed);
      }
    }
    return polygons;
  }

  public static void perturbDiscs() {
    Matrix m = new Matrix(3);
    m.setIdentity();

    //      Matrix.identity(3, null);
    m.translate(-50, -50);
    m.rotate(MyMath.radians(1));
    m.translate(50, 50);

    for (Iterator it = Editor.editObjects(EdDisc.FACTORY, false, false)
        .iterator(); it.hasNext();) {
      EdDisc ed = (EdDisc) it.next();
      FPoint2 loc = ed.getOrigin();
      loc = m.apply(loc, null);
      ed.setPoint(0, loc);
    }
    Editor.unselectAll();
  }

  public static EdDisc[] getDiscs2() {
    getDiscs();
    return discs2;
  }

  public static EdDisc[] getDiscs() {
    if (discs == null) {
      DArray a = Editor.readObjects(EdDisc.FACTORY, false, true);
      DArray b = Editor.readObjects(EdDisc.FACTORY, false, false);

      discs = (EdDisc[]) a.toArray(EdDisc.class);
      discs2 = (EdDisc[]) b.toArray(EdDisc.class);

      for (int i = 0; i < discs2.length; i++) {
        discs2[i].clearFlags(DiscUtil.DISC_OVERLAPPING
            | DiscUtil.DISC_CONTAINED);
      }

    }
    return discs;
  }
  public static EdPoint[] getPoints() {
    if (points == null) {
      DArray a = Editor.readObjects(EdPoint.FACTORY, false, true);
      points = (EdPoint[]) a.toArray(EdPoint.class);
    }
    return points;
  }
  public static EdSegment[] getSegments() {
    if (segs == null) {
      DArray a = Editor.readObjects(EdSegment.FACTORY, false, true);
      segs = (EdSegment[]) a.toArray(EdSegment.class);
    }
    return segs;
  }

  private static EdDisc[] discs;
  private static EdDisc[] discs2;
  private static DArray polygons;
  private static EdPoint[] points;
  private static EdSegment[] segs;

}
