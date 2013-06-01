package gvorn;

import java.awt.*;
import java.util.*;
import base.*;
import testbed.*;

public class VornTestOper implements TestBedOperation, Globals {
  /*! .enum  .private 1850
       scell  testmode filter _ _ type point disc seg poly iso0 iso1
  */

    private static final int SCELL            = 1850;//!
    private static final int TESTMODE         = 1851;//!
    private static final int FILTER           = 1852;//!
    private static final int TYPE             = 1855;//!
    private static final int POINT            = 1856;//!
    private static final int DISC             = 1857;//!
    private static final int SEG              = 1858;//!
    private static final int POLY             = 1859;//!
    private static final int ISO0             = 1860;//!
    private static final int ISO1             = 1861;//!
/*!*/

  public void addControls() {
    Tools.unimp("create composite bisector from segment pieces");
    C.sOpenTab("VTest");
    C.sStaticText("Tests new Vorn diagram functions");

    {
      {
        C.sOpen();
        C.sIntSpinner(SCELL, "Standard cell", "Plot standard Voronoi cell", 0,
            20, 0, 1);
        C.sOpenComboBox(TYPE, "Sites", "Select site type", false);
        C.sChoice(POINT, "point");
        C.sChoice(DISC, "disc");
        C.sChoice(SEG, "segment");
        C.sChoice(POLY, "polygon");
        C.sCloseComboBox();
        C.sCheckBox(TESTMODE, "test", null, false);
        C.sCheckBox(FILTER, "filter", null, false);
        C.sIntSpinner(ISO0, "iso0", null, 0, 20, 0, 1);
        C.sIntSpinner(ISO1, "iso1", null, 0, 20, 0, 1);
        C.sClose();
      }
    }

    C.sCloseTab();
  }

  public void processAction(TBAction a) {
    switch (a.code) {
    case TBAction.HOVER:
      if (C.vb(TESTMODE)) {
        hoverLoc = a.loc;
        V.repaint();
      }
      break;
    }
  }

  public void runAlgorithm() {
    trace = new DArray();
    pa = null;
    sites = null;
    cell = null;
    if (C.vb(TESTMODE)) {
      do {
        EdDisc[] d = Main.getDiscs();
        if (d.length < 3)
          break;
        VornDiscSite s0 = new VornDiscSite(d[0]);
        VornDiscSite s1 = new VornDiscSite(d[1]);
        VornDiscSite s2 = new VornDiscSite(d[2]);

        VornDiscBisector b1 = new VornDiscBisector(s0, s1);
        VornDiscBisector b2 = new VornDiscBisector(s0, s2);

        trace.add(b1);
        trace.add(b2);

        DArray tv = VUtil2.intersect(b1, b2, C.vb(FILTER));
        for (int i = 0; i < tv.size(); i += 2) {
          double t0 = tv.getDouble(i);
          double t1 = tv.getDouble(i + 1);
          trace.add(b1.pointAt(t0));
          trace.add(b2.pointAt(t1));
        }

        //      
        //      EdPoint[] p = Main.getPoints();
        //      if (p.length < 2)
        //        break;
        //      FPoint2 s0 = p[0].getOrigin();
        //      FPoint2 s1 = p[1].getOrigin();
        //
        //      pa = new Hyperb(s0, s1, FPoint2.distance(s0, s1) * .3);
        //      
        //      
        //      EdSegment[] s = Main.getSegments();
        //      if (p.length < 1 || s.length < 1)
        //        return;
        //
        //      pa = new Parabola(s[0].getPoint(0), s[0].getPoint(1), p[0].getPoint(0));
      } while (false);
      return;
    }

    VUtil2.setTraceArray(trace);

    sites = new DArray();

    switch (C.vi(TYPE)) {
    default:
      Tools.unimp();
      break;
    case POINT:
      {
        EdPoint[] d = Main.getPoints();

        for (int i = 0; i < d.length; i++) {
          sites.add(new VornPointSite(d[i].getOrigin()));
        }
      }
      break;
    case DISC:
      {
        EdDisc[] d = Main.getDiscs();

        // filter out discs that lie within others
        {
          DArray f = new DArray();
          iLoop: for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d.length; j++) {
              if (j == i || !EdDisc.contains(d[j], d[i]))
                continue;
              continue iLoop;
            }
            f.add(d[i]);
          }
          d = (EdDisc[]) f.toArray(EdDisc.class);
        }

        for (int i = 0; i < d.length; i++) {
          sites.add(new VornDiscSite(d[i].getOrigin(), d[i].getRadius(), d[i]
              .getLabel()));
        }
      }
      break;
    case SEG:
      {
        EdSegment[] d = Main.getSegments();
        for (int i = 0; i < d.length; i++) {
          sites.add(new VornSegSite(d[i].getPoint(0), d[i].getPoint(1), d[i]
              .getLabel()));
        }
     //   VornSegBisector.setIsolatedSegments(C.vi(ISO0) - 1, C.vi(ISO1) - 1);
      }
      break;
    }

    cell = new DArray();
    DArray[] bLists = VUtil2.buildVornDiagram(sites);
    for (int i = 0; i < bLists.length; i++) {
      VornSite si = (VornSite) sites.get(i);
      DArray bl = bLists[i];
      for (int j = 0; j < bl.size(); j++) {
        VornBisector bs = (VornBisector) bl.get(j);
        if (bs.otherSite(si).id() < si.id())
          continue;
        cell.add(bs);
      }
    }

  }
  public void paintView() {
    Editor.render();
    T.renderAll(sites, MyColor.cBLUE);
    T.renderAll(cell, MyColor.cDARKGREEN, -1, Globals.MARK_NONE);
    T.renderAll(trace, MyColor.cPURPLE);

    if (hoverLoc != null && pa != null) {
      Hyperb h = (Hyperb) pa;

      h.render(-100, 100);
      //      T.render(pa);
      double t = pa.parameterFor(hoverLoc);
      FPoint2 pt = pa.pt(t);

      double eval = h.eval(hoverLoc.x, hoverLoc.y);

      T.show("h=" + hoverLoc + " t=" + Tools.f(t) + "\nh'=" + pt + "\neval="
          + eval, MyColor.cRED, 0, 100, 40 | Globals.TX_CLAMP, .8);
      //          hoverLoc, 40, .7);
      V.mark(pt);
    }

  }

  private DArray sites, cell;
  private DArray trace;
  private IPlaneCurve pa;
  private FPoint2 hoverLoc;
}
