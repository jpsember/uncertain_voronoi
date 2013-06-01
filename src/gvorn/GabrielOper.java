package gvorn;

import java.awt.*;
import base.*;
import testbed.*;
import static gvorn.Main.FULL;

public class GabrielOper implements TestBedOperation, Globals {
  /*! .enum  .private 2100
      _ plotupper plotgabdiscs samples dosamp
      disc1 
      disc2 disc1angle disc2angle
  */

    private static final int PLOTUPPER        = 2101;//!
    private static final int PLOTGABDISCS     = 2102;//!
    private static final int SAMPLES          = 2103;//!
    private static final int DOSAMP           = 2104;//!
    private static final int DISC1            = 2105;//!
    private static final int DISC2            = 2106;//!
    private static final int DISC1ANGLE       = 2107;//!
    private static final int DISC2ANGLE       = 2108;//!
/*!*/
//private static final boolean FULL = false;

  public void addControls() {
    C.sOpenTab("Gabriel");
    C.sStaticText("Constructs Gabriel zones for a pair of uncertain discs; tests if third disc intersects this zone");
    
    {
      C.sOpen();
      if (FULL)
      C.sCheckBox(DOSAMP, "generate samples",
          "Use sampling method to verify Guaranteed Gabriel Zone", false);
      if (FULL)
        C.sCheckBox(PLOTUPPER, "plot upper function only", null, false);

      if (FULL) {
      C.sNewColumn();

      C.sIntSpinner(SAMPLES, "count:", null, 5, 180, 20, 1);
      }
      C.sClose();
    }
    {
      C.sOpen();
      C.sCheckBox(DISC1, "gabriel disc #1", null, false);
      C.sCheckBox(DISC2, "gabriel disc #2", null, false);
      C.sCheckBox(PLOTGABDISCS, "plot gabriel discs", null, false);

      C.sNewColumn();

      C.sIntSlider(DISC1ANGLE, "angle:", null, 0, 360, 45, 5);
      C.sIntSlider(DISC2ANGLE, "angle:", null, 0, 360, 45, 5);

      C.sClose();
    }
    C.sCloseTab();
  }

  public void paintView() {
    Editor.render();
    do {
      EdDisc[] circ = Main.getDiscs();
      if (circ.length < 2)
        break;

      fn = new GabrielZoneFunction(circ[0], circ[1]);
      fn.setUpperOnly(FULL && C.vb(PLOTUPPER));

      T.render(fn);

      if (FULL && C.vb(DOSAMP))
        fn.sample(C.vi(SAMPLES));

      if (circ.length >= 3) {
        EdDisc d = circ[2];

        FPoint2 bndPt = new FPoint2();
        LineEqn[] ln = new LineEqn[1];

        V.pushColor(MyColor.cRED);
        if (fn.discIntersectsZone(d, bndPt, ln)) {
          V.pushColor(new Color(255, 200,200));
          V.fillCircle(d.getOrigin(), d.getRadius());
          V.pop();
        }
  
        LineEqn eqn = ln[0];
        double t = eqn.parameterFor(bndPt);
        final double TANGENT_LEN = 20.0;
        EdSegment.plotDirectedLine(eqn.pt(t-TANGENT_LEN), eqn.pt(t+TANGENT_LEN));
        V.mark(bndPt, MARK_CIRCLE);
        V.pop();
      }

      {
        V.pushColor(MyColor.cDARKGRAY);
        FPoint2 c0 = circ[0].getOrigin();
        FPoint2 c1 = circ[1].getOrigin();
        V.mark(c0);
        V.mark(c1);
        LineEqn line = new LineEqn(c0, c1);

        V.pushStroke(STRK_RUBBERBAND);

        V.drawLine(line.pt(-200), line.pt(200));
        V.pop(2);
      }

      for (int k = 0; k < 2; k++) {
        if (!C.vb(DISC1 + k))
          continue;

        EdDisc d0 = circ[k];
        EdDisc d1 = circ[1 - k];

        fn = new GabrielZoneFunction(d0, d1);

        V.pushColor(k == 0 ? MyColor.cDARKGREEN : MyColor.cRED);

        double theta0 = MyMath.radians(C.vi(DISC1ANGLE + k) - 45);

        FPoint2 pt = fn.pointAt(theta0, false);
        V.mark(pt, MARK_DISC);

        FPoint2[] diam = new FPoint2[2];

        double theta = theta0 + fn.getAxisAngle();

        for (int i = 0; i < 2; i++) {
          EdDisc disc = (i == 0) ? d0 : d1;
          FPoint2 orig = disc.getOrigin();
          double rad = disc.getRadius();
          FPoint2 pa = MyMath.ptOnCircle(orig, theta
              + (i == 1 ? 0 : Math.PI / 2), rad);

          V.mark(pa, MARK_X);
          diam[i] = pa;
          FPoint2 pext = FPoint2.interpolate(pt, pa, 1.5);
          V.pushStroke(STRK_NORMAL);
          EdSegment.plotDirectedLine(pt, pext);
          V.pop();
        }
        if (C.vb(PLOTGABDISCS)) {
          V.pushStroke(STRK_THIN);
          V.drawCircle(FPoint2.midPoint(diam[0], diam[1]), diam[0]
              .distance(diam[1]) / 2);
          V.pop();
        }
        V.pop();
      }
    } while (false);
  }

  public void processAction(TBAction a) {
  }

  public void runAlgorithm() {
  }

  private GabrielZoneFunction fn;
}
