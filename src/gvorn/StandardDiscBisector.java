package gvorn;

import testbed.*;
import base.*;

public class StandardDiscBisector extends PointBisector {

  public Hyperbola getBisector(EdPoint a, EdPoint b) {
    Hyperbola h = null;
    do {
      //    try {
      EdDisc c1 = (EdDisc) a, c2 = (EdDisc) b;
//      if (c1.contains(c2) || c2.contains(c1))
//        break;
      if (c1.getRadius() == c2.getRadius()) {
        h = PointBisector.S.getBisector(c1, c2);
        break;
      }
      double dist = FPoint2.distance(c1.getOrigin(), c2.getOrigin());
      double adj = .5 * (-c2.getRadius() + c1.getRadius() + dist);
      if (adj > 0 && adj < dist) {
        h = new Hyperbola(c1.getOrigin(), c2.getOrigin(), adj);
        h.setLabel(c1.getLabel() + "/" + c2.getLabel());
        h.setData(c1, c2);
      }
      //    } catch (FPError e) {
      //      Tools.warn("FPError caught");
      //    }
    } while (false);
    return h;
  }

  public void clipPair(Hyperbola hij, Hyperbola hik, EdPoint commonSite,
      DArray iPts) {

    if (iPts.size() < 2) {
      super.clipPair(hij, hik, commonSite, iPts);
      return;
    }

    final boolean db = false;

    FPoint2 ipt0 = iPts.getFPoint2(0);
    FPoint2 ipt1 = iPts.getFPoint2(1);
    if (db)
      Streams.out.println("two intersects between " + hij + " and " + hik
          + " at\n " + ipt0 + ", " + ipt1);

    for (int arm = 0; arm < 2; arm++) {
      Hyperbola ha = arm == 0 ? hij : hik;
      Hyperbola hb = arm == 0 ? hik : hij;
      double t0 = ha.calcParameter(ipt0);
      double t1 = ha.calcParameter(ipt1);
      FPoint2 i2 = ha.calcPoint((t0 + t1) * .5);
      if (hb.testPoint(i2) < 0) {
        ha.clip(t0, t1);
        if (db)
          Streams.out.println(" clipping from " + t0 + " to " + t1);
      } else {
        if (db)
          Streams.out.println(" clipping all but " + t0 + " to " + t1);
        ha.clipAllBut(t0, t1);
      }
    }
  }

  public static final ISiteBisector S = new StandardDiscBisector();
}