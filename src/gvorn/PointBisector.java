package gvorn;

import testbed.*;
import base.*;

public class PointBisector implements ISiteBisector {

  public static PointBisector S = new PointBisector();

  protected PointBisector() {
  }
  //  public static Hyperbola constructPointBisector(VornSite c1, VornSite c2) {
  //    Hyperbola h = new Hyperbola(c1.getOrigin(), c2.getOrigin());
  //    h.setLabel("" + c1.getLabel() + "/" + c2.getLabel());
  //    h.setData(c1, c2);
  //    return h;
  //  }

  public Hyperbola getBisector(EdPoint c1, EdPoint c2) {
    Hyperbola h = null;
    try {
      h = new Hyperbola(c1.getOrigin(), c2.getOrigin());
      h.setLabel("" + c1.getLabel() + "/" + c2.getLabel());
      h.setData(c1, c2);
      //      h = VornUtil.constructPointBisector(a, b);
    } catch (FPError e) {
      Tools.warn("FPError caught");
      h = null;
    }
    return h;
  }
//  public boolean guaranteed() {
//    return false;
//  }
//  public boolean overlap(VornSite a, VornSite b) {
//    return false;
//  }
//  public boolean symmetric() {
//    return true;
//  }

  public void clipPair(Hyperbola hij, Hyperbola hik, EdPoint commonSite,
      DArray iPts) {
    final boolean db = false; // || debug;

    int k = 1;
//    Tools.warn("k > 1 no longer supported");
    if (db)
      Streams.out.println("\n\nclipPair " + hij + " " + hik);

    switch (iPts.size()) {
    case 0:
      if (db)
        Streams.out.println("No intersects " + hij + " with " + hik);

      // Throw out the arc farthest from the common site, if
      // we're doing inner clipping; else, throw out the nearest.

      {
        if (k > 1) {
          Tools
              .warn("don't know how to handle no intersections properly, in k>1 case");
          // break;
        }

        FPoint2 pi = hij.calcPoint(0), pj = hik.calcPoint(0);
        FPoint2 origin = commonSite.getOrigin();
        if (

        //(clipType == CLIP_FURTHEST) ^

        FPoint2.distanceSquared(origin, pi) < FPoint2.distanceSquared(origin,
            pj)) {
          if (db) {
            Streams.out.println(" furthest, clipping all of " + hik);
            Streams.out.println("pi=" + pi + " pj=" + pj);

          }
          hik.clipAll();
        } else {
          if (db)
            Streams.out.println(" clipping all of " + hij);
          hij.clipAll();
        }
      }
      break;
    case 1:
      {
        FPoint2 ipt = iPts.getFPoint2(0);
        if (db)
          Streams.out.println("one intersect between " + hij + " and " + hik
              + " at " + ipt);

        for (int arm = 0; arm < 2; arm++) {
          Hyperbola ha = arm == 0 ? hij : hik;
          Hyperbola hb = arm == 0 ? hik : hij;

          double ti = ha.calcParameter(ipt);
          // move a little further along, and compare resulting point
          FPoint2 i2 = ha.calcPoint(ti + (.5));
          if ( //(clipType == CLIP_FURTHEST) ^ 
          hb.testPoint(i2) < 0) {
            ha.clip(ti, Hyperbola.CLIP_MAX);
            if (db)
              Streams.out.println(" clipping " + ha + " above " + ti);

          } else {
            ha.clip(Hyperbola.CLIP_MIN, ti);
            if (db)
              Streams.out.println(" clipping " + ha + " above " + ti);
          }
        }
      }
      break;
    case 2:
      throw new FPError("didn't expect two intersections");
    }

  }

  public boolean cellIsEmpty(EdPoint[] sites, EdPoint site) {
    return false;
  }

};
