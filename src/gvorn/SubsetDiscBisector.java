package gvorn;

import testbed.*;
import base.*;

public class SubsetDiscBisector extends GuaranteedDiscBisector {

  public void clipPair(Hyperbola hij, Hyperbola hik, EdPoint commonSite,
      DArray iPts) {
    final boolean db = false; // || debug;
    if (db)
      Streams.out.println("\n\nclipPair " + hij + " " + hik);

    //      DArray iPts = Hyperbola.findIntersections(hij, hik, null);
    switch (iPts.size()) {
    case 0:
      if (db)
        Streams.out.println("No intersects " + hij + " with " + hik);

      // Throw out the arc farthest from the common site, if
      // we're doing inner clipping; else, throw out the nearest.

      //        if (subsets) 
      {
        if (db)
          Streams.out.println("hij=" + hij + ", hik=" + hik);

//        if (hij.isEmpty()) {
//          Tools.warn(" was empty! " + hij);
//          break;
//        }
//        if (hik.isEmpty()) {
//          Tools.warn(" was empty! " + hik);
//          break;
//        }
//
        FPoint2 pi = hij.getSamplePoint(), pj = hik.getSamplePoint();

        if (hij.testPoint(pj) > 0) {
          hik.clipAll();
        } else if (hik.testPoint(pi) > 0) {
          hij.clipAll();
        }
        //            
        //            
        //            if (FPoint2.distanceSquared(origin, pi) < FPoint2
        //                    .distanceSquared(origin, pj)) {
        //              if (db) {
        //                Streams.out.println(" furthest, clipping all of " + hik);
        //                Streams.out.println("pi=" + pi + " pj=" + pj + " origin="
        //                    + origin);
        //
        //              }
        //              hik.clipAll();
        //            } else {
        //              if (db)
        //                Streams.out.println(" clipping all of " + hij);
        //              hij.clipAll();
        //            }

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
          FPoint2 i2 = ha.calcPoint(ti + (-.5));
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
      {
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
          if ( //(clipType == CLIP_FURTHEST || 
          true ^ hb.testPoint(i2) < 0) {
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
      break;
    }

  }

  //    public Hyperbola getBisector(VornSite a, VornSite b) {
  //      return guaranteedDiscBisector.getBisector(a, b);
  //      // TODO Auto-generated method stub
  //      return null;
  //    }

  //  public boolean guaranteed() {
  //    return true;
  //  }
  //
  //    public boolean overlap(VornSite a, VornSite b) {
  //      // TODO Auto-generated method stub
  //      return false;
  //    }

  //    public boolean symmetric() {
  //      // TODO Auto-generated method stub
  //      return false;
  //    }
  private SubsetDiscBisector() {
  }

  public static final ISiteBisector S = new SubsetDiscBisector();
}
