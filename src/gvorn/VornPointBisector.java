package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class VornPointBisector extends VornBisector {
  public  IPlaneCurve curve(int i){return lineEqn0;}


  public VornPointBisector(VornPointSite sa, VornPointSite sb) {
    super(sa, sb);

    final boolean db = false;

    FPoint2 mid = FPoint2.midPoint(sa.location(), sb.location());
    double theta = MyMath.polarAngle(sa.location(), sb.location()) + Math.PI
        / 2;
    lineEqn0 = new LineCurve(mid, theta);
    if (db && T.update())
      T.msg("constructed bisector between " + sa.id() + " and " + sb.id()
          + T.show(this) + " sa=" + sa + " sb=" + sb);
  }

  @Override
  public FPoint2 pointAt(double t) {
    return lineEqn0.pt(t);
  }

  @Override
  public double parameterFor(FPoint2 pt) {
    double t = lineEqn0.parameterFor(pt);
    return t;
  }

  @Override
  public void render(Color c, int stroke, int markType) {
    for (int i = 0; i < this.nComponents(); i++) {
      double t0 = this.cStart(i);
      double t1 = this.cEnd(i);
      renderLine(pointAt(t0), pointAt(t1) );

    }
  }


  /**
   * @deprecated use VUtil2.renderLine
   * @param p0
   * @param p1
   */public static void renderLine(FPoint2 p0, FPoint2 p1 ) {
//    V.pushStroke(stroke, Globals.STRK_NORMAL);
//    V.pushColor(c, MyColor.cDARKGREEN);

    //    FPoint2 p0 = pointAt(t0);
    //    FPoint2 p1 = pointAt(t1);

    {
      if (MyMath.clipSegmentToRect(p0, p1, V.viewRect)) {
        V.drawLine(p0, p1);
        //        EdSegment.plotDirectedLine(p0, p1);

        if (false) {
          FPoint2 np = MyMath.ptOnCircle(FPoint2.midPoint(p0, p1), MyMath
              .polarAngle(p0, p1)
              + Math.PI / 2, 4);
          V.pushScale(.5);
          //   V.draw(Tools.f(t0) + "..." + Tools.f(t1), np);
          V.pop();
        }

      }
    }

//    V.pop(2);
  }

//  @Override
//  public IPlaneCurve getPlaneCurve() {
//    return lineEqn0;
//  }
  private LineEqn lineEqn0;

}
