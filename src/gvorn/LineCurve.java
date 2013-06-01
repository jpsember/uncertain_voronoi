package gvorn;

import testbed.*;
import base.*;

public class LineCurve extends LineEqn implements IPlaneCurve {

  public LineCurve(FPoint2 p1, FPoint2 p2) {
    super(p1, p2);
  }
  
  public LineCurve(FPoint2 p0, double dir) {
    super(p0, MyMath.ptOnCircle(p0, dir, 1.0));
  }

  @Override
  public void render(double t0, double t1) {
    FPoint2 p0 = this.pt(t0);
    FPoint2 p1 = this.pt(t1);
    {
      if (MyMath.clipSegmentToRect(p0, p1, V.viewRect)) {
        V.drawLine(p0, p1);

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
  }
}
