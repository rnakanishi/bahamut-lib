import numpy
import pylab as plt

def CatmullRomSpline(P0, P1, P2, P3, Px):
  """
  P0, P1, P2, and P3 should be (x,y) point pairs that define the Catmull-Rom spline.
  nPoints is the number of points to include in this curve segment.
  """
  # Convert the points to numpy so that we can do array multiplication
  P0, P1, P2, P3 = map(numpy.array, [P0, P1, P2, P3])

  # Calculate t0 to t4
  alpha = 0.5
  def tj(ti, Pi, Pj):
    xi, yi = Pi
    xj, yj = Pj
    return ( ( (xj-xi)**2 + (yj-yi)**2 )**0.5 )**alpha + ti

  t0 = 0
  t1 = tj(t0, P0, P1)
  t2 = tj(t1, P1, P2)
  t3 = tj(t2, P2, P3)

  t = (Px - P1[0])/(P2[0] - P1[0]) * (t2 - t1) + t1;

  # Reshape so that we can multiply by the points P0 to P3
  # and get a point for each value of t.
#   t = t.reshape(len(Px),1)
#   print(t)
  A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
  A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
  A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
#   print(A1)
#   print(A2)
#   print(A3)
  B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
  B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3

  C  = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2
  return C


# Define a set of points for curve to go through
Points = [[0,1.5],[2,2],[3,1],[4,0.5]]

v = CatmullRomSpline(Points[0], Points[1], Points[2], Points[3], 2.1)
print(v)