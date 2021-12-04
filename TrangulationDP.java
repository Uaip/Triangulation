import java.text.DecimalFormat;
import java.util.Arrays;

public class TrangulationDP {
    static final double EARTH_RADIUS = 6378.137; // equatorial radius
    static final int kmToMeter = 1000; // km to meter

    // static double[][] point = new double[][] { { 102.746, 24.9555 }, { 102.872, 25.01466 }, { 102.879, 25.024444 },
    //         { 102.895, 25.050191 }, { 102.895, 25.095 }, { 102.863, 25.09831 }, { 102.768, 25.10243 },
    //         { 102.736, 25.10088 }, { 102.719, 25.096666 }, { 102.678, 25.05902 }, { 102.677, 25.05719 },
    //         { 102.676, 25.05218 }, { 102.675, 25.03436 }, { 102.675, 25.01318 }, { 102.675, 25.00996 },
    //         { 102.675, 24.99597 }, { 102.676, 24.99248 }, { 102.686, 24.98032 }, { 102.704, 24.97184 },
    //         { 102.74, 24.95721 }, { 102.746, 24.9555 } };

    // static double[][] point = new double[][] { { 0, 0 }, { 3, 0 }, { 4, 4 }, {
    // 2,8 }, {1, 9}, { 0, 4 }};

    static double[][] point = new double[][] { { 102.732, 24.98221 }, { 102.746,
    24.98604 }, { 102.749, 24.988197 },
    { 102.763, 25.0001 }, { 102.769, 25.01036 }, { 102.77, 25.01896 }, { 102.77,
    25.027222 },
    { 102.769, 25.033055 }, { 102.764, 25.04446 }, { 102.754, 25.05378 }, {
    102.744, 25.05965 },
    { 102.74, 25.06162 }, { 102.729, 25.06275 }, { 102.718, 25.06268 }, {
    102.714, 25.061944 },
    { 102.704, 25.05773 }, { 102.695, 25.05141 }, { 102.687, 25.04281 }, {
    102.685, 25.038284 },
    { 102.682, 25.03012 }, { 102.681, 25.0252 }, { 102.682, 25.01405 }, {
    102.686, 25.005124 },
    { 102.689, 25.00074 }, { 102.704, 24.98862 }, { 102.707, 24.9864 }, {
    102.709, 24.98496 },
    { 102.714, 24.9832 }, { 102.732, 24.98221 } };

    public static void main(String[] args) {
        int n = point.length;
        double t[][] = new double[n + 1][n + 1];
        int s[][] = new int[n + 1][n + 1];
        TrangulationDP test = new TrangulationDP();

        //get the perimeter to remove duplicate chord
        double Perimeter = test.getPerimeter();
        System.out.println("<<\t" + n + " sides triangulation\t>>");
        //SubOptimal Process
        test.HeuristicTriangulation(n, t, s);
        double subOptimal = (t[1][n - 1] + Perimeter) / 2;
        DecimalFormat df = new DecimalFormat("#.000000"); // round six places
        double res = Double.parseDouble(df.format(subOptimal));
        System.out.println("The SubOptimal Triangulation:\t" + res  + " meters");
        test.traceBack(1, n - 1, s);

        //Optimal process
        test.Initial(s);
        for(int i = 0;i < t.length;i++)
            Arrays.fill(t[i], 0);
        
        test.MinWeightTriangulation(n, t, s);
        double Optimal = (t[1][n - 1] + Perimeter) / 2; 
        res = Double.parseDouble(df.format(Optimal));
        System.out.println("The Optimal Triangulation:\t" + res + " meters");
        test.traceBack(1, n - 1, s);
    }

    private double getPerimeter() {
        double s = 0;
        for(int i = 0;i < point.length - 1;i++) {
            s += getDistance(point[i][0], point[i][1], point[i + 1][0], point[i + 1][1]);
        }
        s += getDistance(point[0][0], point[0][1], point[point.length - 1][0], point[point.length - 1][1]); //edge between start and end
        return s;
    }

    void Initial(int[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++)
                a[i][j] = -1;
        }
    }

    void traceBack(int i, int j, int[][] s) {
        if (i == j)
            return;
        traceBack(i, s[i][j], s);
        traceBack(s[i][j] + 1, j, s);
        System.out.println("(" + (i - 1) + "," + j + "," + s[i][j] + ")");
    }

    void MinWeightTriangulation(int n, double[][] t, int[][] s) {
        for (int i = 1; i <= n; i++)
            t[i][i] = 0;
        for (int r = 2; r <= n; r++) {
            for (int i = 1; i <= n - r + 1; i++) {
                int j = i + r - 1;
                t[i][j] = t[i + 1][j] + weightFunc(i - 1, i, j);
                s[i][j] = i;
                for (int k = i + 1; k < i + r - 1; k++) {
                    double u = t[i][k] + t[k + 1][j] + weightFunc(i - 1, k, j);
                    if (u < t[i][j]) {
                        t[i][j] = u;
                        s[i][j] = k;
                    }
                }
            }
        }
    }

    void HeuristicTriangulation(int n, double[][] t, int[][] s) {
        for (int i = 1; i <= n; i++)
            t[i][i] = 0;
        for (int r = 2; r <= n; r++) {
            for (int i = 1; i <= n - r + 1; i++) {
                int j = i + r - 1;
                t[i][j] = t[i + 1][j] + weightFunc(i - 1, i, j);
                s[i][j] = i;
                    for (int k = i + 1; k < i + r - 1 && k < i + 4; k++) {
                        double u = t[i][k] + t[k + 1][j] + weightFunc(i - 1, k, j);
                        if (u < t[i][j]) {
                            t[i][j] = u;
                            s[i][j] = k;
                        }
                    }
            }
        }
    }
 
    double weightFunc(int a, int b, int c) {
        if (c == point.length)
            c = 0;
        double sum = 0;
        sum += getDistance(point[a][0], point[a][1], point[b][0], point[b][1]);
        sum += getDistance(point[a][0], point[a][1], point[c][0], point[c][1]);
        sum += getDistance(point[b][0], point[b][1], point[c][0], point[c][1]);
        return sum;
    }

    double dis(double x1, double y1, double x2, double y2) {
        return Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    }

    double getDistance(double lng1, double lat1, double lng2, double lat2) {
        // lng -> longitude
        // lat -> latitude
        double radLng1 = rad(lng1);
        double radLng2 = rad(lng2);
        double radLat1 = rad(lat1);
        double radLat2 = rad(lat2);
        double dis = Math.acos(Math.cos(radLat1) * Math.cos(radLat2) * Math.cos(radLng1 - radLng2)
                + Math.sin(radLat1) * Math.sin(radLat2));
        dis *= EARTH_RADIUS;
        dis *= kmToMeter;
        DecimalFormat df = new DecimalFormat("#.000000"); // round six places
        // System.out.println(df.format(dis));
        double res = Double.parseDouble(df.format(dis));
        return res;
    }

    private double rad(double latOrlon) {
        return latOrlon * Math.PI / 180.0;
    }
}
