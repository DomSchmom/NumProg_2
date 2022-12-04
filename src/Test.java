import java.util.Arrays;
import java.util.List;

import dft.DFT;
import dft.IFFT;
import dft.Complex;

public class Test {

    /**
     * @param args
     */
    public static void main(String[] args) {
        //testLinear();
        testNewton();
        //testSplines();
        //testFFT();
    }

    private static void testLinear() {

        double[] x = { -1, 1, 3 };
        double[] y = { -3, 1, -3 };
        LinearInterpolation l = new LinearInterpolation();
        l.init(x,y);

        System.out.println("Ergebnis Linear: " + l.evaluate(0.5));
        System.out.println("-------------------------------");
    }

    private static void testNewton() {

        double[] x = { -1, 1, 3 };
        double[] y = { -3, 1, -3 };
        NewtonPolynom p = new NewtonPolynom(x, y);

        System.out.println("Ergebnis Newton: " + p.evaluate(1));
        System.out.println("-------------------------------");
    }

    public static void testSplines() {
        CubicSpline spl = new CubicSpline();
        double[] y = { 2, 0, 2, 3 };
        spl.init(-1, 2, 3, y);
        spl.setBoundaryConditions(9, 0);
        System.out.println(Arrays.toString(spl.getDerivatives())
                + " sollte sein: [9.0, -3.0, 3.0, 0.0].");
    }

    public static void testFFT() {
        System.out.println("Teste Fast Fourier Transformation");

        double[] v = new double[4];
        for (int i = 0; i < 4; i++)
            v[i] = i + 1;
        Complex[] c = dft.DFT.dft(v);
        Complex[] v2 = dft.IFFT.ifft(c);

        for (int i = 0; i < 4; i++) {
            System.out.println(v2[i]);
        }
        System.out
                .println("Richtig waeren gerundet: Eigene Beispiele ueberlegen");

        System.out.println("*************************************\n");
    }
}
