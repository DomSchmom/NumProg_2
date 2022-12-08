import java.util.Arrays;

/**
 * Die Klasse CubicSpline bietet eine Implementierung der kubischen Splines. Sie
 * dient uns zur effizienten Interpolation von aequidistanten Stuetzpunkten.
 *
 * @author braeckle
 *
 */
public class CubicSpline implements InterpolationMethod {

    /** linke und rechte Intervallgrenze x[0] bzw. x[n] */
    double a, b;

    /** Anzahl an Intervallen */
    int n;

    /** Intervallbreite */
    double h;

    /** Stuetzwerte an den aequidistanten Stuetzstellen */
    double[] y;

    /** zu berechnende Ableitunge an den Stuetzstellen */
    double yprime[];

    /**
     * {@inheritDoc} Zusaetzlich werden die Ableitungen der stueckweisen
     * Polynome an den Stuetzstellen berechnet. Als Randbedingungen setzten wir
     * die Ableitungen an den Stellen x[0] und x[n] = 0.
     */
    @Override
    public void init(double a, double b, int n, double[] y) {
        this.a = a;
        this.b = b;
        this.n = n;
        h = ((double) b - a) / (n);

        this.y = Arrays.copyOf(y, n + 1);

        /* Randbedingungen setzten */
        yprime = new double[n + 1];
        yprime[0] = 0;
        yprime[n] = 0;

        /* Ableitungen berechnen. Nur noetig, wenn n > 1 */
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * getDerivatives gibt die Ableitungen yprime zurueck
     */
    public double[] getDerivatives() {
        return yprime;
    }

    /**
     * Setzt die Ableitungen an den Raendern x[0] und x[n] neu auf yprime0 bzw.
     * yprimen. Anschliessend werden alle Ableitungen aktualisiert.
     */
    public void setBoundaryConditions(double yprime0, double yprimen) {
        yprime[0] = yprime0;
        yprime[n] = yprimen;
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * Berechnet die Ableitungen der stueckweisen kubischen Polynome an den
     * einzelnen Stuetzstellen. Dazu wird ein lineares System Ax=c mit einer
     * Tridiagonalen Matrix A und der rechten Seite c aufgebaut und geloest.
     * Anschliessend sind die berechneten Ableitungen y1' bis yn-1' in der
     * Membervariable yprime gespeichert.
     *
     * Zum Zeitpunkt des Aufrufs stehen die Randbedingungen in yprime[0] und yprime[n].
     * Speziell bei den "kleinen" Faellen mit Intervallzahlen n = 2
     * oder 3 muss auf die Struktur des Gleichungssystems geachtet werden. Der
     * Fall n = 1 wird hier nicht beachtet, da dann keine weiteren Ableitungen
     * berechnet werden muessen.
     */
    public void computeDerivatives() {
        /* TODO: diese Methode ist zu implementieren */
        if(n == 2){

        }

        h = (b - a) / n;
        double[] yNew = new double[y.length - 2];
        for (int i = 3; i < y.length; i++) {
            yNew[i-2] = (y[i] + y[i-2]) * 3/h;
        }
        yNew[0] = (y[2] - y[0] - yprime[0] * h/3) * 3/h;
        yNew[yNew.length - 1] = (y[y.length-1] - y[y.length-3] - yprime[n] * h/3) * 3/h;

        double lower[] = new double[yNew.length-1];
        double diag[] = new double[yNew.length];
        double upper[] = new double[yNew.length-1];
        Arrays.fill(lower, 1.0);
        Arrays.fill(diag, 4.0);
        Arrays.fill(upper, 1.0);
        TridiagonalMatrix tri = new TridiagonalMatrix(lower, diag, upper);

        double[] solved = tri.solveLinearSystem(yNew);
        for(int i = 0; i < solved.length; i++) {
            yprime[i+1] = solved[i];
        }
    }

    /**
     * {@inheritDoc} Liegt z ausserhalb der Stuetzgrenzen, werden die
     * aeussersten Werte y[0] bzw. y[n] zurueckgegeben. Liegt z zwischen den
     * Stuetzstellen x_i und x_i+1, wird z in das Intervall [0,1] transformiert
     * und das entsprechende kubische Hermite-Polynom ausgewertet.
     */
    @Override
    public double evaluate(double z) {
        /* TODO: diese Methode ist zu implementieren */
        if(z < a)
            return a;
        if(z > b)
            return b;
        int i = (int) ((z-a) / h);
        double t = (z - getXi(i)) / h;
        double q = y[i] * h0(t) + y[i+1] * h1(t) + h * yprime[i] * h2(t) + h * yprime[i+1] * h3(t);
        return q;
    }

    private double h0(double t) {
        return (1 - 3 * t * t + 2 * t * t * t);
    }

    private double h1(double t) {
        return (3 * t * t - 2 * t * t * t);
    }

    private double h2(double t) {
        return (t - 2 * t * t + t * t * t);
    }

    private double h3(double t) {
        return (-1 * t * t + t * t * t);
    }

    private double getXi(int i) {
        return a + h * i;
    }
}
