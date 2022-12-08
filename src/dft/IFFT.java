package dft;

import java.util.Arrays;

/**
 * Schnelle inverse Fourier-Transformation
 *
 * @author Sebastian Rettenberger
 */
public class IFFT {
    /**
     * Schnelle inverse Fourier-Transformation (IFFT).
     *
     * Die Funktion nimmt an, dass die Laenge des Arrays c immer eine
     * Zweierpotenz ist. Es gilt also: c.length == 2^m fuer ein beliebiges m.
     */
    public static Complex[] ifft(Complex[] c) {
        // TODO: diese Methode ist zu implementieren
        return ifftRec(c, c.length);
    }

    private static Complex[] ifftRec(Complex[] c, int n) {
        Complex[] v = new Complex[n];
        if(n == 1){
            v[0] = c[0];
        }
        else{
            int m = n/2;
            Complex[] z1Array = new Complex[n/2];
            Complex[] z2Array = new Complex[n/2];
            for(int i = 0; i < n; i+=2) {
                z1Array[i/2] = c[i];
            }
            for(int i = 1; i < n; i+=2) {
                z2Array[i/2] = c[i];
            }

            Complex[] z1 = ifftRec(z1Array, n/2);
            Complex[] z2 = ifftRec(z2Array, n/2);
            Complex omega = new Complex(Math.cos(2*Math.PI/n), Math.sin(2*Math.PI/n));
            for (int j = 0; j < m; j++) {
                v[j] = z1[j].add(omega.power(j).mul(z2[j]));
                v[m+j] = z1[j].sub(omega.power(j).mul(z2[j]));
            }
        }
        return v;
    }
}
