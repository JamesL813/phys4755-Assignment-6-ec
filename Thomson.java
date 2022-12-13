///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Program file name: Thomson.java                                       //
//                                                                       //
// Tao Pang 2006                                                         //
//                                                                       //
// Last modified: January 18, 2006                                       //
//                                                                       //
// (1) This Java program is part of the book, "An Introduction to        //
//     Computational Physics, 2nd Edition," written by Tao Pang and      //
//     published by Cambridge University Press on January 19, 2006.      //
//                                                                       //
// (2) No warranties, express or implied, are made for this program.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// This program solves the Thomson problem, that is, to find the geometric
// configuration of the lowest electrostatic energy for a given number of
// identical charges confined on a unit sphere by applying the Genetic
// Algorithm.

import java.lang.*;
import java.io.*;
import java.util.*;

public class Thomson {
    final static int nc = 4; // number of charges
    final static int nv = 2 * nc; // number of variables
    final static int nb = 64; // number of bits used in encoding
    final static int nd = nv * nb;
    final static int ng = 128; // size of gene pool
    final static int ni = 256; // initial size of gene pool
    final static int nr = ng / 2; // number of parents or kids
    final static int ne = 1; // number of elite states
    // final static int gmax = 60000; // max number of generations
    final static int gmax = 600; // max number of generations
    final static double pm = 1.0 / 100; // mutation percentage
    static boolean c[][] = new boolean[ng][nd]; // chromosomes
    static double f[] = new double[ng]; // cost function

    public static void main(String argv[]) throws FileNotFoundException {
        initiate();
        int g = 0;
        boolean contd = true;
        // while (contd && g < gmax) {
        while (g < gmax) {
            System.out.println("g: " + g);
            ++g;
            select();
            cross();
            rank();
            mutate();
            // contd = converge();
        }
        export();
    }

    // Method to initialize the simulation by creating an
    // zeroth generation of the gene population.

    public static void initiate() {
        // System.out.println("Initiate");
        Random rnd = new Random();
        boolean d[][] = new boolean[ni][nd];
        boolean w[] = new boolean[nd];
        double r[] = new double[nv];
        double e[] = new double[ni];
        int index[] = new int[ni];

        for (int i = 0; i < ni; ++i) {
            for (int j = 0; j < nv; ++j)
                r[j] = rnd.nextDouble();
            e[i] = cost(r);
            index[i] = i;
            w = encode(r, nb);
            for (int j = 0; j < nd; ++j)
                d[i][j] = w[j];
        }
        sort(e, index);
        for (int i = 0; i < ng; ++i) {
            f[i] = e[i];
            for (int j = 0; j < nd; ++j)
                c[i][j] = d[index[i]][j];
        }
    }

    // Method to run tournaments in selecting the parents.

    public static void select() {
        // System.out.println("Select");
        int index[] = new int[ng];
        boolean d[][] = new boolean[ng][nd];
        double e[] = new double[ng];

        for (int i = 0; i < ng; ++i) {
            for (int l = 0; l < nd; ++l)
                d[i][l] = c[i][l];
            e[i] = f[i];
            index[i] = i;
        }
        shuffle(index);
        int k = 0;
        for (int i = 0; i < nr; ++i) {
            if (e[index[k]] < e[index[k + 1]]) {
                for (int l = 0; l < nd; ++l)
                    c[i][l] = d[index[k]][l];
                f[i] = e[index[k]];
            } else {
                for (int l = 0; l < nd; ++l)
                    c[i][l] = d[index[k + 1]][l];
                f[i] = e[index[k + 1]];
            }
            k += 2;
        }
    }

    // Method to perform single-point crossover operations.

    public static void cross() {
        // System.out.println("Cross");
        Random rnd = new Random();

        int k = 0;
        for (int i = nr; i < nr + nr / 2; ++i) {
            int nx = 1 + (int) (nd * rnd.nextDouble());
            for (int l = 0; l < nx; ++l) {
                c[i][l] = c[k][l];
                c[i + nr / 2][l] = c[k + 1][l];
            }
            for (int l = nx; l < nd; ++l) {
                c[i][l] = c[k + 1][l];
                c[i + nr / 2][l] = c[k][l];
            }
            k += 2;
        }
    }

    // Method to mutate a percentage of bits in the selected
    // chromosomes except the best one.

    public static void mutate() {
        // System.out.println("Mutate");
        Random rnd = new Random();
        double r[] = new double[nv];
        boolean w[] = new boolean[nd];

        // Mutation in the elite configurations
        for (int i = 0; i < ne; ++i) {
            for (int l = 0; l < nd; ++l)
                w[l] = c[i][l];
            int mb = (int) (nd * pm + 1);
            for (int j = 0; j < mb; ++j) {
                int ib = (int) (nd * rnd.nextDouble());
                w[ib] = !w[ib];
            }
            r = decode(w, nb);
            double e = cost(r);
            if (e < f[i]) {
                for (int l = 0; l < nd; ++l)
                    c[i][l] = w[l];
                f[i] = e;
            }
        }

        // Mutation in other configurations
        int mmax = (int) ((ng - ne) * nd * pm + 1);
        for (int i = 0; i < mmax; ++i) {
            int ig = (int) ((ng - ne) * rnd.nextDouble() + ne);
            int ib = (int) (nd * rnd.nextDouble());
            c[ig][ib] = !c[ig][ib];
        }

        // Rank the chromosomes in the population
        rank();
    }

    // Method to save the coordinates and to output the cost function
    // of the best chromosome in the final generation.

    public static void export() throws FileNotFoundException {
        // System.out.println("Export");
        double r[] = new double[nv];
        boolean w[] = new boolean[nd];
        double theta[] = new double[nc];
        double phi[] = new double[nc];
        double x[] = new double[nc];
        double y[] = new double[nc];
        double z[] = new double[nc];
        // double rescale = Math.sqrt(nc*Math.E/8);
        double rescale = Math.sqrt(nc / 4.0);
        // double rescale = 2.5;

        // Write out the coordinates of the best chromosome
        PrintWriter q = new PrintWriter(new FileOutputStream("xyz.data"), true);
        q.println("x y z");
        for (int i = 0; i < nd; ++i)
            w[i] = c[0][i];
        r = decode(w, nb);
        int k = 0;
        for (int i = 0; i < nc; ++i) {
            theta[i] = Math.PI * r[k];
            phi[i] = 2 * Math.PI * r[k + 1];
            double ri = rescale * Math.sin(theta[i]);
            x[i] = ri * Math.cos(phi[i]);
            y[i] = ri * Math.sin(phi[i]);
            z[i] = rescale * Math.cos(theta[i]);
            q.println(x[i] + " " + y[i] + " " + z[i]);
            k += 2;
        }
        // Print out the energy of the best chromosome
        for (int i = 0; i < ng; ++i)
            System.out.println("Energy: " + f[i]);
    }

    // Method to rank chromosomes in the population.

    public static void rank() {
        // System.out.println("Rank");
        boolean d[][] = new boolean[ng][nd];
        double r[] = new double[nv];
        double e[] = new double[ng];
        int index[] = new int[ng];

        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < nd; ++j)
                d[i][j] = c[i][j];
            r = decode(d[i], nb);
            e[i] = cost(r);
            index[i] = i;
        }
        sort(e, index);
        for (int i = 0; i < ng; ++i) {
            f[i] = e[i];
            for (int j = 0; j < nd; ++j)
                c[i][j] = d[index[i]][j];
        }
    }

    // Method to evaluate the cost for a given variable array.

    public static double cost(double r[]) {
        // System.out.println("Cost");
        double g = 0;
        double theta[] = new double[nc];
        double phi[] = new double[nc];

        int k = 0;
        for (int i = 0; i < nc; ++i) {
            theta[i] = Math.PI * r[k];
            phi[i] = 2 * Math.PI * r[k + 1];
            k += 2;
        }

        for (int i = 0; i < nc - 1; ++i) {
            double ri = Math.sin(theta[i]);
            double xi = ri * Math.cos(phi[i]);
            double yi = ri * Math.sin(phi[i]);
            double zi = Math.cos(theta[i]);
            double minDist = ri;
            int minCount = 0;
            for (int j = i + 1; j < nc; ++j) {
                double rj = Math.sin(theta[j]);
                double dx = xi - rj * Math.cos(phi[j]);
                double dy = yi - rj * Math.sin(phi[j]);
                double dz = zi - Math.cos(theta[j]);

                double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);

                if (Math.abs(dist - minDist) < 0.0001) {
                    minCount++;
                } else if (dist < minDist) {
                    minDist = dist;
                    minCount = 0;
                }

                g += 1 / dist;
            }
            g -= minCount;
        }
        return g;
    }

    // Method to encode an array of n real numbers r[i] in
    // [0,1] into an n*m binary representation w[j].

    public static boolean[] encode(double r[], int m) {
        // System.out.println("Encode");
        int n = r.length;
        boolean w[] = new boolean[n * m];
        for (int i = 0; i < n; ++i) {
            double sum = r[i];
            w[i * m] = false;
            if ((int) (0.5 + sum) == 1)
                w[i * m] = true;
            double d = 2;
            for (int j = 1; j < m; ++j) {
                if (w[i * m + j - 1])
                    sum -= 1 / d;
                w[i * m + j] = false;
                if ((int) (0.5 + d * sum) == 1)
                    w[i * m + j] = true;
                d *= 2;
            }
        }
        return w;
    }

    // Method to decode an array of n*m binary numbers w[j]
    // into an array of n real numbers r[i].

    public static double[] decode(boolean w[], int m) {
        // System.out.println("Decode");
        int n = w.length / m;
        double r[] = new double[n];
        for (int i = 0; i < n; ++i) {
            double d = 2;
            double sum = 0;
            for (int j = 0; j < m; ++j) {
                if (w[i * m + j])
                    sum += 1 / d;
                d *= 2;
            }
            r[i] = sum + 1 / d;
        }
        return r;
    }

    // Method to sort an array x[i] from the lowest to the
    // highest with the original order stored in index[i].

    public static void sort(double x[], int index[]) {
        // System.out.println("Sort");
        int m = x.length;
        for (int i = 0; i < m; ++i) {
            for (int j = i + 1; j < m; ++j) {
                if (x[i] > x[j]) {
                    double xtmp = x[i];
                    x[i] = x[j];
                    x[j] = xtmp;
                    int itmp = index[i];
                    index[i] = index[j];
                    index[j] = itmp;
                }
            }
        }
    }

    // Method to shuffle the index array.

    public static void shuffle(int index[]) {
        // System.out.println("Shuffle");
        int k = index.length;
        Random rnd = new Random();
        for (int i = 0; i < k; ++i) {
            int j = (int) (k * rnd.nextDouble());
            if (j != i) {
                int itmp = index[i];
                index[i] = index[j];
                index[j] = itmp;
            }
        }
    }

    // Method to check the convergence.

    // public static boolean converge(){ }

}
