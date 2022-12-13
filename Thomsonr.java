///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Program file name: Thomsonr.java                                      //
//                                                                       //
// © Tao Pang 2006                                                       //
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
// identical charges confined on a unit sphere by applying the real version
// of the Genetic Algorithm.

import java.lang.*;
import java.io.*;
import java.util.*;
public class Thomsonr {
  final static int nc = 60; // number of charges
  final static int nv = 2*nc; // number of variables
  final static int ng = 128; // size of gene pool
  final static int ni = 256; // initial size of gene pool
  final static int nr = ng/2; // number of parents or kids
  final static int ne = 1; // number of elite states
  final static int gmax = 60000; // max number of generations
  final static double pm = 1.0/100; // mutation percentage
  static double c[][] = new double[ng][nv]; // chromosomes
  static double f[] = new double[ng]; // cost function

  public static void main(String argv[]) throws
    FileNotFoundException {
    initiate();
    int g = 0;
    boolean contd = true;
    while (contd && g<gmax){
      ++g;
      select();
      cross();
      rank();
      mutate();
//    contd = converge();
    }
    export();
  }

// Method to initialize the simulation by creating an
// zeroth generation of the gene population.

  public static void initiate(){
    Random rnd = new Random();
    double d[][] = new double[ni][nv];
    double e[] = new double[ni];
    int index[] = new int[ni];

    for (int i=0; i<ni; ++i) {
      for (int j=0; j<nv; ++j) d[i][j] = rnd.nextDouble();
      e[i] = cost(d[i]);
      index[i] = i;
    }
    sort(e,index);

    for (int i=0; i<ng; ++i){
      f[i] = e[i];
      for (int j=0; j<nv; ++j) c[i][j] = d[index[i]][j];
    }
  }

// Method to run tournaments in selecting the parents.

  public static void select() {
    int index[] = new int[ng];
    double d[][] = new double[ng][nv];
    double e[] = new double[ng];

    for (int i=0; i<ng; ++i){
      for (int l=0; l<nv; ++l) d[i][l] = c[i][l];
      e[i] = f[i];
      index[i] = i;
    }
    shuffle(index);
    int k = 0;
    for (int i=0; i<nr; ++i) {
      if (e[index[k]] < e[index[k+1]]){
        for (int l=0; l<nv; ++l) c[i][l]=d[index[k]][l];
        f[i] = e[index[k]];
      }
      else {
        for (int l=0; l<nv; ++l) c[i][l]=d[index[k+1]][l];
        f[i] = e[index[k+1]];
      }
      k += 2;
    }
  }

// Method to perform single-point crossover operations.

  public static void cross() {
    Random rnd = new Random();

    int k = 0;
    for (int i=nr; i<nr+nr/2; ++i) {
      int nx = 1 + (int)(nv*rnd.nextDouble());
      for (int l=0; l<nx; ++l){
        c[i][l] = c[k][l];
        c[i+nr/2][l] = c[k+1][l];
      }
      for (int l=nx; l<nv; ++l){
        c[i][l] = c[k+1][l];
        c[i+nr/2][l] = c[k][l];
      }
      k += 2;
    }
  }

// Method to mutate a percentage of bits in the selected
// chromosomes except the best one.

  public static void mutate() {
    Random rnd = new Random();
    double r[] = new double[nv];

 // Mutation in the elite configurations
    for (int i=0; i<ne; ++i) {
      for (int l=0; l<nv; ++l) r[l] = c[i][l];
      int mb = (int)(nv*pm+1);
      for (int j=0; j<mb; ++j){
        int ib = (int)(nv*rnd.nextDouble());
        r[ib] = rnd.nextDouble();
      }
      double e = cost(r);
      if (e<f[i]){
        for (int l=0; l<nv; ++l) c[i][l] = r[l];
        f[i] = e;
      }
    }

 // Mutation in other configurations 
    int mmax = (int)((ng-ne)*nv*pm+1);
    for (int i=0; i<mmax; ++i) {
      int ig = (int)((ng-ne)*rnd.nextDouble()+ne);
      int ib = (int)(nv*rnd.nextDouble());
      c[ig][ib] = rnd.nextDouble();
    }

 // Rank the chromosomes in the population
    rank();
  }

// Method to save the coordinates and to output the cost function
// of the best chromosome in the final generation.

  public static void export() throws FileNotFoundException {
    double r[] = new double[nv];
    double theta[] = new double[nc];
    double phi[] = new double[nc];
    double x[] = new double[nc];
    double y[] = new double[nc];
    double z[] = new double[nc];
//  double rescale = Math.sqrt(nc*Math.E/8);
    double rescale = Math.sqrt(nc/4);

 // Write out the coordinates of the best chromosome
    PrintWriter q = new PrintWriter 
      (new FileOutputStream("xyz.data"), true);
    for (int i=0; i<nv; ++i) r[i] = c[0][i];
    int k = 0;
    for (int i=0; i<nc; ++i){
      theta[i] = Math.PI*r[k];
      phi[i] = 2*Math.PI*r[k+1];
      double ri = rescale*Math.sin(theta[i]);
      x[i] = ri*Math.cos(phi[i]);
      y[i] = ri*Math.sin(phi[i]);
      z[i] = rescale*Math.cos(theta[i]);
      q.println(x[i] + "  " + y[i] + "  " + z[i]);
      k += 2;
    }
 // Print out the energy of the best chromosome
    for (int i=0; i<ng; ++i)
      System.out.println("Energy: " + f[i]);
  }

// Method to rank chromosomes in the population.

  public static void rank() {
    double d[][] = new double[ng][nv];
    double e[] = new double[ng];
    int index[] = new int[ng];

    for (int i=0; i<ng; ++i) {
      for (int j=0; j<nv; ++j) d[i][j] = c[i][j];
      e[i] = cost(d[i]);
      index[i] = i;
    }
    sort(e,index);
    for (int i=0; i<ng; ++i){
      f[i] = e[i];
      for (int j=0; j<nv; ++j) c[i][j] = d[index[i]][j];
    }
  }

// Method to evaluate the cost for a given variable array.

  public static double cost(double r[]) {
    double g = 0;
    double theta[] = new double[nc];
    double phi[] = new double[nc];

    int k = 0;
    for (int i=0; i<nc; ++i){
      theta[i] = Math.PI*r[k];
      phi[i] = 2*Math.PI*r[k+1];
      k += 2;
    }

    for (int i=0; i<nc-1; ++i){
      double ri = Math.sin(theta[i]);
      double xi = ri*Math.cos(phi[i]);
      double yi = ri*Math.sin(phi[i]);
      double zi = Math.cos(theta[i]);
      for (int j=i+1; j<nc; ++j){
        double rj = Math.sin(theta[j]);
        double dx = xi - rj*Math.cos(phi[j]);
        double dy = yi - rj*Math.sin(phi[j]);
        double dz = zi - Math.cos(theta[j]);
        g += 1/Math.sqrt(dx*dx+dy*dy+dz*dz);
      }
    }
    return g;
  }

// Method to sort an array x[i] from the lowest to the
// highest with the original order stored in index[i].

  public static void sort(double x[], int index[]){
    int m = x.length;
    for (int i = 0; i<m; ++i) {
      for (int j = i+1; j<m; ++j) {
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

  public static void shuffle(int index[]){
    int k = index.length;
    Random rnd = new Random();
    for (int i = 0; i<k; ++i) {
      int j = (int)(k*rnd.nextDouble());
      if (j!=i) {
        int itmp = index[i];
        index[i] = index[j];
        index[j] = itmp;
      }
    }
  }
}
