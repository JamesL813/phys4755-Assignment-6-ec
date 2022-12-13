---
geometry:
  - top=10mm
  - left=10mm
  - right=10mm
  - bottom=10mm
---

# Problem 1

Consider the Poisson equation
$$\nabla^2 \phi(x, y) = -\rho(x, y)/\epsilon_0$$

from electrostatics on a rectangular geometry with $x \in [0, L_x ]$ and $y \in
[0, L_y ]$. Write a program that solves this equation using the relaxation
method. Test your program with:

### (a)

$\rho(x, y) = 0$

$\phi(0, y) = \phi(L_x , y) = \phi(x, 0) = 0$

$\phi(x, L_y ) = 1 V, L_x = 1 m$

$L_y = 1.5 m$

### (b)

$\rho(x, y)/\epsilon_0 = 1 V/m^2$

$\phi(0, y) = \phi(L_x , y) = \phi(x, 0) = \phi(x, L_y ) = 0$

$L_x = L_y = 1 m$

# Problem 2

Consider a nuclear waste rod, similar to the situation described in section 7.8, but let's model it along a one dimensional line instead of in cylindrical coordinates. The relevant diffusion equation for temperature becomes:

$$
    \frac{1}{\kappa}
    \frac{\partial T(x,t)}{\partial t} -
    \frac{\partial^2 T(x,t)}{\partial^x x} =
    S(x,t)
$$

Let's use different boundary conditions than what's discussed in section 7.8: assume the rod is still 100 cm long but has the left and right ends held at the temperature of the earth (T_E=300K). Also, let's use a different source function. Let the source function add heat to just the middle 30 cm of the rod :

$S(x,t) = \frac{T_0}{a^2}e^{-t/\tau_0}$, for $x > 35$cm, and $x <  65$cm, but $S(x,t)=0$ otherwise

where T_0=1K, $\tau_0 = 100$ years, and $a = 35$cm, similar to parameters from section 7.8. See the figure below:

![image](./p2/plot_rod.png)

Assume that the temperature of the whole rod starts at TE at time t=0. The goal is to solve the diffusion equation to find the temperature difference $$ as a function of position ( for x between 0 and 100cm) and time. See Figure 7.3 in the text as an example way to plot results. Let's solve the problem using two methods:

Method 1-Use the Crank-Nicolson method to find the temperature.  You can write your own Crank-Nicolson code or use the one in section 7.8. If you're going to use the template code in Section 7.8, note that the diffusion equation above and 7.102 are almost the same, but not quite. You'll have to make changes.

Method 2-Write a program to solve for the temperature using the Euler method, Equation 7.90 . The Euler method should start to give wrong answers for large.

For your plots make sure to choose $\kappa$ to illustrate two scenarios as you compare methods: small $\kappa$, where both methods nearly agree; and large $\kappa$, where they start to disagree. Recall from the text that a realistic value should be around $\kappa=2.00*10^7cm^2/100yr$. In your write up, discuss both methods, describe your results, and be sure to discuss why and where the methods start to disagree.

# Problem 3

Do problem 10.3 from your text. Plot the answer (with error bars) as a function of the number of Monte Carlo steps. Make sure to vary the number of steps so that you can see the trend where the Monte Carlo becomes more and more accurate as you increase the number of steps.

### 10.3

Calculate the integral

$$
S=
\int_{-\infty}^{\infty}
e^{-r^2/2} (xyz)^2 dr
$$

with the Metropolis algorithm and compare the Monte Carlo result with the
exact result. Does the Monte Carlo errorbar decrease with the total number
of points as expected?

# Extra Credit (worth 50 points added to this problem set):

In this project you'll teach yourself the basics of genetic algorithms that use the laws of evolution to solve problems. Before starting, please study Chapter 11 in the text. Your project is to use a genetic algorithm to find several of the Platonic solids. Include schematic plots of the Platonic solids that you found. Plot the distance between vertices for the 4 vertex case as a function of the maximum number of generations. Hints: There are only 5 Platonic solids. Also, solutions of the Thompson problem of like charges sitting on a sphere can be used to search for some of the Platonic solids.
