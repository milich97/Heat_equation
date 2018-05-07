import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.io.IOException;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;

/**
 * Created by Миша on 06.05.2018.
 */
/*програма соответствует теории, описанной на с. 205-210 в http://www.math.spbu.ru/user/pan/site_full.pdf.
   Применение дискретного ряда Фурье для решения методом сеток однородного уравнения теплопроводности
 */
public class Main {
    private static DoubleUnaryOperator phi = x -> Math.sin(Math.PI * x * (x - 1));
//    private static DoubleUnaryOperator phi = x -> (x - 1) * Math.sin(Math.PI * x / 3);
//    private static DoubleUnaryOperator phi = x -> Math.sin(Math.PI * x);

    final private static double pi = Math.PI;
    private static double T = 0.1;
    private static DoubleBinaryOperator u_acc;
    private static double[] c;

    public static void main(String[] args) throws IOException {

        int n_acc = 300;
        u_acc = get_u_acc(n_acc);
//        u_acc=get_ud(20);
        System.out.print("1)");
        table(u_acc);
//        table((x,t)->Math.exp(-pi*pi*t)*Math.sin(pi*x));
        System.out.println();
        System.out.print("2)");
        DoubleBinaryOperator ud = get_ud(5);
        System.out.println();
        System.out.println("||uf-udsf(5)|| = " + max(ud));
        ud = get_ud(10);
        System.out.println("||uf-udsf(10)|| = " + max(ud));
        ud = get_ud(20);
        System.out.println("||uf-udsf(20)|| = " + max(ud));
        System.out.println();
        System.out.println("3)");
        double sigma = 0;
        printRow(sigma);
        sigma = 1;
        printRow(sigma);
        sigma = 0.5;
        printRow(sigma);

        double h = 0.2, tau = 0.02;
        sigma = 0.5 - h * h / (12 * tau);
        System.out.print(norm(sigma, h, tau) + "    ");
        h = 0.1;
        tau = 0.005;
        sigma = 0.5 - h * h / (12 * tau);
        System.out.print(norm(sigma, h, tau) + "    ");
        h = 0.05;
        tau = 0.00125;
        sigma = 0.5 - h * h / (12 * tau);
        System.out.print(norm(sigma, h, tau) + "    ");
        h = 0.05;
        tau = 0.01;
        sigma = 0.5 - h * h / (12 * tau);
        System.out.print(norm(sigma, h, tau) + "    ");

    }

    private static void printRow(double sigma) {
        double h = 0.2, tau = 0.02;
        System.out.print(norm(sigma, h, tau) + "    ");
        h = 0.1;
        tau = 0.005;
        System.out.print(norm(sigma, h, tau) + "    ");
        h = 0.05;
        tau = 0.00125;
        System.out.print(norm(sigma, h, tau) + "    ");
        h = 0.05;
        tau = 0.01;
        System.out.print(norm(sigma, h, tau) + "    ");
        System.out.println();
    }

    private static double norm(double sigma, double h, double tau) {
        double[][] u_grid = new double[6][6];
        int N = (int) (1 / h);
        for (int ii = 0; ii < 6; ii++) {
            for (int kk = 0; kk < 6; kk++) {
                int i = (int) (0.2 * ii * N);
                int k = (int) (0.2 * kk * T / tau);
                double term = 0;
                for (int p = 1; p < N; p++) {
                    double mu = -4 * N * N * Math.pow(Math.sin(p * pi * h / 2), 2);
                    double lambda = (1 + tau * (1 - sigma) * mu) / (1 - tau * sigma * mu);
                    term += c[p - 1] * Math.pow(lambda, k) * Math.sin(p * pi * i * h);
                }
                u_grid[ii][kk] = term * Math.sqrt(2);
            }
        }

        double ans = 0;
        double max;
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                max = Math.abs(u_acc.applyAsDouble(j * 0.2, i * 0.2 * T) - u_grid[j][i]);
                if (max > ans) ans = max;
            }
        }
//        if (sigma>=0.5) System.out.println("stability");
//        else if (tau<h*h/2/(1-2*sigma)) System.out.println("stability");
//        else System.out.println("no stability!!!");
        return ans;
    }

    private static void table(DoubleBinaryOperator u) {
        System.out.println();
        System.out.println("t/x    0   0.2   0.4   0.6   0.8   1");
        for (int i = 0; i < 6; i++) {
            switch (i) {
                case 0:
                    System.out.print("0      ");
                    break;
                case 1:
                    System.out.print("T/5    ");
                    break;
                case 2:
                    System.out.print("2T/5   ");
                    break;
                case 3:
                    System.out.print("3T/5   ");
                    break;
                case 4:
                    System.out.print("4T/5   ");
                    break;
                case 5:
                    System.out.print("T      ");
                    break;
            }
            for (int j = 0; j < 6; j++) {
                System.out.print((u.applyAsDouble(j * 0.2, i * 0.2 * T)) + "   ");
            }
            System.out.println();
        }
    }


    private static DoubleBinaryOperator get_ud(int N) {
        double[] c = new double[N - 1];
        double h = 1.0 / N;
        for (int k = 1; k < N; k++) {
            c[k - 1] = 0;
            for (int i = 1; i < N; i++) {
                c[k - 1] += phi.applyAsDouble(i * h) * Math.sin(k * pi * i * h);
            }
            c[k - 1] *= h * Math.sqrt(2);
        }
        return get_u(c);
    }

    public static DoubleBinaryOperator get_u_acc(int n_acc) throws IOException {

        c = new double[n_acc];
        for (int i = 1; i < n_acc + 1; i++) {
            IterativeLegendreGaussIntegrator I =
                    new IterativeLegendreGaussIntegrator(10, 1, 100);
            int finalI = i;
            UnivariateFunction g = x -> Math.sin(pi * finalI * x) * phi.applyAsDouble(x);
            c[i - 1] = Math.sqrt(2) * I.integrate(Integer.MAX_VALUE, g, 0, 1);
        }


        return get_u(c);
    }

    private static double r(double value, int k) {                 // функция для округления
        return (double) Math.round((Math.pow(10, k) * value)) / Math.pow(10, k);
    }

    private static double max(DoubleBinaryOperator ud) {
        double ans = 0;
        double help;
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                help = Math.abs(ud.applyAsDouble(j * 0.2, i * 0.2 * T) - u_acc.applyAsDouble(j * 0.2, i * 0.2 * T));
                if (help > ans) ans = help;
            }
        }
        return ans;
    }

    private static DoubleBinaryOperator get_u(double[] c) {
        DoubleBinaryOperator help = (x, t) -> c[0] * Math.exp(-pi * pi * t) * Math.sin(pi * x);
        for (int i = 2; i < c.length + 1; i++) {
            int p = i;
            DoubleBinaryOperator finalHelp = help;
            help = (x, t) -> finalHelp.applyAsDouble(x, t) +
                    c[p - 1] * Math.exp(-pi * pi * p * p * t) * Math.sin(pi * p * x);
        }
        DoubleBinaryOperator finalHelp1 = help;

        return (x, t) -> finalHelp1.applyAsDouble(x, t) * Math.sqrt(2);
    }
}
