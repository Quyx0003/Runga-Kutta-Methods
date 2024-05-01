#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

std::vector<double> operator*(double scalar, const std::vector<double> &v) {
    std::vector<double> result(v.size());
    for (int i = 0; i < v.size(); i++) {
        result[i] = v[i] * scalar;
    }
    return result;
}

std::vector<double> operator+(const std::vector<double> &v1, const std::vector<double> &v2) {
    std::vector<double> result(v1.size());
    for (int i = 0; i < v1.size(); i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

std::vector<double> F(double x, const std::vector<double> &y) {
    std::vector<double> f(2);
    f[0] = std::cos(-1 + x + y[0] + 3 * y[1]);
    f[1] = -y[0] * y[0] + 2 * std::sin(y[1]);
    return f;
}

void Euler(double n, std::vector<double> (*func)(double, const std::vector<double> &), std::vector<double> &a, double x0, double xEnd, double accuracy, int eq) {
    std::vector<std::vector<double>> container;
    std::cout << std::setw(2)  << "i"; 
    std::cout << std::setw(15) << "A(h_i)";
    std::cout << std::setw(20) << "A(h_(i-1))-A(h_i)"; 
    std::cout << std::setw(15) << "alpha_k";
    std::cout << std::setw(15) << "rich-error"; 
    std::cout << std::setw(15) << "order";
    std::cout << std::setw(10) << "n";
    std::cout << std::endl;

    double prev_prev_integral = 0.0;
    double prev_integral = 0.0;
    for (int i = 0; i < 50; i++) {
        std::vector<double> yn = a;
        double h = (xEnd - x0) / n;
        for (int j = 0; j < n; j++) {
            yn = yn + h * func(x0 + j * h, yn);
        }
        
        container.push_back(yn);

        double integral = container[i][eq];
        std::cout << std::setw(2) << i << std::setw(15) << integral;
        
        if (i > 0) {
            double diff = prev_integral - integral;
            std::cout << std::setw(20) << diff;
            if (i > 1) {
                double diff1 = prev_prev_integral - prev_integral; 
                double alpha_k = diff1 / diff;           
                double richardson = (diff) / (alpha_k - 1.0); 
                std::cout << std::setw(15) << alpha_k; 
                std::cout << std::setw(15) << richardson;
                std::cout << std::setw(15) << std::log2(diff1 / diff);
                if (std::abs(richardson) < accuracy) {
                    std::cout << std::setw(10) << n << std::endl;
                    std::cout << std::endl;
                    break;
                }
            } else {
                std::cout << std::setw(45) << " ";
            }
        } else {
            std::cout << std::setw(65) << " ";
        }
        std::cout << std::setw(10) << n;
        std::cout << std::endl;
        prev_prev_integral = prev_integral;
        prev_integral = integral;
        n *= 2;
    }
    std::cout << std::endl;
}

void Midpoint(double n, std::vector<double> (*func)(double, const std::vector<double> &), std::vector<double> &a, double x0, double xEnd, double accuracy, int eq) {
    std::vector<std::vector<double>> container;
    std::cout << std::setw(2)  << "i"; 
    std::cout << std::setw(15) << "A(h_i)";
    std::cout << std::setw(20) << "A(h_(i-1))-A(h_i)"; 
    std::cout << std::setw(15) << "alpha_k";
    std::cout << std::setw(15) << "rich-error"; 
    std::cout << std::setw(15) << "order";
    std::cout << std::setw(10) << "n";
    std::cout << std::endl;

    double prev_prev_integral = 0.0;
    double prev_integral = 0.0;
    for (int i = 0; i < 50; i++) {
        std::vector<double> yn = a;
        double h = (xEnd - x0) / n;
        for (int j = 0; j < n; j++) {
            std::vector<double> k1 = func(x0 + j * h, yn);
            std::vector<double> k2 = func(x0 + j * h + h * 0.5, yn + 0.5 * h * k1);
            yn = yn + h * k2;
        }
        
        container.push_back(yn);

        double integral = container[i][eq];
        std::cout << std::setw(2) << i << std::setw(15) << integral;
        
        if (i > 0) {
            double diff = prev_integral - integral;
            std::cout << std::setw(20) << diff;
            if (i > 1) {
                double diff1 = prev_prev_integral - prev_integral; 
                double alpha_k = diff1 / diff;           
                double richardson = (diff) / (alpha_k - 1.0); 
                std::cout << std::setw(15) << alpha_k; 
                std::cout << std::setw(15) << richardson;
                std::cout << std::setw(15) << std::log2(diff1 / diff);
                if (std::abs(richardson) < accuracy) {
                    std::cout << std::setw(10) << n << std::endl;
                    std::cout << std::endl;
                    break;
                }
            } else {
                std::cout << std::setw(45) << " ";
            }
        } else {
            std::cout << std::setw(65) << " ";
        }
        std::cout << std::setw(10) << n;
        std::cout << std::endl;
        prev_prev_integral = prev_integral;
        prev_integral = integral;
        n *= 2;
    }
    std::cout << std::endl;
}

void FourthOrder(double n, std::vector<double> (*func)(double, const std::vector<double> &), std::vector<double> &a, double x0, double xEnd, double accuracy, int eq) {
    std::vector<std::vector<double>> container;
    std::cout << std::setw(2)  << "i"; 
    std::cout << std::setw(15) << "A(h_i)";
    std::cout << std::setw(20) << "A(h_(i-1))-A(h_i)"; 
    std::cout << std::setw(15) << "alpha_k";
    std::cout << std::setw(15) << "rich-error"; 
    std::cout << std::setw(15) << "order";
    std::cout << std::setw(10) << "n";
    std::cout << std::endl;

    double prev_prev_integral = 0.0;
    double prev_integral = 0.0;
    for (int i = 0; i < 50; i++) {
        std::vector<double> yn = a;
        double h = (xEnd - x0) / n;
        for (int j = 0; j < n; j++) {
            std::vector<double> k1 = func(x0 + j * h, yn);
            std::vector<double> k2 = func(x0 + j * h + h * 0.5, yn + 0.5 * h * k1);
            std::vector<double> k3 = func(x0 + j * h + h * 0.5, yn + 0.5 * h * k2);
            std::vector<double> k4 = func(x0 + j * h + h, yn + h * k3);
            yn = yn + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        
        container.push_back(yn);

        double integral = container[i][eq];
        std::cout << std::setw(2) << i << std::setw(15) << integral;
        
        if (i > 0) {
            double diff = prev_integral - integral;
            std::cout << std::setw(20) << diff;
            if (i > 1) {
                double diff1 = prev_prev_integral - prev_integral; 
                double alpha_k = diff1 / diff;           
                double richardson = (diff) / (alpha_k - 1.0); 
                std::cout << std::setw(15) << alpha_k; 
                std::cout << std::setw(15) << richardson;
                std::cout << std::setw(15) << std::log2(diff1 / diff);
                if (std::abs(richardson) < accuracy) {
                    std::cout << std::setw(10) << n << std::endl;
                    std::cout << std::endl;
                    break;
                }
            } else {
                std::cout << std::setw(45) << " ";
            }
        } else {
            std::cout << std::setw(65) << " ";
        }
        std::cout << std::setw(10) << n;
        std::cout << std::endl;
        prev_prev_integral = prev_integral;
        prev_integral = integral;
        n *= 2;
    }
    std::cout << std::endl;
}

void Trapezoidal(double n, std::vector<double> (*func)(double, const std::vector<double> &), std::vector<double> &a, double x0, double xEnd, double accuracy, int eq) {
    std::vector<std::vector<double>> container;
    std::cout << std::setw(2)  << "i"; 
    std::cout << std::setw(15) << "A(h_i)";
    std::cout << std::setw(20) << "A(h_(i-1))-A(h_i)"; 
    std::cout << std::setw(15) << "alpha_k";
    std::cout << std::setw(15) << "rich-error"; 
    std::cout << std::setw(15) << "order";
    std::cout << std::setw(10) << "n";
    std::cout << std::endl;

    double prev_prev_integral = 0.0;
    double prev_integral = 0.0;
    for (int i = 0; i < 50; i++) {
        std::vector<double> yn = a;
        double h = (xEnd - x0) / n;
        for (int j = 0; j < n; j++) {
            std::vector<double> f_val = func(x0 + j * h, yn);
			std::vector<double> f_val_next = func(x0 + (j + 1) * h, yn + h * f_val);
			yn = yn + h * 0.5 * (f_val + f_val_next);
        }
        
        container.push_back(yn);

        double integral = container[i][eq];
        std::cout << std::setw(2) << i << std::setw(15) << integral;
        
        if (i > 0) {
            double diff = prev_integral - integral;
            std::cout << std::setw(20) << diff;
            if (i > 1) {
                double diff1 = prev_prev_integral - prev_integral; 
                double alpha_k = diff1 / diff;           
                double richardson = (diff) / (alpha_k - 1.0); 
                std::cout << std::setw(15) << alpha_k; 
                std::cout << std::setw(15) << richardson;
                std::cout << std::setw(15) << std::log2(diff1 / diff);
                if (std::abs(richardson) < accuracy) {
                    std::cout << std::setw(10) << n << std::endl;
                    std::cout << std::endl;
                    break;
                }
            } else {
                std::cout << std::setw(45) << " ";
            }
        } else {
            std::cout << std::setw(65) << " ";
        }
        std::cout << std::setw(10) << n;
        std::cout << std::endl;
        prev_prev_integral = prev_integral;
        prev_integral = integral;
        n *= 2;
    }
    std::cout << std::endl;
}

int main() {
    double x0 = 0.0;
    double xEnd = 1.0;
    std::vector<double> a(2);
    a[0] = 1;
    a[1] = 0;
    double n = 1.0;
    double accuracy = 1e-6;

    std::cout << "Euler Method results for u(x):" << std::endl;
    Euler(n, F, a, x0, xEnd, accuracy, 0);

    std::cout << "Euler Method results for v(x):" << std::endl;
    Euler(n, F, a, x0, xEnd, accuracy, 1);

    std::cout << "Midpoint Method results for u(x):" << std::endl;
    Midpoint(n, F, a, x0, xEnd, accuracy, 0);

    std::cout << "Midpoint Method results for v(x):" << std::endl;
    Midpoint(n, F, a, x0, xEnd, accuracy, 1);

    std::cout << "Fourth Order Method results for u(x):" << std::endl;
    FourthOrder(n, F, a, x0, xEnd, accuracy, 0);

    std::cout << "Fourth Order Method results for v(x):" << std::endl;
    FourthOrder(n, F, a, x0, xEnd, accuracy, 1);

    std::cout << "Trapezoidal Method results for u(x):" << std::endl;
    Trapezoidal(n, F, a, x0, xEnd, accuracy, 0);

    std::cout << "Trapezoidal Method results for v(x):" << std::endl;
    Trapezoidal(n, F, a, x0, xEnd, accuracy, 1);

    return 0;
}


