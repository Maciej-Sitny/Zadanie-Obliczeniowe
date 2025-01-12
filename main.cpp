#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
#include <fstream>

double E(double x)
{
    return (x <= 1.0) ? 3.0 : 5.0;
}

double e(int n, int i, double x)
{
    double h = 2.0 / n;
    double xi = i * h;
    double xiMinus1 = (i - 1) * h;
    double xiPlus1 = (i + 1) * h;

    if (xiMinus1 <= x && x < xi)
    {
        return (x - xiMinus1) / (xi - xiMinus1);
    }
    else if (xi <= x && x <= xiPlus1)
    {
        return (xiPlus1 - x) / (xiPlus1 - xi);
    }
    else
    {
        return 0.0;
    }
}

double ePrim(int n, int i, double x)
{
    double h = 2.0 / n;
    double xi = i * h;
    double xiMinus1 = (i - 1) * h;
    double xiPlus1 = (i + 1) * h;

    if (xiMinus1 <= x && x < xi)
    {
        return 1.0 / (xi - xiMinus1);
    }
    else if (xi <= x && x <= xiPlus1)
    {
        return -1.0 / (xiPlus1 - xi);
    }
    else
    {
        return 0.0;
    }
}

double gauss(double a, double b, std::function<double(double)> interior)
{
    const double sqrt3 = std::sqrt(3.0) / 3.0;
    const double x1 = -sqrt3, x2 = sqrt3;
    const double w1 = 1.0, w2 = 1.0;

    double xMapped1 = (b - a) / 2.0 * x1 + (a + b) / 2.0;
    double xMapped2 = (b - a) / 2.0 * x2 + (a + b) / 2.0; // przeskalowanie punktow z [-1,1] do [0,2]

    double result = w1 * interior(xMapped1) + w2 * interior(xMapped2);
    return result * (b - a) / 2.0;
}

double computeIntegral(int n, int i, int j)
{
    if (std::abs(i - j) > 1)
    {
        return 0.0;
    }

    double h = 2.0 / n;
    double a = std::max(0.0, (std::max(i, j) - 1) * h);
    double b = std::min(2.0, (std::min(i, j) + 1) * h); // granice calkowania

    auto interior = [n, i, j](double x) // srodek calki
    {
        return E(x) * ePrim(n, i, x) * ePrim(n, j, x);
    };

    return gauss(a, b, interior);
}

void createMatrix(int n, std::vector<std::vector<double>> &B, std::vector<double> &L)
{
    for (int j = 0; j < n; ++j)
    {
        L[j] = -30.0 * e(n, j, 0.0);
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            B[i][j] = computeIntegral(n, i, j) - 3.0 * e(n, i, 0.0) * e(n, j, 0.0);
        }
    }
}

std::vector<double> solveMatrix(std::vector<std::vector<double>> &B, std::vector<double> &L) // metoda eliminacji Gaussa
{
    int n = L.size();
    std::vector<double> x(n, 0.0);

    for (int i = 0; i < n; ++i)
    {
        for (int k = i + 1; k < n; ++k)
        {
            double factor = B[k][i] / B[i][i];
            for (int j = 0; j < n; ++j)
            {
                B[k][j] -= factor * B[i][j];
            }
            L[k] -= factor * L[i];
        }
    }

    for (int i = n - 1; i >= 0; --i)
    {
        x[i] = L[i];
        for (int j = i + 1; j < n; ++j)
        {
            x[i] -= B[i][j] * x[j];
        }
        x[i] /= B[i][i];
    }

    return x;
}

void toFile(int n, const std::vector<double> &uVector, const std::string &filename)
{
    double h = 2.0 / n;
    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Nie można otworzyć pliku: " << filename << "\n";
        return;
    }
    for (int i = 0; i <= n; ++i)
    {
        double x = i * h;
        double uX = 0.0;
        for (int j = 0; j < n; ++j)
        {
            uX += uVector[j] * e(n, j, x);
        }
        file << x << " " << uX << "\n";
    }
    file.close();
    std::cout << "Dane zapisane do pliku: " << filename << "\n";
}

int main()
{
    int n;
    std::cout << "Podaj liczbe elementow n: ";
    std::cin >> n;

    std::vector<std::vector<double>> B(n, std::vector<double>(n, 0.0));
    std::vector<double> L(n, 0.0);

    createMatrix(n, B, L);

    std::vector<double> u_vector = solveMatrix(B, L);

    toFile(n, u_vector, "data.txt");
    return 0;
}