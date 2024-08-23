#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <AMReX_REAL.H>

class CubicSplineInterpolator {
public:
	CubicSplineInterpolator(){}
    CubicSplineInterpolator(const std::string& filename, int x_column_index,int y_column_index) {
        define(filename, x_column_index,y_column_index);
    }

    void define(const std::string& filename,int x_column_index, int y_column_index) {
        // Read x and y values from the file
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file." << std::endl;
            std::exit(1);
        }

        // Read the file line by line
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            amrex::Real value;
            int column_index = 0;
            amrex::Real x = 0.0, y = 0.0;

            // Parse the line to extract x and y values
            while (iss >> value) {
                if (column_index == x_column_index) {
                    x = value;  // Store value as x (first column)
                }
                if (column_index == y_column_index) {
                    y = value;  // Store value as y
                }
                ++column_index;
            }
            x_values.push_back(x);
            y_values.push_back(y);
        }

        file.close();

        // Perform initial steps
        n = x_values.size();
        h.resize(n - 1);
        alpha.resize(n - 1);
        l.resize(n);
        mu.resize(n - 1);
        z.resize(n);
        c.resize(n);
        b.resize(n - 1);
        d.resize(n - 1);

        for (int i = 0; i < n - 1; ++i) {
            h[i] = x_values[i + 1] - x_values[i];
            alpha[i] = (3.0 / h[i]) * (y_values[i + 1] - y_values[i]) - (3.0 / h[std::max(i - 1, 0)]) * (y_values[i] - y_values[std::max(i - 1, 0)]);
        }

        l[0] = 1.0;
        mu[0] = 0.0;
        z[0] = 0.0;

        for (int i = 1; i < n - 1; ++i) {
            l[i] = 2.0 * (x_values[i + 1] - x_values[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n - 1] = 1.0;
        z[n - 1] = 0.0;
        c[n - 1] = 0.0;

        for (int j = n - 2; j >= 0; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y_values[j + 1] - y_values[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
    }

    amrex::Real interpolate(amrex::Real x_interpolate) const {
        if (x_interpolate <= x_values.front()) {
			//std::cerr << "Caution no datapoints for "<<x_interpolate<<" returning closest "<<y_values.front()<< std::endl; 
            return y_values.front();
        } else if (x_interpolate >= x_values.back()) {
			//std::cerr << "Caution no datapoints for "<<x_interpolate<<" returning closest "<<y_values.back()<< std::endl; 
            return y_values.back();
        }

        int i = findInterval(x_interpolate);
        amrex::Real dx = x_interpolate - x_values[i];
        
        amrex::Real ans = y_values[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        if(ans < y_values[i] || ans > y_values[i+1])
        {
           
           amrex::Real ans2 = y_values[i] + dx*(y_values[i+1] - y_values[i])/(x_values[i+1] - x_values[i]);
           if (ans2 < 1.0e-18L)
           {
           // std::cerr << "Caution overshooting, include more data points, returning the average value" << std::endl; 
           // std::cerr << "temp" <<x_interpolate<<"\t"<<"value"<< 0.5*(y_values[i] + y_values[i+1])<<std::endl; 
            return 0.5*(y_values[i] + y_values[i+1]) ;
           }
           else{
            //std::cerr << "Caution overshooting, include more data points, preforming linear interpolation." << std::endl; 
            //std::cerr << "temp" <<x_interpolate<<"\t"<<"value"<< ans2<<std::endl; 
            return ans2;
           }

        }
        else
        {
            return ans;
        }
    }

private:
    int findInterval(amrex::Real x_interpolate) const {
        auto it = std::lower_bound(x_values.begin(), x_values.end(), x_interpolate);
        return std::max(static_cast<int>(it - x_values.begin()) - 1, 0);
    }

private:
    int n;
    std::vector<amrex::Real> x_values;
    std::vector<amrex::Real> y_values;
    std::vector<amrex::Real> h;
    std::vector<amrex::Real> alpha;
    std::vector<amrex::Real> l;
    std::vector<amrex::Real> mu;
    std::vector<amrex::Real> z;
    std::vector<amrex::Real> c;
    std::vector<amrex::Real> b;
    std::vector<amrex::Real> d;
};

//int main() {
    //// Example usage
    //std::string filename = "data.txt";
    //CubicSplineInterpolator interpolator(filename,0, 11);

    //// Set precision for std::cout
    //std::cout << std::fixed << std::setprecision(std::numeric_limits<amrex::Real>::digits10 + 1);

    //// Interpolate at a specific x value
    //amrex::Real x_to_interpolate = 1456.6497845218;
    //amrex::Real result = interpolator.interpolate(x_to_interpolate);
    //std::cout << "Interpolated value at x = " << x_to_interpolate << ": " << result << "\n";

    //return 0;
//}
