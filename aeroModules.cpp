#include <iostream>
#include <cmath>
#include <string>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <sstream>
// #include "xtensor/xarray.hpp"
// #include "xtensor/xio.hpp"
// #include "xtensor/xview.hpp"

using namespace std;

#define PI 3.141592653589793238463

tuple<float, float> invRotation(float x, float y, float angle)
{
    float xt =  x*cos(angle) + y*sin(angle);
    float yt = -x*sin(angle) + y*cos(angle);
    
    return make_tuple(xt, yt);
}

tuple<float, float> rotation(float x, float y, float angle)
{
    float xt = x*cos(angle) - y*sin(angle);
    float yt = x*sin(angle) + y*cos(angle);
    
    return make_tuple(xt, yt);
}

tuple<float, float> panelCoords(float x, float y, float x0, float y0, float angle)
{
    return rotation(x - x0, y - y0, angle);
}

class Vortex2D
{
    public:
        float strength, x0, y0;
    
    Vortex2D(float gamma, float x_0, float y_0)
    {
        strength = gamma, x0 = x_0, y0 = y_0;
    }

    tuple<float, float> velocity(float x, float y)
    {
        float u = -strength/(2*PI)*(y - y0)/(pow(x - x0, 2) + pow(y - y0, 2));
        float v = strength/(2*PI)*(x - x0)/(pow(x - x0, 2) + pow(y - y0, 2));
        
        return make_tuple(u, v);
    }
};

class Doublet2D
{
    public:
        float strength, x0, y0;
    
    Doublet2D(float mu, float x_0, float y_0)
    {
        strength = mu, x0 = x_0, y0 = y_0;
    }

    tuple<float, float> velocity(float x, float y)
    {
        float u = -strength/(2*PI)*(x - x0)/pow((pow(x - x0, 2) + pow(y - y0, 2)), 2);
        float v = -strength/(2*PI)*(y - y0)/pow((pow(x - x0, 2) + pow(y - y0, 2)), 2);
        
        return make_tuple(u, v);
    }
};

class Uniform2D
{
    public:
        float magnitude = 1.0, angle = 0.0;
        tuple<float, float> velocity();
    
    Uniform2D() {}
    Uniform2D(float mag, float ang)
    {
        magnitude = mag, angle = ang*PI/180.0;
    }

    tuple<float, float> velocity(float x, float y)
    {
        float u = magnitude*cos(angle);
        float v = magnitude*sin(angle);
        
        return make_tuple(u, v);
    }
};

class Panel2D
{
    public:
        float xs, ys, xe, ye, xc, yc, length, angle;
        string loc;
    
    Panel2D() {}
    Panel2D(float xs1, float ys1, float xe1, float ye1)
    {
        xs = xs1, ys = ys1, xe = xe1, ye = ye1;
        xc = (xe + xs)/2.0, yc = (ye + ys)/2.0;
        length = sqrt(pow(xe - xs, 2) + pow(ye - ys, 2));
        angle = atan((ye - ys)/(xe -xs));
        loc = angle <= PI ? "upper" : "lower";
    }
};

class VortexSourcePanel2D: public Panel2D
{
    public:
        float strength = 0, vt = 0, cp = 0;

    VortexSourcePanel2D(float xs1, float ys1, float xe1, float ye1):Panel2D(xs1, ys1, xe1, ye1){}

};

class DoubletSourcePanel2D: public Panel2D
{
    public:
        float doublet_strength = 0, source_strength = 0, vt = 0, cp = 0;
        float xsl, ysl, xel, yel, xcl, ycl;

    DoubletSourcePanel2D() {}
    DoubletSourcePanel2D(float xs1, float ys1, float xe1, float ye1):Panel2D(xs1, ys1, xe1, ye1)
    {
        float xsl = 0, ysl = 0, xel = length, yel = 0;
        xcl = length/2, ycl = 0;
    }
   
    float sourceInfluence(float x, float y)
    {
        float phi_s = -1/(4*PI)*((x)*log(pow((x), 2) + pow(y, 2)) - (x - xel)*log(pow((x - xel), 2) + pow(y, 2)) + 2*y*(atan2(y, x - xel) - atan2(y, x)));

        return phi_s;
    }
    float doubletInfluence(float x, float y)
    {
        float phi_d = -1/(2*PI)*(atan2(y, x - xel) - atan2(y, x));

        return phi_d;
    }

};

class DoubletSourceSolver2D
{
    public:
        DoubletSourcePanel2D *panels;
        int num_panels;
        DoubletSourcePanel2D woke_panel;
        Uniform2D uniform;
        // DoubletSourceSolver2D() : uniform(1.0, 0.0)
    DoubletSourceSolver2D() {}

    DoubletSourceSolver2D(DoubletSourcePanel2D *slenap, int number_of_panels, Uniform2D mrofinu)
    {
        panels = slenap;
        num_panels = number_of_panels;
        woke_panel = DoubletSourcePanel2D(panels[num_panels].xe, panels[num_panels].ye, 10*panels[num_panels].xe, panels[num_panels].ye);
        uniform = mrofinu;
    }

    float *solveStrengths()
    {
        // Calculates the influence of panel j on panel i.
        float doublet_matrix[num_panels][num_panels];
        float source_matrix[num_panels][num_panels];
        float source_vector[num_panels];
        float woke_vector[num_panels];
        float b[num_panels];
        float x, y;
        for (int i = 0; i < num_panels; i++)
        {
            cout << "[ ";
            for (int j = 0; j < num_panels; j++)
            {
                // Creating doublet and source matrices
                tie(x, y) = panelCoords(panels[i].xc, panels[i].yc, panels[j].xc, panels[j].yc, panels[j].angle);
                if (i == j)
                {
                    doublet_matrix[i][j] = 0.5;
                    source_matrix[i][j] = 0.5*panels[j].doubletInfluence(x, y);
                }
                else
                    doublet_matrix[i][j] = panels[j].doubletInfluence(x, y); 
                source_matrix[i][j] = panels[j].doubletInfluence(x, y);
                cout << doublet_matrix[i][j] << " ";
            }
            cout << "] \n";
            // Source and wake vectors
            tie(x, y) = panelCoords(panels[i].xc, panels[i].yc, woke_panel.xc, woke_panel.yc, woke_panel.angle);
            woke_vector[i] = panels[i].sourceInfluence(x, y);
            source_vector[i] = uniform.magnitude*sin(uniform.angle - panels[i].angle);
        }
        cout << "[ ";
        for (int i = 0; i < num_panels; i++)
            cout << woke_vector[i] << ", ";
        cout << " ]";
        // Building LHS
        for (int i = 0; i < num_panels; i++)
        {
            doublet_matrix[i][0] = doublet_matrix[i][0] - source_vector[i];
            doublet_matrix[i][num_panels - 1] = doublet_matrix[i][num_panels - 1] + source_vector[i];
        }
        // Building RHS
        for (int i = 0; i < num_panels; i++)
            for(int j = 0; j < num_panels; j++)
                b[i] = source_matrix[i][j]*source_vector[j];
    }
};

void cosinePanels(DoubletSourcePanel2D *panels, float *x, float *y, int num_lines, int n = 40)
{
    float r = (*max_element(x,x+num_lines) - *min_element(x,x+num_lines))/2.;
    cout << *min_element(x,x+num_lines) << '\n';
    float x_center = (*max_element(x,x+num_lines) + *min_element(x,x+num_lines))/2.;
    // cout << x_center << ", " << r << "Radius" << '\n';
    float x_circ[n];
    float num_panels = n;

    for (int i = 0; i <= n; i++)
    {
        float x = i;
        x_circ[i] = x_center + r*cos(2*PI*x/(num_panels + 1));
    }
       
    float x_ends[n+1], y_ends[n+1], x_new[num_lines+1], y_new[num_lines+1];

    for(int i = 0; i <= n; i++)
    {
        x_ends[i] = x_circ[i];
        y_ends[i] = 0;
    }
    
    for(int i = 0; i <= num_lines; i++)
    {
        if (i == num_lines)
        {
            x_new[i] = x[0];
            y_new[i] = y[0];
        }
        else
        {
            x_new[i] = x[i];
            y_new[i] = y[i];
        }
    }

    int j = 1;
    float m, c;
    for (int i = 0; i < n; i++)
    {
        while (j < num_lines)
        {
            if ((x_new[j] <= x_ends[i] <= x_new[j+1]) || (x_new[j+1] <= x_ends[i] <= x_new[j]))
                break;
            else
                j++;
        }
        m = (y_new[j+1] - y_new[j])/(x_new[j+1] - x_new[j]);
        c = y_new[j+1] - m*x_new[j+1];
        y_ends[i] = m*x_ends[i] + c;
        y_ends[n] = y_ends[0];
    }

    for (int i = 0; i < n; i++)
    {
        panels[i] = DoubletSourcePanel2D(x_ends[i], y_ends[i], x_ends[i+1], y_ends[i+1]);
        cout << "(" << panels[i].xs << ", " << panels[i].ys << ", " << panels[i].xe << ", " << panels[i].ye << ")" << '\n';
    }
}

int readCoords(string filename, float **x, float **y)
{
    ifstream file;
    string line;
    file.open(filename);
    int number_of_lines = 0;
    if (file.is_open())
    {   
        while (getline(file, line))
            ++number_of_lines;
        cout << number_of_lines << '\n';
        file.close();
    }

    *x = new float[number_of_lines];
    *y = new float[number_of_lines];
    file.open(filename);
    if (file.is_open())
    {
        int i = 0;
        while (getline(file,line))
        {
            float val1, val2;
            file >> val1 >> val2;
            (*x)[i] = val1;
            (*y)[i] = val2;
            i++;
        }
        file.close();
    }
    return number_of_lines;
}

int main()
{       
    float *x, *y;
    int number_of_lines = readCoords("resources/ClarkY.dat", &x, &y);
    // Panel2D first(0.0, 0.0, 1.0, 1.0);
    // cout << first.length << endl << first.angle << endl << first.loc << endl;
    // VortexSourcePanel2D second(1.0, 1.0, 1.0, 2.0);
    // cout << second.length << endl << second.angle << endl << second.loc << endl << second.cp << endl;   
    Uniform2D uniform(1.0, 5.0);
    // tuple<float, float> tie;
    // tie = uniform.velocity(2.0, 3.0);
    // cout << get<0>(tie) << " " << get<1>(tie) << endl;
    // Doublet2D doublet(1.0, 1.0, 1.0);
    // tuple<float, float> doubled;
    // doubled = doublet.velocity(2.0, 3.0);
    // cout << get<0>(doubled) << " " << get<1>(doubled) << endl;

    int n_panels = 10;
    DoubletSourcePanel2D panels[n_panels];
    cosinePanels(panels, x, y, number_of_lines, n_panels);
    // float n = n_panels;
    // for (int i = 0; i < n_panels; i++)
    // {
    //     float x = i;
    //     // panels[i] = DoubletSourcePanel2D(x/n, 0, 1 - (x + 1)/n, 0);
    //     cout << "Panel "<< i + 1 << ": (" << panels[i].xs << ", " << panels[i].xe << ")" << "\n";
    // }
    DoubletSourceSolver2D airfoil;
    airfoil = DoubletSourceSolver2D(panels, n_panels, uniform);
    airfoil.solveStrengths();
}