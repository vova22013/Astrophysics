#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cassert>

// Task 1
// Periodic orbits
 
// Global constants

const int m1 = 1;
const int m2 = 499;
const int rOrbit = 8200;
const double omega = 0.02829174;
const double G = 0.004535;
const double alpha = -0.001976995;
const double mu = 1 + double(m1) / m2;
const double absKsi_t = 10.44666776;
const double eps_t = -(3. / 2) * G * m2 / (pow(mu, 2) * absKsi_t);
const double eps0 = -0.67 * eps_t; 
const double precision = 1e-10;
const double dt = 1e-5;

// Override operators

std::vector<double> operator+(const std::vector<double>& vec1,
    const std::vector<double>& vec2) {
    int size = vec1.size();
    if (vec1.size() > vec2.size()) size = vec2.size();
    std::vector<double> res;
    for (int i = 0; i < size; ++i) {
        res.push_back(vec1[i] + vec2[i]);
    }
    return res;
}

std::vector<double> operator*(const std::vector<double>& vec1,
    const std::vector<double>& vec2) {
    int size = vec1.size();
    if (vec1.size() > vec2.size()) size = vec2.size();
    std::vector<double> res;
    for (int i = 0; i < size; ++i) {
        res.push_back(vec1[i] * vec2[i]);
    }
    return res;
}

std::vector<double> operator*(const std::vector<double>& vec,
    double scalar) {
    std::vector<double> res = vec;
    for (auto& elem : res) {
        elem *= scalar;
    }
    return res;
}

// Create start point

std::vector<double> startPoint(double C1, double start = 0) {
    double ksi = -absKsi_t + C1;
    double eta = start;
    double v_ksi = 0;
    double x1 = eps0 + (G * m2) / (pow(mu, 2) * abs(ksi)),
        x2 = -alpha * pow(ksi, 2);
    double v_eta = sqrt(2 * x1 + x2);
    std::vector<double> coords = { ksi, eta, v_ksi, v_eta };
    return coords;
}

// Distance to the starting point

double r(double ksi, double eta) {
    return sqrt(pow(eta, 2) + pow(ksi, 2));
}

// Acceleration

double aKsi(double v_eta, double ksi, double eta) {
    double mult = G * m2 / (pow(mu, 2) * pow(r(ksi, eta), 3));
    return 2 * omega * v_eta - ksi * (mult + alpha);
}

double aEta(double v_ksi, double ksi, double eta) {
    double mult = G * m2 / (pow(mu, 2) * pow(r(ksi, eta), 3));
    return -2 * omega * v_ksi - mult * eta;
}

// Check critical energy

double energy(std::vector<double> coords) {
    double T = 0.5 * (pow(coords[2], 2) + pow(coords[3], 2));
    double U = -G * m2 / (pow(mu, 2) * r(coords[0], coords[1]));
    double W = 0.5 * alpha * pow(coords[0], 2);
    return T + U + W;
}

bool checkEnergy(double eps) {
    return abs((eps - eps0) / eps0) < precision;
}

// System of equations 1

std::vector<double> f1(const std::vector<double>& X) {
    // X[0] = ksi, X[1] = eta, X[2] = v_ksi, X[3] = v_eta
    std::vector<double> vec(4);
    vec[0] = X[2];
    vec[1] = X[3];
    vec[2] = aKsi(X[3], X[0], X[1]);
    vec[3] = aEta(X[2], X[0], X[1]);
    return vec;
}

// Algorithm Runge-Kutta

std::vector<double> R_K(const std::vector<double> &X, double _dt) {
    auto k1 = f1(X) * _dt;
    auto k2 = f1(X + k1 * 0.5) * _dt;
    auto k3 = f1(X + k2 * 0.5) * _dt;
    auto k4 = f1(X + k3) * _dt;
    auto K = k1 + ((k2 + k3) * 2.) + k4;
    double mult = 1. / 6;
    auto X_dt = X + K * mult;
    return X_dt;
}

// Linear approximation 

std::pair<double, double> coeffsOfLine(double x1, double y1,
    double x2, double y2) {
    std::pair<double, double> coeffs;
    coeffs.first = (y2 - y1) / (x2 - x1);
    coeffs.second = (y2 * x1 - y1 * x2) / (x1 - x2);
    return coeffs;
}

const double yOnLine(double x, double a, double b) {
    return a * x + b;
}

const double xOnLine(double y, double a, double b) {
    return (y - b) / a;
}

// Search minimal distance 

std::vector<double> searchMinDist(std::vector<std::vector<double>> vec) {
    // vec[0] = ksi, vec[1] = eta, vec[2] = v_ksi, vec[3] = v_eta,
    // vec[4] = C1, vec[5] = count * dt, vec[6] = ksi_prev, vec[7] = eta_prev,
    // vec[8] = v_ksi_prev, vec[9] = v_eta_prev
    std::vector<double> minDistVec = vec[0]; 
    double minDist;
    for (auto elem : vec) {
        auto startCoords = startPoint(elem[4]);
        double dist1 = pow(elem[6] - startCoords[0], 2) + pow(elem[7] - startCoords[1], 2);
        double dist2 = pow(elem[0] - startCoords[0], 2) + pow(elem[1] - startCoords[1], 2);
        double dist3 = pow(elem[0] - elem[6], 2) + pow(elem[1] - elem[7], 2);

        double dist = dist1 + dist2 + dist3;
        if (elem == vec[0]) minDist = dist;
        if (minDist > dist) {
            minDist = dist;
            minDistVec = elem;
        }
    }
    return minDistVec;
};

// Search minimal velocities

std::vector<double> searchMinVel(std::vector<std::vector<double>> vec) {
    // vec[0] = ksi, vec[1] = eta, vec[2] = v_ksi, vec[3] = v_eta,
    // vec[4] = C1, vec[5] = count * dt, vec[6] = ksi_prev, vec[7] = eta_prev,
    // vec[8] = v_ksi_prev, vec[9] = v_eta_prev
    std::vector<double> minVelVec = vec[0];
    double minVel;
    for (auto elem : vec) {
        auto startCoords = startPoint(elem[4]);
        double vel1 = pow(elem[8] - startCoords[2], 2) + pow(elem[9] - startCoords[3], 2);
        double vel2 = pow(elem[2] - startCoords[2], 2) + pow(elem[3] - startCoords[3], 2);
        double vel3 = pow(elem[2] - elem[8], 2) + pow(elem[3] - elem[9], 2);

        double vel = vel1 + vel2 + vel3;
        if (elem == vec[0]) minVel = vel;
        if (minVel > vel) {
            minVel = vel;
            minVelVec = elem;
        }
    }
    return minVelVec;
};

// Search minimal R

std::vector<double> searchMinVelAndDist(std::vector<std::vector<double>> vec) {
    // vec[0] = ksi, vec[1] = eta, vec[2] = v_ksi, vec[3] = v_eta,
    // vec[4] = C1, vec[5] = count * dt, vec[6] = ksi_prev, vec[7] = eta_prev,
    // vec[8] = v_ksi_prev, vec[9] = v_eta_prev
    std::vector<double> minVec = vec[0];
    double minScore, minDist, minVel;
    for (auto elem : vec) {
        auto startCoords = startPoint(elem[4]);
        double dist1 = pow(elem[6] - startCoords[0], 2) + pow(elem[7] - startCoords[1], 2);
        double dist2 = pow(elem[0] - startCoords[0], 2) + pow(elem[1] - startCoords[1], 2);
        double dist3 = 0; // pow(elem[0] - elem[6], 2) + pow(elem[1] - elem[7], 2);
        double deviderCoords = 1; // pow(startCoords[0], 2) + pow(startCoords[1], 2);

        double vel1 = pow(elem[8] - startCoords[2], 2) + pow(elem[9] - startCoords[3], 2);
        double vel2 = pow(elem[2] - startCoords[2], 2) + pow(elem[3] - startCoords[3], 2);
        double vel3 = 0; // pow(elem[2] - elem[8], 2) + pow(elem[3] - elem[9], 2);
        double deviderVels = pow(startCoords[2], 2) + pow(startCoords[3], 2);

        double vel = (vel1 + vel2 + vel3) / deviderVels;
        double dist = (dist1 + dist2 + dist3) / deviderCoords;
        double score = vel + dist;
        if (elem == vec[0]) {
            minScore = score;
        }
        if (minScore > score) {
            minScore = score;
            minVec = elem;
        }
    }
    return minVec;
};

std::vector<double> searchMinDistWithLinearApprx(std::vector<std::vector<double>> vec) {
    std::vector<double> minVec = vec[0];
    double minScore;
    for (auto elem : vec) {
        auto startCoords = startPoint(elem[4]);

        auto coeffsCoords = coeffsOfLine(elem[6], elem[7], elem[0], elem[1]);
        auto ksi0 = startCoords[0];
        auto eta0 = startCoords[1];
        auto ksiExpect = xOnLine(eta0, coeffsCoords.first, coeffsCoords.second);
        
        auto coeffsVels = coeffsOfLine(elem[8], elem[9], elem[2], elem[3]);
        auto v_ksi0 = startCoords[2];
        auto v_eta0 = startCoords[3];
        auto v_etaExpect = yOnLine(v_ksi0, coeffsVels.first, coeffsVels.second);
        
        double distScore = pow(ksiExpect - ksi0, 2) / pow(ksi0, 2);
        double velScore = pow(v_etaExpect - v_eta0, 2) / pow(v_eta0, 2);
        auto score = distScore + velScore;

        if (elem == vec[0]) minScore = score;

        if (score < minScore) {
            minScore = score;
            minVec = elem;
        }
    }
    return minVec;
};

std::vector<double> searchMinRRatio(std::vector<std::vector<double>> vec) {
    // vec[0] = ksi, vec[1] = eta, vec[2] = v_ksi, vec[3] = v_eta,
    // vec[4] = C1, vec[5] = count * dt, vec[6] = ksi_prev, vec[7] = eta_prev,
    // vec[8] = v_ksi_prev, vec[9] = v_eta_prev
    double minScore;
    std::vector<double> minVec = vec[0];
    for (auto elem: vec) {
        auto startCoords = startPoint(elem[4]);
        double coordsRatio = pow(elem[0], 2) + pow(elem[1], 2);
        double velsRatio = pow(elem[2], 2) + pow(elem[3], 2);
        double coordsRatio0 = pow(startCoords[0], 2) + pow(startCoords[1], 2);
        double velsRatio0 = pow(startCoords[2], 2) + pow(startCoords[3], 2);

        double score = coordsRatio / coordsRatio0 + 
            velsRatio / velsRatio0;
        if (elem == vec[0]) minScore = score;
        if (minScore > score) {
            minScore = score;
            minVec = elem;
        }
    }
    return minVec;
};

// Other functions

int choosePowerOfScale(std::string C1) {
    int indPoint = C1.find('.');
    if (indPoint == -1) return 0;
    return C1.substr(indPoint + 1).length();
}

std::vector<double> calculationOrbit(std::vector<double> coords, double C1,
    bool second_cond = false) {
    double numOfHundreds = 1.9, scale = 60.;
    uint32_t count = 0;
    uint32_t minTime = uint32_t(scale * pow(dt, -1));
    uint32_t maxTime = uint32_t((100. / scale) * numOfHundreds * double(minTime));
    
    std::vector<double> orbit;
    while (true) {
        double ksi_prev = coords[0];
        double eta_prev = coords[1];
        double v_ksi_prev = coords[2];
        double v_eta_prev = coords[3];
        
        ++count;
        if (count == maxTime) break;

        coords = R_K(coords, dt);
        double eps = energy(coords);
        if (!checkEnergy(eps)) break;

        bool endIterCoords = eta_prev < 0 && coords[1] > 0 && coords[0] < 0;
        bool endIterVels = v_ksi_prev < 0 && coords[2] > 0 && coords[3] > 0;
        bool endIter = endIterVels;
        if (second_cond) endIter = endIterVels && endIterCoords;
        if (endIter) {
            double period = count * dt;
            orbit = { coords[0], coords[1], coords[2], coords[3],
                C1, period, ksi_prev, eta_prev, v_ksi_prev, v_eta_prev };
            break;
        }
    }
    return orbit;
}

const void writeResultInFile(std::string C1_str, int& powOfScale,
    std::vector<std::vector<double>>& vecMinR, std::vector<double>& res) {
    std::string address = "C:\\Programs\\code\\MinDist_" + C1_str + "_.txt";
    std::ofstream fout;
    if (powOfScale == 1) fout.open(address, std::ios::trunc);
    else fout.open(address, std::ios::app);
    fout << "Step #" << powOfScale - 1 << ".\n";
    fout << std::showpos << std::showpoint << std::setprecision(8)
        << "Result: C1: " << res[4] << " Ksi: " << res[0]
        << " Eta: " << res[1] << std::endl;
    for (auto minR : vecMinR)
    {
        auto startCoords = startPoint(minR[4]);
        fout << std::showpos << std::showpoint
            << std::setprecision(8) << std::scientific
            << "\nC1:" << minR[4] << " P:" << minR[5] << "m.y."
            << " Ksi:" << minR[0] << " Eta:" << minR[1]
            << " dKsiPrev:" << abs(minR[6] - startCoords[0]) 
            << " dEtaPrev:" << abs(minR[7] - startCoords[1]) 
            << " dKsiNext:" << abs(minR[0] - startCoords[0])
            << " dEtaNext:" << abs(minR[1] - startCoords[1]);
    }
    fout << "\n\n";
    fout.close();
}


// Task 2
// Research on the stability of the orbit

// Auxiliary formulas

double mult1(double x1, double x2, double betta) {
    double num = 3 * pow(x1, 2);
    double denum = pow(r(x1, x2), 2);
    return betta * (1 - num / denum);
}

double mult2(double x1, double x2, double betta) {
    double num = 3 * betta * x1 * x2;
    double denum = pow(r(x1, x2), 2);
    return num / denum;
}

double uVariation(double dksi,  double deta, double dv,
    double ksi, double eta, double betta) {
    double sum_u_1 = 2 * omega * dv;
    double sum_u_2 = -(mult1(ksi, eta, betta) + alpha) * dksi;
    double sum_u_3 = mult2(ksi, eta, betta) * deta;
    return sum_u_1 + sum_u_2 + sum_u_3;
}

double vVariation(double dksi, double deta, double du,
    double ksi, double eta, double betta) {
    double sum_v_1 = -2 * omega * du;
    double sum_v_2 = mult2(ksi, eta, betta) * dksi;
    double sum_v_3 = -mult1(eta, ksi, betta) * deta;
    return sum_v_1 + sum_v_2 + sum_v_3;
}

// System of equations 2

std::vector<double> f2(const std::vector<double>& vars,
    const std::vector<double>& coords) {
    std::vector<double> X(4);
    X[0] = vars[2];
    X[1] = vars[3];
    double betta = G * m2 / (pow(mu, 2) * pow(r(coords[0], coords[1]), 3));
    X[2] = uVariation(vars[0], vars[1], vars[3], coords[0], coords[1], betta);
    X[3] = vVariation(vars[0], vars[1], vars[2], coords[0], coords[1], betta);
    return X;
}

// Algorithm Runge-Kutta for the variations

std::vector<double> R_K(const std::vector<double>& vars, 
    const std::vector<double>& coords, double _dt) {
    auto k1 = f2(vars, coords) * _dt;
    auto k2 = f2(vars + k1 * 0.5, coords) * _dt;
    auto k3 = f2(vars + k2 * 0.5, coords) * _dt;
    auto k4 = f2(vars + k3, coords) * _dt;
    auto K = k1 + ((k2 + k3) * 2.) + k4;
    double mult = 1. / 6;
    auto vars_dt = vars + K * mult;
    return vars_dt;
}

int main()
{
    // std::string C1_str = "-10";
    std::string C1_str = "-8.4730036";
    bool searchIsComplete = true;
    
    int const maxScale = 8;
    double C1 = std::stod(C1_str);
    std::vector<double> res = {0., 0., 0., 0., C1};
    if (!searchIsComplete) {
        int powOfScale = choosePowerOfScale(C1_str);
        while (powOfScale < maxScale) {
            std::vector<std::vector<double>> vecMinR, startCoords;
            
            double scale = pow(0.1, powOfScale - 1);
            double step = pow(0.1, powOfScale);
            ++powOfScale;
            double start = res[4] - scale;
            double end = res[4] + scale;
            C1 = start;

            bool matchCoords = false;
            if (powOfScale > 6) matchCoords = true;
            while (C1 < end) {
                C1 += step;
                std::cout << std::showpoint << std::setprecision(8)
                    << "C1: " << C1 << std::endl;
                auto coords = startPoint(C1);
                auto minR = calculationOrbit(coords, C1, matchCoords);
                if (minR.size() == 0) continue;
                vecMinR.push_back(minR);
            }
            if (vecMinR.size() == 0) return 1;
            res = searchMinDistWithLinearApprx(vecMinR);
            writeResultInFile(C1_str, powOfScale, vecMinR, res);
            std::cout << std::endl;
        }
    }
    C1 = res[4];
    int period;
    if (searchIsComplete) {
        std::ifstream fin("C:\\Programs\\code\\Period.txt", std::ios::in);
        assert(fin.is_open()); fin >> period; fin.close();
    }
    else {
        period = int(res[5] / dt);
        std::ofstream fout("C:\\Programs\\code\\Period.txt", 'w');
        fout << period; fout.close();
    }
    auto coords = startPoint(C1);
    uint32_t count = 0;
    std::cout << std::showpoint << std::setprecision(10)
        << "\nC1 = " << C1 << std::endl;
    
    auto address1 = "C:\\Programs\\code\\Coordinates_" + C1_str + "_.txt";
    std::ofstream fout1(address1, 'w');
    auto address2 = "C:\\Programs\\code\\DifEnergy_" + C1_str + "_.txt";
    std::ofstream fout2(address2, 'w');
    fout1 << "Ksi: \t       " << "Eta:\n";
    fout2 << "Period:       " << "Eps - Eps0:\n";
    while (count != period) {
        ++count;
        coords = R_K(coords, dt);
        double eps = energy(coords);
        double difEps = abs(eps - eps0);
        
        fout1 << std::showpoint << std::setprecision(7) << std::scientific
            << std::showpos << coords[0] << ' ' << coords[1] << '\n';
        fout2 << std::showpoint << std::setprecision(7) << std::scientific
            << count * dt << ' ' << difEps << '\n';
    }
    fout1.close();
    fout2.close();


    // Task 2
    // Research on the stability of the orbit
    /*auto coords1 = startPoint(C1);
    
    std::vector<std::vector<double>> mtrxVars = {
        {1., 0., 0., 0.},
        {0., 1., 0., 0.},
        {0., 0., 1., 0.},
        {0., 0., 0., 1.}
    };
    count = 0;
    while (count != period) {
        coords1 = R_K(coords1, dt);
        mtrxVars[0] = R_K(mtrxVars[0], coords1, dt);
        mtrxVars[1] = R_K(mtrxVars[1], coords1, dt);
        mtrxVars[2] = R_K(mtrxVars[2], coords1, dt);
        mtrxVars[3] = R_K(mtrxVars[3], coords1, dt);
        ++count;
    }*/
    

};
        

