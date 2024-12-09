#include "function_1task.h"

// Task 1
// Periodic orbits
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
};
