#ifndef POOL_BOILING_LBM_HPP
#define POOL_BOILING_LBM_HPP

#include "main.hpp"

class PoolBoilingLBM
{
public:
    static constexpr int NX = 160;
    static constexpr int NY = 100;
    static constexpr int Q = 9;
    static constexpr int N = NX * NY;

    PoolBoilingLBM();

    void initialize();
    void step();
    void writeFieldCSV(int step_id) const;
    void writeDiagnosticsHeader() const;
    void appendDiagnostics(int step_id) const;

private:
    inline static constexpr std::array<int, Q> cx{0, 1, 0, -1, 0, 1, -1, -1, 1};
    inline static constexpr std::array<int, Q> cy{0, 0, 1, 0, -1, 1, 1, -1, -1};
    inline static constexpr std::array<int, Q> opp{0, 3, 4, 1, 2, 7, 8, 5, 6};
    inline static constexpr std::array<double, Q> w{
        4.0 / 9.0,
        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};

    std::vector<double> f;
    std::vector<double> f_new;
    std::vector<double> ux;
    std::vector<double> uy;
    std::vector<double> phi;
    std::vector<double> temp;
    std::vector<int> solid;
    int step_count;

    static int id(int x, int y);
    static int px(int x);
    static int clampY(int y);

    double feq(int k, double rho, double ux, double uy) const;
    double interfaceIndicator(double phi) const;

    void computeVelocity();
    void evolveTemperature();
    void evolvePhaseField();
    void collideAndStream();
    void applyBoundaryConditions();
};

#endif
