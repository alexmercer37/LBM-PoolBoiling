#include "../inc/PoolBoilingLBM.hpp"

namespace
{

    template <typename T>
    T clampValue(T value, T lower, T upper)
    {
        if (value < lower)
            return lower;
        if (value > upper)
            return upper;
        return value;
    }

    constexpr double RHO0 = 1.0;
    constexpr double TAU = 0.78;
    constexpr double OMEGA = 1.0 / TAU;
    constexpr double VMAX = 0.10;

    constexpr double PHI_LIQUID = 1.0;

    constexpr double MOBILITY = 0.010;
    constexpr double DETACHED_MOBILITY = 0.0065;
    constexpr double EPS2 = 2.18;

    constexpr double SURFACE_TENSION = 8.6e-4;

    constexpr double T_SAT = 1.0;
    constexpr double T_WALL = 1.010;
    constexpr double T_TOP = 0.985;

    constexpr double THERMAL_DIFF = 0.045;
    constexpr double THERMAL_ADV = 0.018;
    constexpr double LATENT_COOLING = 0.016;

    constexpr double TEMP_MIN = 0.965;
    constexpr double TEMP_MAX = 1.155;

    constexpr int NUCLEATION_X = 80;

    constexpr double INITIAL_BUBBLE_RADIUS = 8.2;
    constexpr double INITIAL_BUBBLE_Y = 5.8;
    constexpr double INITIAL_INTERFACE_WIDTH = 1.7;

    constexpr int HOT_HALF_WIDTH = 8;
    constexpr int HEATER_LATERAL_HALF_WIDTH = 18;
    constexpr int EVAP_HEIGHT = 12;

    constexpr double HEATER_BOOST_Y1 = 0.090;
    constexpr double HEATER_BOOST_Y2 = 0.058;
    constexpr double HEATER_BOOST_Y3 = 0.034;

    constexpr int GROWTH_END = 2600;
    constexpr int PINCH_END = 5200;

    constexpr int RISE_DAMPING_START = 3600;
    constexpr double RISE_VELOCITY_DAMPING = 0.993;

    constexpr double BUOYANCY_GROWTH = 2.4e-5;
    constexpr double BUOYANCY_PINCH = 5.2e-5;
    constexpr double BUOYANCY_RISE = 2.5e-5;

    constexpr double EVAP_GROWTH = 0.082;
    constexpr double EVAP_PINCH = 0.018;
    constexpr double EVAP_RISE = 0.020;

    constexpr double WALL_EVAP_BOOST_GROWTH = 2.35;
    constexpr double WALL_EVAP_BOOST_PINCH = 1.15;
    constexpr double WALL_EVAP_BOOST_RISE = 1.00;

    constexpr double COND_GROWTH = 0.0028;
    constexpr double COND_PINCH = 0.0038;
    constexpr double COND_RISE = 0.00030;

    constexpr int COND_UPPER_Y = 78;
    constexpr double COND_UPPER_FACTOR = 1.00;
    constexpr double COND_LOWER_FACTOR = 0.22;
    constexpr double COND_VAPOR_REDUCTION = 0.25;

    constexpr int DETACHED_SUPPORT_START = 3800;
    constexpr double DETACHED_SUPPORT = 0.0022;
    constexpr int DETACHED_SUPPORT_Y_MIN = 12;
    constexpr int DETACHED_SUPPORT_Y_MAX = 82;
    constexpr double DETACHED_SUPPORT_PHI_LIMIT = 0.60;

    constexpr int NECK_ZONE_Y_MIN = 2;
    constexpr int NECK_ZONE_Y_MAX = 8;
    constexpr int NECK_HALF_WIDTH = 11;

    constexpr double NECK_REWETTING = 0.00110;
    constexpr double STEM_SUPPRESSION = 0.00040;

    constexpr double RHO_LIQUID_VIEW = 1.0;
    constexpr double RHO_VAPOR_VIEW = 0.12;

    enum class Stage
    {
        Growth,
        Pinch,
        Rise
    };

    struct PhaseChangeParams
    {
        double evaporation;
        double condensation;
        double wall_boost;
    };

    Stage currentStage(int step)
    {
        if (step < GROWTH_END)
            return Stage::Growth;

        if (step < PINCH_END)
            return Stage::Pinch;

        return Stage::Rise;
    }

    PhaseChangeParams phaseChangeParams(Stage stage)
    {
        switch (stage)
        {
        case Stage::Growth:
            return {EVAP_GROWTH, COND_GROWTH, WALL_EVAP_BOOST_GROWTH};

        case Stage::Pinch:
            return {EVAP_PINCH, COND_PINCH, WALL_EVAP_BOOST_PINCH};

        case Stage::Rise:
        default:
            return {EVAP_RISE, COND_RISE, WALL_EVAP_BOOST_RISE};
        }
    }

    double buoyancyForStage(Stage stage)
    {
        switch (stage)
        {
        case Stage::Growth:
            return BUOYANCY_GROWTH;

        case Stage::Pinch:
            return BUOYANCY_PINCH;

        case Stage::Rise:
        default:
            return BUOYANCY_RISE;
        }
    }
}

PoolBoilingLBM::PoolBoilingLBM()
    : f(N * Q, 0.0),
      f_new(N * Q, 0.0),
      ux(N, 0.0),
      uy(N, 0.0),
      phi(N, PHI_LIQUID),
      temp(N, T_TOP),
      solid(N, 0),
      step_count(0)
{
}

int PoolBoilingLBM::id(int x, int y)
{
    return y * NX + x;
}

int PoolBoilingLBM::px(int x)
{
    if (x < 0)
        return NX - 1;

    if (x >= NX)
        return 0;

    return x;
}

int PoolBoilingLBM::clampY(int y)
{
    return clampValue(y, 1, NY - 2);
}

double PoolBoilingLBM::feq(int k, double rho, double ux_value, double uy_value) const
{
    const double cu = 3.0 * (cx[k] * ux_value + cy[k] * uy_value);
    const double u2 = ux_value * ux_value + uy_value * uy_value;

    return w[k] * rho * (1.0 + cu + 0.5 * cu * cu - 1.5 * u2);
}

double PoolBoilingLBM::interfaceIndicator(double phi_value) const
{
    const double v = clampValue(phi_value, -1.0, 1.0);
    return std::max(0.0, 1.0 - v * v);
}

void PoolBoilingLBM::initialize()
{
    _mkdir("csv_out");

    for (int y = 0; y < NY; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            const int p = id(x, y);

            solid[p] = (y == 0 || y == NY - 1) ? 1 : 0;
            ux[p] = 0.0;
            uy[p] = 0.0;
            phi[p] = PHI_LIQUID;

            temp[p] = T_TOP +
                      (T_WALL - T_TOP) *
                          (1.0 - static_cast<double>(y) / static_cast<double>(NY - 1));
        }
    }

    const double x0 = static_cast<double>(NUCLEATION_X);

    for (int y = 1; y < NY - 1; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            const double dx = static_cast<double>(x) - x0;
            const double dy = static_cast<double>(y) - INITIAL_BUBBLE_Y;
            const double r = std::sqrt(dx * dx + dy * dy);

            const double profile =
                std::tanh((r - INITIAL_BUBBLE_RADIUS) / INITIAL_INTERFACE_WIDTH);

            phi[id(x, y)] = clampValue(profile, -1.0, 1.0);
        }
    }

    for (int dx = -HOT_HALF_WIDTH; dx <= HOT_HALF_WIDTH; ++dx)
    {
        const int x = px(NUCLEATION_X + dx);

        temp[id(x, 1)] = T_WALL + HEATER_BOOST_Y1;
        temp[id(x, 2)] = T_WALL + HEATER_BOOST_Y2;
        temp[id(x, 3)] = T_WALL + HEATER_BOOST_Y3;
    }

    for (int p = 0; p < N; ++p)
    {
        for (int k = 0; k < Q; ++k)
        {
            f[p * Q + k] = feq(k, RHO0, 0.0, 0.0);
        }
    }
}

void PoolBoilingLBM::computeVelocity()
{
    for (int y = 1; y < NY - 1; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            const int p = id(x, y);

            if (solid[p])
                continue;

            double rho = 0.0;
            double jx = 0.0;
            double jy = 0.0;

            for (int k = 0; k < Q; ++k)
            {
                const double fk = f[p * Q + k];

                rho += fk;
                jx += fk * cx[k];
                jy += fk * cy[k];
            }

            rho = std::max(1.0e-12, rho);

            ux[p] = clampValue(jx / rho, -VMAX, VMAX);
            uy[p] = clampValue(jy / rho, -VMAX, VMAX);

            if (step_count >= RISE_DAMPING_START)
            {
                const double vapor_fraction =
                    clampValue(0.5 * (1.0 - phi[p]), 0.0, 1.0);

                const double damping =
                    1.0 - vapor_fraction * (1.0 - RISE_VELOCITY_DAMPING);

                ux[p] *= damping;
                uy[p] *= damping;
            }
        }
    }
}

void PoolBoilingLBM::evolveTemperature()
{
    std::vector<double> next = temp;

    for (int y = 1; y < NY - 1; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            const int p = id(x, y);

            if (solid[p])
                continue;

            const double lap =
                temp[id(px(x + 1), y)] +
                temp[id(px(x - 1), y)] +
                temp[id(x, y + 1)] +
                temp[id(x, y - 1)] -
                4.0 * temp[p];

            const int xu = px(x - (ux[p] > 0.0 ? 1 : -1));
            const int yu = clampY(y - (uy[p] > 0.0 ? 1 : -1));

            const double adv =
                -THERMAL_ADV * (temp[p] - temp[id(xu, yu)]);

            const double cooling =
                LATENT_COOLING *
                std::max(0.0, T_SAT - temp[p]) *
                interfaceIndicator(phi[p]);

            next[p] = temp[p] + THERMAL_DIFF * lap + adv - cooling;
        }
    }

    for (int x = 0; x < NX; ++x)
    {
        next[id(x, 1)] = T_WALL;
        next[id(x, NY - 2)] = T_TOP;
    }

    for (int dx = -HOT_HALF_WIDTH; dx <= HOT_HALF_WIDTH; ++dx)
    {
        const int x = px(NUCLEATION_X + dx);

        next[id(x, 1)] = T_WALL + HEATER_BOOST_Y1;
        next[id(x, 2)] = T_WALL + HEATER_BOOST_Y2;
        next[id(x, 3)] = T_WALL + HEATER_BOOST_Y3;
    }

    for (double &t : next)
    {
        t = clampValue(t, TEMP_MIN, TEMP_MAX);
    }

    temp.swap(next);
}

void PoolBoilingLBM::evolvePhaseField()
{
    std::vector<double> next = phi;

    const Stage stage = currentStage(step_count);
    const PhaseChangeParams params = phaseChangeParams(stage);

    for (int y = 1; y < NY - 1; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            const int p = id(x, y);

            if (solid[p])
                continue;

            const double ph = phi[p];

            const double lap =
                phi[id(px(x + 1), y)] +
                phi[id(px(x - 1), y)] +
                phi[id(x, y + 1)] +
                phi[id(x, y - 1)] -
                4.0 * ph;

            const double gradx =
                0.5 * (phi[id(px(x + 1), y)] -
                       phi[id(px(x - 1), y)]);

            const double grady =
                0.5 * (phi[id(x, y + 1)] -
                       phi[id(x, y - 1)]);

            const double adv = -(ux[p] * gradx + uy[p] * grady);
            const double double_well = ph * ph * ph - ph;

            const double mobility =
                step_count >= DETACHED_SUPPORT_START ? DETACHED_MOBILITY : MOBILITY;

            const double ac = mobility * (EPS2 * lap - double_well);

            const double I = interfaceIndicator(ph);
            const double superheat = std::max(0.0, temp[p] - (T_SAT + 0.010));
            const double subcool = std::max(0.0, T_SAT - temp[p]);

            const bool near_wall = y <= 6;
            const double dx_center =
                std::abs(static_cast<double>(x - NUCLEATION_X));

            double evaporation = 0.0;

            if (y <= EVAP_HEIGHT && dx_center <= HEATER_LATERAL_HALF_WIDTH)
            {
                evaporation =
                    -params.evaporation *
                    superheat *
                    I *
                    (near_wall ? params.wall_boost : 1.0);
            }

            if (step_count >= DETACHED_SUPPORT_START &&
                y >= DETACHED_SUPPORT_Y_MIN &&
                y <= DETACHED_SUPPORT_Y_MAX &&
                ph < DETACHED_SUPPORT_PHI_LIMIT)
            {
                const double vapor_fraction =
                    clampValue(0.5 * (1.0 - ph), 0.0, 1.0);

                evaporation +=
                    -DETACHED_SUPPORT *
                    vapor_fraction *
                    (0.65 * I + 0.08 * vapor_fraction);
            }

            double condensation_factor =
                y > COND_UPPER_Y ? COND_UPPER_FACTOR : COND_LOWER_FACTOR;

            if (ph < 0.0)
            {
                condensation_factor *= COND_VAPOR_REDUCTION;
            }

            const double condensation =
                params.condensation * subcool * I * condensation_factor;

            double neck_rewetting = 0.0;
            double stem_suppression = 0.0;

            if (stage == Stage::Pinch &&
                y >= NECK_ZONE_Y_MIN &&
                y <= NECK_ZONE_Y_MAX &&
                dx_center <= static_cast<double>(NECK_HALF_WIDTH))
            {
                neck_rewetting = NECK_REWETTING * I;

                if (ph < 0.20)
                {
                    stem_suppression = STEM_SUPPRESSION * I;
                }
            }

            next[p] = clampValue(
                ph + adv + ac + evaporation + condensation + neck_rewetting + stem_suppression,
                -1.0,
                1.0);
        }
    }

    for (int x = 0; x < NX; ++x)
    {
        next[id(x, 1)] = std::max(next[id(x, 1)], -0.82);
        next[id(x, NY - 2)] = PHI_LIQUID;
    }

    if (step_count < GROWTH_END)
    {
        for (int dx = -4; dx <= 4; ++dx)
        {
            const int p = id(px(NUCLEATION_X + dx), 2);
            next[p] = std::min(next[p], -0.34);
        }
    }

    phi.swap(next);
}

void PoolBoilingLBM::collideAndStream()
{
    std::fill(f_new.begin(), f_new.end(), 0.0);

    const Stage stage = currentStage(step_count);
    const double buoyancy = buoyancyForStage(stage);

    for (int y = 1; y < NY - 1; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            const int p = id(x, y);

            if (solid[p])
                continue;

            const double gradx =
                0.5 * (phi[id(px(x + 1), y)] -
                       phi[id(px(x - 1), y)]);

            const double grady =
                0.5 * (phi[id(x, y + 1)] -
                       phi[id(x, y - 1)]);

            const double lap =
                phi[id(px(x + 1), y)] +
                phi[id(px(x - 1), y)] +
                phi[id(x, y + 1)] -
                4.0 * phi[p] +
                phi[id(x, y - 1)];

            const double curvature_proxy = -lap;

            double Fx = SURFACE_TENSION * curvature_proxy * gradx;
            double Fy = SURFACE_TENSION * curvature_proxy * grady;

            const double vapor_fraction = 0.5 * (1.0 - phi[p]);
            Fy += buoyancy * vapor_fraction;

            const double uxf =
                clampValue(ux[p] + 0.5 * Fx / RHO0, -VMAX, VMAX);

            const double uyf =
                clampValue(uy[p] + 0.5 * Fy / RHO0, -VMAX, VMAX);

            for (int k = 0; k < Q; ++k)
            {
                const double eiF = cx[k] * Fx + cy[k] * Fy;
                const double force_term = 3.0 * w[k] * eiF;

                const double post =
                    f[p * Q + k] -
                    OMEGA * (f[p * Q + k] - feq(k, RHO0, uxf, uyf)) +
                    force_term;

                const int xn = px(x + cx[k]);
                const int yn = y + cy[k];

                if (yn <= 0)
                {
                    f_new[p * Q + opp[k]] += post;
                }
                else if (yn >= NY - 1)
                {
                    f_new[id(xn, NY - 2) * Q + k] += post;
                }
                else
                {
                    const int pn = id(xn, yn);

                    if (solid[pn])
                    {
                        f_new[p * Q + opp[k]] += post;
                    }
                    else
                    {
                        f_new[pn * Q + k] += post;
                    }
                }
            }
        }
    }

    f.swap(f_new);
}

void PoolBoilingLBM::applyBoundaryConditions()
{
    for (int x = 0; x < NX; ++x)
    {
        const int p_top = id(x, NY - 2);

        ux[p_top] *= 0.5;
        uy[p_top] = std::max(0.0, uy[p_top]);
    }
}

void PoolBoilingLBM::step()
{
    computeVelocity();
    evolveTemperature();
    evolvePhaseField();
    collideAndStream();
    applyBoundaryConditions();

    ++step_count;
}

void PoolBoilingLBM::writeFieldCSV(int step_id) const
{
    std::ostringstream filename;
    filename << "csv_out/fields_"
             << std::setw(7)
             << std::setfill('0')
             << step_id
             << ".csv";

    std::ofstream out(filename.str());

    out << "x,y,rho,rho_lbm,phi,T,dx,dy,speed,solid\n";
    out << std::setprecision(12);

    for (int y = 0; y < NY; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            const int p = id(x, y);

            double rho_lbm = 0.0;

            for (int k = 0; k < Q; ++k)
            {
                rho_lbm += f[p * Q + k];
            }

            const double phi_clamped = clampValue(phi[p], -1.0, 1.0);

            const double rho_phase =
                0.5 * (RHO_LIQUID_VIEW + RHO_VAPOR_VIEW) +
                0.5 * (RHO_LIQUID_VIEW - RHO_VAPOR_VIEW) * phi_clamped;

            const double dx = ux[p];
            const double dy = uy[p];
            const double speed = std::sqrt(dx * dx + dy * dy);

            out << x << ','
                << y << ','
                << rho_phase << ','
                << rho_lbm << ','
                << phi[p] << ','
                << temp[p] << ','
                << dx << ','
                << dy << ','
                << speed << ','
                << solid[p] << '\n';
        }
    }
}

void PoolBoilingLBM::writeDiagnosticsHeader() const
{
    std::ofstream out("diagnostics.csv");
    out << "step,vapor_cells,vapor_centroid_y,wall_vapor_cells,mean_temperature,max_speed\n";
}

void PoolBoilingLBM::appendDiagnostics(int step_id) const
{
    int vapor_cells = 0;
    int wall_vapor_cells = 0;
    int fluid_cells = 0;

    double y_sum = 0.0;
    double t_sum = 0.0;
    double max_speed = 0.0;

    for (int y = 1; y < NY - 1; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            const int p = id(x, y);

            ++fluid_cells;
            t_sum += temp[p];

            const double speed = std::sqrt(ux[p] * ux[p] + uy[p] * uy[p]);
            max_speed = std::max(max_speed, speed);

            if (phi[p] < 0.0)
            {
                ++vapor_cells;
                y_sum += static_cast<double>(y);

                if (y <= 3)
                {
                    ++wall_vapor_cells;
                }
            }
        }
    }

    const double vapor_centroid_y =
        vapor_cells > 0 ? y_sum / vapor_cells : 0.0;

    std::ofstream out("diagnostics.csv", std::ios::app);
    out << step_id << ','
        << vapor_cells << ','
        << vapor_centroid_y << ','
        << wall_vapor_cells << ','
        << t_sum / fluid_cells << ','
        << max_speed << '\n';
}