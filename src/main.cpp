#include "../inc/PoolBoilingLBM.hpp"
#include "../inc/main.hpp"

int main()
{
    PoolBoilingLBM sim;
    sim.initialize();
    sim.writeDiagnosticsHeader();

    constexpr int TOTAL_STEPS = 10000;
    constexpr int OUTPUT_EVERY = 200;

    // 输出初始时刻 0 步全场数据
    sim.writeFieldCSV(0);
    sim.appendDiagnostics(0);

    for (int step = 1; step <= TOTAL_STEPS; ++step)
    {
        sim.step();

        if (step % OUTPUT_EVERY == 0)
        {
            sim.writeFieldCSV(step);
            sim.appendDiagnostics(step);

            std::cout << "Output CSV at step = " << step << std::endl;
        }
    }

    std::cout << "Simulation finished." << std::endl;
    return 0;
}